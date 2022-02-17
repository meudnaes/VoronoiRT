using Unitful

import PhysicalConstants.CODATA2018: c_0, h, k_B, m_u, m_e, R_∞, ε_0, e
const E_∞ = R_∞ * c_0 * h
const hc = h * c_0

@derived_dimension PerLength Unitful.𝐋^-1
@derived_dimension PerArea Unitful.𝐋^-2
@derived_dimension NumberDensity Unitful.𝐋^-3
@derived_dimension ColumnDensity Unitful.𝐋^-2
@derived_dimension Volume Unitful.𝐋^3
@derived_dimension UnitsIntensity_λ Unitful.𝐋^-1 * Unitful.𝐌 * Unitful.𝐓^-3

"""
    Structure containing atmospheric grid and physical values at grid point
"""
struct Atmosphere
    z::Vector{typeof(1.0u"m")}
    x::Vector{typeof(1.0u"m")}
    y::Vector{typeof(1.0u"m")}
    temperature::Array{typeof(1.0u"K"), 3}
    electron_density::Array{typeof(1.0u"m^-3"), 3}
    hydrogen_populations::Array{typeof(1.0u"m^-3"), 3}
    velocity_z::Array{typeof(1.0u"m*s^-1"), 3}
    velocity_x::Array{typeof(1.0u"m*s^-1"), 3}
    velocity_y::Array{typeof(1.0u"m*s^-1"), 3}
    # Would be nice to include these things in the structure, but I couldn't get it to work
    # nz::Integer
    # nx::Integer
    # ny::Integer
    # Δx::Unitful.Length
    # Δy::Unitful.Length
    # function Atmosphere(z::Vector{<:Unitful.Length{T}},
    #                     x::Vector{<:Unitful.Length{T}},
    #                     y::Vector{<:Unitful.Length{T}}) where T <: AbstractFloat
    #                     # temperature::Array{<:Unitful.Temperature{T}, 3},
    #                     # electron_density::Array{<:NumberDensity{T}, 3},
    #                     # hydrogen_populations::Array{<:NumberDensity{T}, 3}) where T <: AbstractFloat
    #     #Funcy func
    #     nz = length(z)
    #     nx = length(x)
    #     ny = length(y)
    #     Δx = x[2] - x[1]
    #     Δy = y[2] - y[1]
    #     new{T}(z, x, y, nz, nx, ny, Δx, Δy)
    #            # temperature, electron_density, hydrogen_populations,
    #            # nz, nx, ny, Δx, Δy)
    # end
end

"""
    get_atmos()
Reads and slices atmosphere parameters from a bifrost atmosphere. Atmosphere
has to be stored in a hdf5 file. Returns atmosphere dimensions, velocity,
temperature, electron_density and hydrogen populations.

Original author: Ida Risnes Hansen
"""
function get_atmos(file_path; periodic=true, skip=1)
    println("---Extracting atmospheric data---")

    local x, y, z, temperature, electron_density, hydrogen_populations, velocity_z, velocity_x, velocity_y

    h5open(file_path, "r") do atmos
        z = read(atmos, "z")[1:skip:end]*u"m"
        x = read(atmos, "x")[1:skip:end]*u"m"
        y = read(atmos, "y")[1:skip:end]*u"m"

        velocity_z = read(atmos, "velocity_z")[1:skip:end, 1:skip:end, 1:skip:end]*u"m/s"
        velocity_x = read(atmos, "velocity_x")[1:skip:end, 1:skip:end, 1:skip:end]*u"m/s"
        velocity_y = read(atmos, "velocity_y")[1:skip:end, 1:skip:end, 1:skip:end]*u"m/s"

        temperature = read(atmos, "temperature")[1:skip:end, 1:skip:end, 1:skip:end]u"K"
        electron_density = read(atmos, "electron_density")[1:skip:end, 1:skip:end, 1:skip:end]u"m^-3"
        hydrogen_populations = read(atmos, "hydrogen_populations")[1:skip:end, 1:skip:end, 1:skip:end, 1, 1]u"m^-3"
    end

    #=
    if length(size(z)) == 2
        z = z[:,1]
        velocity_x = velocity_x[:,:,:,1]
        velocity_y = velocity_y[:,:,:,1]
        velocity_z = velocity_z[:,:,:,1]
        temperature = temperature[:,:,:,1]
        electron_density = electron_density[:,:,:,1]
        hydrogen_populations = hydrogen_populations[:,:,:,1,1]
    end
    =#

    if z[1] > z[end]
        reverse!(z)
        reverse!(velocity_x, dims=1)
        reverse!(velocity_y, dims=1)
        reverse!(velocity_z, dims=1)
        reverse!(temperature, dims=1)
        reverse!(electron_density, dims=1)
        reverse!(hydrogen_populations, dims=1)
    end

    if x[1] > x[end]
        reverse!(x)
        reverse!(velocity_x, dims=2)
        reverse!(velocity_y, dims=2)
        reverse!(velocity_z, dims=2)
        reverse!(temperature, dims=2)
        reverse!(electron_density, dims=2)
        reverse!(hydrogen_populations, dims=2)
    end

    if y[1] > y[end]
        reverse!(y)
        reverse!(velocity_x, dims=3)
        reverse!(velocity_y, dims=3)
        reverse!(velocity_z, dims=3)
        reverse!(temperature, dims=3)
        reverse!(electron_density, dims=3)
        reverse!(hydrogen_populations, dims=3)
    end

    if periodic
        println("---Periodic boundaries in x and y---")
        # Fix periodic boundaries with 'ghost' values

        # Fix grid
        x_periodic = periodic_borders(x)
        y_periodic = periodic_borders(y)

        # Add ghost layers on each side in x and y
        size_add = (0, 2, 2)

        # Temperature
        temperature_periodic = periodic_borders(temperature)

        # electron_density
        electron_density_periodic = periodic_borders(electron_density)

        # Hydrogen_populations
        hydrogen_populations_periodic = periodic_borders(hydrogen_populations)

        # velocity_z
        velocity_z_periodic = periodic_borders(velocity_z)

        # velocity_z
        velocity_x_periodic = periodic_borders(velocity_x)

        # velocity_z
        velocity_y_periodic = periodic_borders(velocity_y)

        return z, x_periodic, y_periodic, temperature_periodic, electron_density_periodic, hydrogen_populations_periodic, velocity_z_periodic, velocity_x_periodic, velocity_y_periodic
    end

    return z, x, y, temperature, electron_density, hydrogen_populations, velocity_z, velocity_x, velocity_y
end

"""
    periodic(vec::Vector{T}) where T<:Unitful.Quantity

Function to make boundaries of a vector correspond to periodic quantities in the
grid. Returns a new 'periodic' vector which has lenght 2+length(vec)
"""
function periodic_borders(vec::Vector{T}) where T<:Unitful.Quantity

    # Allocate space for periodic vector with same type as original vector
    periodic_vec = Vector{typeof(vec[1])}(undef, length(vec)+2)

    # Step size
    Δl = vec[2] - vec[1]

    # Fix boundaries
    periodic_vec[1] = vec[1] - Δl
    periodic_vec[end] = vec[end] + Δl

    # Fill in the 'body' of the vector
    periodic_vec[2:end-1] .= vec[:]

    return periodic_vec
end


"""
    periodic(arr::Array{T, 3}) where T<:Unitful.Quantity

Function to make boundaries of an array periodic. Returns a new periodic
array which has size (0,2,2)+size(arr)
"""
function periodic_borders(arr::Array{T, 3}) where T<:Unitful.Quantity

    # Allocate space for periodic vector with same type as original vector
    periodic_arr = Array{typeof(arr[1,1,1]), 3}(undef, size(arr).+(0,2,2))

    # Fix inner box
    periodic_arr[:, 2:end-1, 2:end-1] .= arr

    # x-direction
    periodic_arr[:,1,2:end-1] .= arr[:,end,:]
    periodic_arr[:,end,2:end-1] .= arr[:,1,:]

    # y-direction
    periodic_arr[:,2:end-1,end] .= arr[:,:,1]
    periodic_arr[:,2:end-1,1] .= arr[:,:,end]

    # fix corners
    periodic_arr[:,1,1] .= arr[:,end,end]
    periodic_arr[:,1,end] .= arr[:,end,1]
    periodic_arr[:,end,1] .= arr[:,1,end]
    periodic_arr[:,end,end] .= arr[:,1,1]

    return periodic_arr
end

function periodic_borders(arr::Array{T, 4}) where T<:Unitful.Quantity

    # Allocate space for periodic vector with same type as original vector
    periodic_arr = Array{typeof(arr[1,1,1,1]), 4}(undef, size(arr).+(0,0,2,2))

    # Fix inner box
    periodic_arr[:, :, 2:end-1, 2:end-1] .= arr

    # x-direction
    periodic_arr[:,:,1,2:end-1] .= arr[:,:,end,:]
    periodic_arr[:,:,end,2:end-1] .= arr[:,:,1,:]

    # y-direction
    periodic_arr[:,:,2:end-1,end] .= arr[:,:,:,1]
    periodic_arr[:,:,2:end-1,1] .= arr[:,:,:,end]

    # fix corners
    periodic_arr[:,:,1,1] .= arr[:,:,end,end]
    periodic_arr[:,:,1,end] .= arr[:,:,end,1]
    periodic_arr[:,:,end,1] .= arr[:,:,1,end]
    periodic_arr[:,:,end,end] .= arr[:,:,1,1]

    return periodic_arr
end

function periodic_pops(arr::Array{T, 4}) where T<:NumberDensity

    # Allocate space for periodic vector with same type as original vector
    periodic_arr = Array{typeof(arr[1,1,1,1]), 4}(undef, size(arr).+(0,0,2,2))

    # Fix inner box
    periodic_arr[:, 2:end-1, 2:end-1, :] .= arr

    # x-direction
    periodic_arr[:,1,2:end-1,:] .= arr[:,end,:,:]
    periodic_arr[:,end,2:end-1,:] .= arr[:,1,:,:]

    # y-direction
    periodic_arr[:,2:end-1,end,:] .= arr[:,:,1,:]
    periodic_arr[:,2:end-1,1,:] .= arr[:,:,end,:]

    # fix corners
    periodic_arr[:,1,1,:] .= arr[:,end,end,:]
    periodic_arr[:,1,end,:] .= arr[:,end,1,:]
    periodic_arr[:,end,1,:] .= arr[:,1,end,:]
    periodic_arr[:,end,end,:] .= arr[:,1,1,:]

    return periodic_arr
end
