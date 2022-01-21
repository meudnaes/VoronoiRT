"""
    Structure containing atmospheric grid and physical values at grid point
"""
struct Atmosphere
    z::Vector{<:Unitful.Length}
    x::Vector{<:Unitful.Length}
    y::Vector{<:Unitful.Length}
    temperature::Array{<:Unitful.Temperature, 3}
    electron_density::Array{<:NumberDensity, 3}
    hydrogen_populations::Array{<:NumberDensity, 3}
    velocity_z::Array{Unitful.Velocity, 3}
    velocity_x::Array{Unitful.Velocity, 3}
    velocity_y::Array{Unitful.Velocity, 3}
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
        Δx = x[2] - x[1]
        Δy = y[2] - y[1]

        # Fix grid
        x_periodic = Vector{Unitful.Length}(undef, length(x)+2)
        y_periodic = Vector{Unitful.Length}(undef, length(x)+2)

        x_periodic[2:end-1] = x
        x_periodic[1] = x[1] - Δx
        x_periodic[end] = x[end] + Δx

        y_periodic[2:end-1] = y
        y_periodic[1] = y[1] - Δy
        y_periodic[end] = y[end] + Δy

        # Add ghost layers on each side in x and y
        size_add = (0, 2, 2)

        # Temperature
        temperature_periodic = Array{Unitful.Temperature, 3}(undef, size(temperature) .+ size_add)
        temperature_periodic[:, 2:end-1, 2:end-1] = temperature
        # x-direction
        temperature_periodic[:,1,2:end-1] = temperature[:,end,:]
        temperature_periodic[:,end,2:end-1] = temperature[:,1,:]
        # y-direction
        temperature_periodic[:,2:end-1,end] = temperature[:,:,1]
        temperature_periodic[:,2:end-1,1] = temperature[:,:,end]
        # fix corners
        temperature_periodic[:,1,1] .= temperature[:,end,end]
        temperature_periodic[:,1,end] .= temperature[:,end,1]
        temperature_periodic[:,end,1] .= temperature[:,1,end]
        temperature_periodic[:,end,end] .= temperature[:,1,1]

        # Temperature
        electron_density_periodic = Array{NumberDensity, 3}(undef, size(electron_density) .+ size_add)
        electron_density_periodic[:, 2:end-1, 2:end-1] = electron_density
        # x-direction
        electron_density_periodic[:,1,2:end-1] = electron_density[:,end,:]
        electron_density_periodic[:,end,2:end-1] = electron_density[:,1,:]
        # y-direction
        electron_density_periodic[:,2:end-1,end] = electron_density[:,:,1]
        electron_density_periodic[:,2:end-1,1] = electron_density[:,:,end]
        # fix corners
        electron_density_periodic[:,1,1] .= electron_density[:,end,end]
        electron_density_periodic[:,1,end] .= electron_density[:,end,1]
        electron_density_periodic[:,end,1] .= electron_density[:,1,end]
        electron_density_periodic[:,end,end] .= electron_density[:,1,1]

        # Temperature
        hydrogen_populations_periodic = Array{NumberDensity, 3}(undef, size(hydrogen_populations) .+ size_add)
        hydrogen_populations_periodic[:, 2:end-1, 2:end-1] = hydrogen_populations
        # x-direction
        hydrogen_populations_periodic[:,1,2:end-1] = hydrogen_populations[:,end,:]
        hydrogen_populations_periodic[:,end,2:end-1] = hydrogen_populations[:,1,:]
        # y-direction
        hydrogen_populations_periodic[:,2:end-1,end] = hydrogen_populations[:,:,1]
        hydrogen_populations_periodic[:,2:end-1,1] = hydrogen_populations[:,:,end]
        # fix corners
        hydrogen_populations_periodic[:,1,1] .= hydrogen_populations[:,end,end]
        hydrogen_populations_periodic[:,1,end] .= hydrogen_populations[:,end,1]
        hydrogen_populations_periodic[:,end,1] .= hydrogen_populations[:,1,end]
        hydrogen_populations_periodic[:,end,end] .= hydrogen_populations[:,1,1]

        # velocity_z
        velocity_z_periodic = Array{Unitful.Velocity, 3}(undef, size(velocity_z) .+ size_add)
        velocity_z_periodic[:, 2:end-1, 2:end-1] = velocity_z
        # x-direction
        velocity_z_periodic[:,1,2:end-1] = velocity_z[:,end,:]
        velocity_z_periodic[:,end,2:end-1] = velocity_z[:,1,:]
        # y-direction
        velocity_z_periodic[:,2:end-1,end] = velocity_z[:,:,1]
        velocity_z_periodic[:,2:end-1,1] = velocity_z[:,:,end]
        # fix corners
        velocity_z_periodic[:,1,1] .= velocity_z[:,end,end]
        velocity_z_periodic[:,1,end] .= velocity_z[:,end,1]
        velocity_z_periodic[:,end,1] .= velocity_z[:,1,end]
        velocity_z_periodic[:,end,end] .= velocity_z[:,1,1]

        # velocity_z
        velocity_x_periodic = Array{Unitful.Velocity, 3}(undef, size(velocity_x) .+ size_add)
        velocity_x_periodic[:, 2:end-1, 2:end-1] = velocity_x
        # x-direction
        velocity_x_periodic[:,1,2:end-1] = velocity_x[:,end,:]
        velocity_x_periodic[:,end,2:end-1] = velocity_x[:,1,:]
        # y-direction
        velocity_x_periodic[:,2:end-1,end] = velocity_x[:,:,1]
        velocity_x_periodic[:,2:end-1,1] = velocity_x[:,:,end]
        # fix corners
        velocity_x_periodic[:,1,1] .= velocity_x[:,end,end]
        velocity_x_periodic[:,1,end] .= velocity_x[:,end,1]
        velocity_x_periodic[:,end,1] .= velocity_x[:,1,end]
        velocity_x_periodic[:,end,end] .= velocity_x[:,1,1]

        # velocity_z
        velocity_y_periodic = Array{Unitful.Velocity, 3}(undef, size(velocity_y) .+ size_add)
        velocity_y_periodic[:, 2:end-1, 2:end-1] = velocity_y
        # x-direction
        velocity_y_periodic[:,1,2:end-1] = velocity_y[:,end,:]
        velocity_y_periodic[:,end,2:end-1] = velocity_y[:,1,:]
        # y-direction
        velocity_y_periodic[:,2:end-1,end] = velocity_y[:,:,1]
        velocity_y_periodic[:,2:end-1,1] = velocity_y[:,:,end]
        # fix corners
        velocity_y_periodic[:,1,1] .= velocity_y[:,end,end]
        velocity_y_periodic[:,1,end] .= velocity_y[:,end,1]
        velocity_y_periodic[:,end,1] .= velocity_y[:,1,end]
        velocity_y_periodic[:,end,end] .= velocity_y[:,1,1]

        return z, x_periodic, y_periodic, temperature_periodic, electron_density_periodic, hydrogen_populations_periodic, velocity_z_periodic, velocity_x_periodic, velocity_y_periodic
    end

    return z, x, y, temperature, electron_density, hydrogen_populations, velocity_z, velocity_x, velocity_y
end
