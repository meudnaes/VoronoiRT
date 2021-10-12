using HDF5
using Unitful
using Distances
using Transparency
using LinearAlgebra
using NumericalIntegration

import PhysicalConstants.CODATA2018: h, c_0, k_B

@derived_dimension NumberDensity Unitful.𝐋^-3
@derived_dimension UnitsIntensity_λ Unitful.P * Unitful.L^-3

#=
    Structure containing atmospheric grid and physical values at grid point
=#
struct Atmosphere
    z::Vector{<:Unitful.Length}
    x::Vector{<:Unitful.Length}
    y::Vector{<:Unitful.Length}
    temperature::Array{<:Unitful.Temperature, 3}
    electron_density::Array{<:NumberDensity, 3}
    hydrogen_populations::Array{<:NumberDensity, 3}
end

#=
    function B_ν(ν, T)

Planck's law! Radiation in LTE. Takes frequency and temperature, returns
specific intensity
=#
function B_ν(ν, T)::AbstractFloat
    return 2*h*ν^3/c_0^2 * 1/(exp(h*ν/(k_B*T)) - 1)
end

#=
    function B_ν(λ, T)

Planck's law! Radiation in LTE. Takes wavelength and temperature, returns
specific intensity
=#
function B_λ(λ, T)
    return 2*h*c_0^2/λ^5 * 1/(exp(h*c_0/(λ*k_B*T)) - 1)
end

#=
    get_atmos()
Reads and slices atmosphere parameters from a bifrost atmosphere. Atmosphere
has to be stored in a hdf5 file. Returns atmosphere dimensions, velocity,
temperature, electron_density and hydrogen populations.

Original author: Ida Risnes Hansen
=#
function get_atmos(file_path; periodic=true)
    local x, y, z, hydrogen_populations
    h5open(file_path, "r") do atmos
        z = read(atmos, "z")u"m"
        x = read(atmos, "x")u"m"
        y = read(atmos, "y")u"m"

        velocity_x = read(atmos, "velocity_x")u"m/s"
        velocity_y = read(atmos, "velocity_y")u"m/s"
        velocity_z = read(atmos, "velocity_z")u"m/s"

        temperature = read(atmos, "temperature")u"K"
        electron_density = read(atmos, "electron_density")u"m^-3"
        hydrogen_populations = read(atmos, "hydrogen_populations")u"m^-3"
    end

    if length(size(z)) == 2
        z = z[:,1]
        velocity_x = velocity_x[:,:,:,1]
        velocity_y = velocity_y[:,:,:,1]
        velocity_z = velocity_z[:,:,:,1]
        temperature = temperature[:,:,:,1]
        electron_density = electron_density[:,:,:,1]
        hydrogen_populations = hydrogen_populations[:,:,:,1,1]
    end

    if z[1] > z[end]
        reverse!(z)
        reverse!(velocity_x, dims=1)
        reverse!(velocity_y, dims=1)
        reverse!(velocity_z, dims=1)
        reverse!(temperature, dims=1)
        reverse!(electron_density, dims=1)
        reverse!(hydrogen_populations, dims=1)
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

        # Fix quantities
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
        temperature_periodic[:,1,1] .= NaN*u"K"
        temperature_periodic[:,1,end] .= NaN*u"K"
        temperature_periodic[:,end,1] .= NaN*u"K"
        temperature_periodic[:,end,end] .= NaN*u"K"

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
        electron_density_periodic[:,1,1] .= NaN*u"m^-3"
        electron_density_periodic[:,1,end] .= NaN*u"m^-3"
        electron_density_periodic[:,end,1] .= NaN*u"m^-3"
        electron_density_periodic[:,end,end] .= NaN*u"m^-3"

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
        hydrogen_populations_periodic[:,1,1] .= NaN*u"m^-3"
        hydrogen_populations_periodic[:,1,end] .= NaN*u"m^-3"
        hydrogen_populations_periodic[:,end,1] .= NaN*u"m^-3"
        hydrogen_populations_periodic[:,end,end] .= NaN*u"m^-3"

        return z, x_periodic, y_periodic, temperature_periodic, electron_density_periodic, hydrogen_populations_periodic
    end

    return z, x, y, temperature, electron_density, hydrogen_populations
end

#=
    rejection_sampling(n_sites::Int)

Sample 3D vectors according to hydrogen distribution. Samples are obtained
through the acceptance-rejection method. Read
https://towardsdatascience.com/what-is-rejection-sampling-1f6aff92330d
for a nice explanation. Takes number of samples to generate. These are the
voronoi sites.

Make this method better by comparing with another distribution??? This will
increase the rate for acceptance.
Right now the reference distribution is the uniform pdf, scaled between N_H_max
and N_H_min.
=#
function rejection_sampling(n_sites::Int, atmos::Atmosphere)
    # Find max and min to convert random number between 0 and 1 to coordinate
    z_min = minimum(atmos.z); z_max = maximum(atmos.z)
    x_min = minimum(atmos.x); x_max = maximum(atmos.x)
    y_min = minimum(atmos.y); y_max = maximum(atmos.y)

    # Find max and min populations to scale uniform distribution
    N_H_min = minimum(atmos.hydrogen_populations)
    N_H_max = maximum(atmos.hydrogen_populations)

    # allocate arrays for new sites
    p_vec = Matrix{Unitful.Length}(undef, (3, n_sites))
    #z_new = Vector{Unitful.Length}(undef, n_sites)
    #x_new = Vector{Unitful.Length}(undef, n_sites)
    #y_new = Vector{Unitful.Length}(undef, n_sites)
    N_H_new = Vector{NumberDensity}(undef, n_sites)

    for i in 1:n_sites
        while true
            ref_vec = rand(Float64, 3)
            z_ref = ref_vec[1]*(z_max - z_min) + z_min
            x_ref = ref_vec[2]*(x_max - x_min) + x_min
            y_ref = ref_vec[3]*(y_max - y_min) + y_min

            density_ref = trilinear(z_ref, x_ref, y_ref, atmos, atmos.hydrogen_populations)
            density_ran = rand(Float32)*(N_H_max - N_H_min) + N_H_min
            if density_ran < density_ref
                # z_new[i] = z_ref
                # x_new[i] = x_ref
                # y_new[i] = y_ref
                p_vec[:, i] .= (z_ref, x_ref, y_ref)
                N_H_new[i] = trilinear(z_ref, x_ref, y_ref, atmos, atmos.hydrogen_populations)
                break
            end
        end
    end
    return p_vec, N_H_new #z_new, x_new, y_new, N_H_new
end

#=
    find_sites_sorted(z_new::Array, x_new::Array, y_new::Array,
                      z_bounds::Tuple{Float64, Float64},
                      x_bounds::Tuple{Float64, Float64},
                      y_bounds::Tuple{Float64, Float64})

Identifies and counts number of sites in (z_new, x_new, y_new) that are inside
the cube with bounds in (z, x, y) defined by z_bounds, x_bounds and y_bounds.
Assumes the coordinate z_new to be in increasing order.
=#
function find_sites_sorted(z_new::Array, x_new::Array, y_new::Array,
                           z_bounds::Tuple,
                           x_bounds::Tuple,
                           y_bounds::Tuple)

    hits = 0
    k_lower = searchsortedfirst(z_new, z_bounds[1])
    for k in k_lower:length(z_new)
        if z_new[k] > z_bounds[2]
            break
        elseif x_bounds[1] < x_new[k] < x_bounds[2] && y_bounds[1] < y_new[k] < y_bounds[2]
            hits = hits+1
        end
    end
    return hits::Int64
end

#=
    find_sites(z_new::Array, x_new::Array, y_new::Array,
               z_bounds::Tuple{Float64, Float64},
               x_bounds::Tuple{Float64, Float64},
               y_bounds::Tuple{Float64, Float64})

Same as `find_sites_sorted`, but doesn't assume any array to be sorted
=#
function find_sites(z_new::Array, x_new::Array, y_new::Array,
                    z_bounds::Tuple{Float64, Float64},
                    x_bounds::Tuple{Float64, Float64},
                    y_bounds::Tuple{Float64, Float64})

    hits = 0
    for k in 1:length(z_new)
        if z_bounds[1] < z_new[k] < z_bounds[2] && x_bounds[1] < x_new[k] < x_bounds[2] && y_bounds[1] < y_new[k] < y_bounds[2]
            hits = hits+1
        end
    end
    return hits::Int64
end

#=
    trilinear(x, y, z, hydrogen_populations)

Three-dimensional linear interpolation. Takes a function
$f: \mathbb{R}^3 -> \mathbb{R}$ and returns the trilinear interpolation in the
coordinates (x, y, z). Assumes an original cartesian grid (x, y, z), and
values for each grid point (hydrogen_populations) are defined
=#

function trilinear(z_mrk, x_mrk, y_mrk,
                   atmos::Atmosphere, vals::AbstractArray)
    # Returns the index of the first value in a greater than or equal to x
    # Subtract by 1 to get coordinate of lower corner
    idz = searchsortedfirst(atmos.z, z_mrk) - 1
    idx = searchsortedfirst(atmos.x, x_mrk) - 1
    idy = searchsortedfirst(atmos.y, y_mrk) - 1

    # bounding corner coordinates
    z0 = atmos.z[idz]; z1 = atmos.z[idz+1]
    x0 = atmos.x[idx]; x1 = atmos.x[idx+1]
    y0 = atmos.y[idy]; y1 = atmos.y[idy+1]

    # difference between coordinates and interpolation point
    x_d = (x_mrk - x0)/(x1 - x0)
    y_d = (y_mrk - y0)/(y1 - y0)
    z_d = (z_mrk - z0)/(z1 - z0)

    # values at each corner (z is first index in data array)
    c000 = vals[idz, idx, idy]
    c010 = vals[idz, idx, idy+1]
    c100 = vals[idz, idx+1, idy]
    c110 = vals[idz, idx+1, idy+1]
    c001 = vals[idz+1, idx, idy]
    c011 = vals[idz+1, idx, idy+1]
    c101 = vals[idz+1, idx+1, idy]
    c111 = vals[idz+1, idx+1, idy+1]

    # interpolate in x direction
    c00 = c000*(1 - x_d) + c100*x_d
    c01 = c001*(1 - x_d) + c101*x_d
    c10 = c010*(1 - x_d) + c100*x_d
    c11 = c011*(1 - x_d) + c111*x_d

    # interpolate in y direction
    c0 = c00*(1 - y_d) + c10*y_d
    c1 = c01*(1 - y_d) + c11*y_d

    # intepolate in z direction
    c = c0*(1 - z_d) + c1*z_d
    return c
end

function trilinear_no_search(z_mrk, x_mrk, y_mrk, idz, idx, idy,
                   atmos::Atmosphere, vals::AbstractArray)

    # This function fulfills what I wanted it to do, from the comments in
    # function bilinear(...). It's periodic in x and y. NB it should never be
    # periodic in z, so it isn't, as that wouldn't make sense at all. Used
    # internally for the short characteristics calculations.

    nx = length(atmos.x)
    ny = length(atmos.y)

    # bounding corner coordinates
    z0 = atmos.z[idz]; z1 = atmos.z[idz+1]
    x0 = atmos.x[idx]; x1 = atmos.x[idx%nx+1]
    y0 = atmos.y[idy]; y1 = atmos.y[idy%ny+1]

    # difference between coordinates and interpolation point
    x_d = (x_mrk - x0)/(x1 - x0)
    y_d = (y_mrk - y0)/(y1 - y0)
    z_d = (z_mrk - z0)/(z1 - z0)

    # values at each corner (z is first index in data array)
    c000 = vals[idz, idx, idy]
    c010 = vals[idz, idx, idy%ny+1]
    c100 = vals[idz, idx%nx+1, idy]
    c110 = vals[idz, idx%nx+1, idy%ny+1]
    c001 = vals[idz+1, idx, idy]
    c011 = vals[idz+1, idx, idy%ny+1]
    c101 = vals[idz+1, idx%nx+1, idy]
    c111 = vals[idz+1, idx%nx+1, idy%ny+1]

    # interpolate in x direction
    c00 = c000*(1 - x_d) + c100*x_d
    c01 = c001*(1 - x_d) + c101*x_d
    c10 = c010*(1 - x_d) + c100*x_d
    c11 = c011*(1 - x_d) + c111*x_d

    # interpolate in y direction
    c0 = c00*(1 - y_d) + c10*y_d
    c1 = c01*(1 - y_d) + c11*y_d

    # intepolate in z direction
    c = c0*(1 - z_d) + c1*z_d
    return c
end

function bilinear(x_mrk, y_mrk, x_bounds, y_bounds, vals)

    # corner coordinates
    x1 = x_bounds[1]; x2 = x_bounds[2]
    y1 = y_bounds[1]; y2 = y_bounds[2]

    # function value for each corner. 11 is lower left corner, 22 is upper right
    Q11 = vals[1,1]     # (x1, y1)
    Q12 = vals[1,2]     # (x1, y2)
    Q21 = vals[2,1]     # (x2, y1)
    Q22 = vals[2,2]     # (x2, y2)

    dx = x2 - x1
    dy = y2 - y1

    f1 = ((x2 - x_mrk)*Q11 + (x_mrk - x1)*Q21)/dx
    f2 = ((x2 - x_mrk)*Q12 + (x_mrk - x1)*Q22)/dx

    f = ((y2 - y_mrk)*f1 + (y_mrk - y1)*f2)/dy

    return f
end

#=
    avg_mass(k::Int64, i::Int64, j::Int64)

Finds the mass of a cell by taking the average mass density over all cell
corners and multiplies with the volume of the cell. Assumes an original
cartesian grid (x, y, z) and a corresponding value (hydrogen_populations) for
each grid point
=#
function mass_function(k::Int64, i::Int64, j::Int64, atmos::Atmosphere)
    volume = Float32((atmos.z[k+1] - atmos.z[k])
                    *(atmos.x[i+1] - atmos.x[i])
                    *(atmos.y[j+1] - atmos.y[j]))
    avg_density = (atmos.hydrogen_populations[k, i, j] +
                   atmos.hydrogen_populations[k, i, j+1] +
                   atmos.hydrogen_populations[k, i+1, j+1] +
                   atmos.hydrogen_populations[k, i+1, j] +
                   atmos.hydrogen_populations[k+1, i, j] +
                   atmos.hydrogen_populations[k+1, i, j+1] +
                   atmos.hydrogen_populations[k+1, i+1, j+1] +
                   atmos.hydrogen_populations[k+1, i+1, j])/8
    mass = volume*avg_density
    return mass::Float32
end

function write_arrays(z::AbstractArray, x::AbstractArray, y::AbstractArray,
                      fname::String)

    if length(z) != length(y) || length(y) != length(x)
        println("Wrong shapes of input data")
        exit()
    end

    open(fname, "w") do io
        for i in 1:length(z)
            println(io, "$i\t$(x[i])\t$(y[i])\t$(z[i])")
        end
    end
end

# Physics functions
function αcont(λ::Unitful.Length, temperature::Unitful.Temperature,
               electron_density::NumberDensity, h_ground_density::NumberDensity,
               proton_density::NumberDensity)
    α = Transparency.hminus_ff_stilley(λ, temperature, h_ground_density, electron_density)
    # Wait with this
    #α = Transparency.hminus_bf_geltman(λ, temperature, h_ground_density, electron_density)
    α += hydrogenic_ff(c_0 / λ, temperature, electron_density, proton_density, 1)
    α += h2plus_ff(λ, temperature, h_ground_density, proton_density)
    α += h2plus_bf(λ, temperature, h_ground_density, proton_density)
    α += thomson(electron_density)
    α += rayleigh_h(λ, h_ground_density)
    return α
end

function α_scattering(λ::Unitful.Length, temperature::Unitful.Temperature,
               electron_density::NumberDensity, h_ground_density::NumberDensity,
               proton_density::NumberDensity)
   α = thomson(electron_density)
   α += rayleigh_h(λ, h_ground_density)
   return α
end

function J_ν(weights, intensities)
    # Gaussian quadrature to calculate mean intensity
    J = 0u"kW*m^-2*nm^-1"
    for (weight, intensity) in zip(weights, intensitites)
        J = J + weight*intensity
    end
    return J::Float64
end

#=
    read_neighbours(fname::String, n_sites::Int)

Reads a file containing neighbouring cells for each grid point in the voronoi
tesselation.
=#
function read_neighbours(fname::String, n_sites::Int)::AbstractMatrix
    ID = Vector{Int64}(undef, n_sites)
    neighbours = zeros(Int64, n_sites, 40)
    open(fname, "r") do io
        for (i, l) in enumerate(eachline(io))
            ID[i] = parse(Int64, split(l)[1])
            for j in 2:length(split(l))
                neighbours[i, j-1] = parse(Int64, split(l)[j])
            end
        end
    end
    return neighbours[sortperm(ID), :]
end

function trapezoidal(Δx, a, b)
    area = Δx*(a + b)/2
    return area
end

function xy_intersect(ϕ::AbstractFloat)
    local sign_x, sign_y
    if ϕ < π/2
        # 1st quadrant. Positive x, positive y
        print("1st, ")
        sign_x = 1
        sign_y = 1
    elseif π/2 < ϕ < π
        # 2nd quadrant. Negative x, positive y
        print("2nd, ")
        sign_x = -1
        sign_y = 1
    elseif π < ϕ < 3π/2
        # 3rd quadrant. Negative x, negative y
        print("3rd, ")
        sign_x = -1
        sign_y = -1
    elseif 3π/2 < ϕ < 2π
        # 4th quadrant. Positive x, negative y
        print("4th, ")
        sign_x = 1
        sign_y = -1
    end
    return sign_x::Int, sign_y::Int
end

function z_intersect(θ)
    if θ < π/2
        sign_z = 0
    elseif θ > π/2
        sign_z = -1
    end
    return sign_z
end

function read_quadrature(fname::String)
    # elaborate (bad) scheme to extract number of points from filename
    n_points = ""
    switch = false
    for (i, char) in enumerate(fname)
        if char == 'n'
            switch = true
        elseif switch == true
            try parse(Int, char)
                n_points = string(n_points, char)
            catch
                break
            end
        end
    end

    n_points = parse(Int, n_points)

    weights = zeros(n_points)
    θ_array = zeros(n_points)
    ϕ_array = zeros(n_points)

    open(fname, "r") do io
        for (i, line) in enumerate(eachline(fname))
            weights[i] = parse(Float64, split(line)[1])
            θ_array[i] = parse(Float64, split(line)[2])
            ϕ_array[i] = parse(Float64, split(line)[3])
        end
    end

    return weights::AbstractArray, θ_array::AbstractArray, ϕ_array::AbstractArray, n_points::Int
end

function range_bounds(sign::Int, bound::Int)
    if sign == -1
        start = 2
        stop = bound-1
    elseif sign == 1
        start = bound-1
        stop = 2
    end
    return start::Int, stop::Int
end
