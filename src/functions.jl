using HDF5
using Random
using Unitful
using Distances
using Transparency
using LinearAlgebra

import PhysicalConstants.CODATA2018: h, c_0, k_B, m_p

@derived_dimension PerLength Unitful.𝐋^-1
@derived_dimension PerArea Unitful.𝐋^-2
@derived_dimension NumberDensity Unitful.𝐋^-3
@derived_dimension ColumnDensity Unitful.𝐋^-2
@derived_dimension Volume Unitful.𝐋^3
@derived_dimension UnitsIntensity_λ Unitful.P * Unitful.L^-3

"""
    function B_ν(ν, T)

Planck's law! Radiation in LTE. Takes frequency and temperature, returns
specific intensity
"""
function B_ν(ν, T)::AbstractFloat
    return 2*h*ν^3/c_0^2 * 1/(exp(h*ν/(k_B*T)) - 1)
end

"""
    function B_ν(λ, T)

Planck's law! Radiation in LTE. Takes wavelength and temperature, returns
specific intensity
"""
function B_λ(λ, T)
    return 2*h*c_0^2/λ^5 * 1/(exp(h*c_0/(λ*k_B*T)) - 1)
end

"""
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
"""
function rejection_sampling(n_sites::Int, atmos::Atmosphere, quantity::AbstractArray)
    # Find max and min to convert random number between 0 and 1 to coordinate
    println("---Sampling new sites---")
    z_min = minimum(atmos.z); z_max = maximum(atmos.z)
    x_min = minimum(atmos.x); x_max = maximum(atmos.x)
    y_min = minimum(atmos.y); y_max = maximum(atmos.y)

    Δz = z_max - z_min
    Δx = x_max - x_min
    Δy = y_max - y_min

    # Find max and min populations to scale uniform distribution
    q_min = minimum(quantity)
    q_max = maximum(quantity)

    Δq = q_max - q_min

    # allocate arrays for new sites
    p_vec = Matrix{Unitful.Length}(undef, (3, n_sites))

    for i in 1:n_sites
        print("site $i/$n_sites \r")
        while true
            ref_vec = rand(Float64, 3)
            z_ref = ref_vec[1]*Δz + z_min
            x_ref = ref_vec[2]*Δx + x_min
            y_ref = ref_vec[3]*Δy + y_min

            # acceptance criterion, "reference"
            density_ref = trilinear(z_ref, x_ref, y_ref, atmos, quantity)
            # random sample, compare to reference
            density_ran = rand(Float64)*Δq + q_min
            if density_ref > density_ran
                # a point is accepted, store position and move on
                p_vec[:, i] .= (z_ref, x_ref, y_ref)
                # break to find next site
                break
            end
        end
    end
    print("                                                                 \r")
    return p_vec
end

function rejection_sampling(n_sites::Int, boundaries::Matrix, quantity::AbstractArray)
    # Find max and min to convert random number between 0 and 1 to coordinate
    println("---Sampling new sites---")
    z_min = boundaries[1,1]; z_max = boundaries[1,2]
    x_min = boundaries[2,1]; x_max = boundaries[2,2]
    y_min = boundaries[3,1]; y_max = boundaries[3,2]

    Δz = z_max - z_min
    Δx = x_max - x_min
    Δy = y_max - y_min

    nz = length(quantity[:,1,1])
    nx = length(quantity[1,:,1])
    ny = length(quantity[1,1,:])

    z = Vector{Float64}(undef, nz)
    x = Vector{Float64}(undef, nx)
    y = Vector{Float64}(undef, ny)

    δz = (z_max - z_min)/(nz + 1)
    δx = (x_max - x_min)/(nx + 1)
    δy = (y_max - y_min)/(ny + 1)

    z[1] = z_min
    x[1] = x_min
    y[1] = y_min

    for i in 2:nz
        z[i] = z[i-1] + δz
    end

    for i in 2:nx
        x[i] = x[i-1] + δx
    end

    for i in 2:ny
        y[i] = y[i-1] + δy
    end

    z[end] = z_max
    x[end] = x_max
    y[end] = y_max

    # Find max and min populations to scale uniform distribution
    q_min = minimum(quantity)
    q_max = maximum(quantity)

    Δq = q_max - q_min

    # allocate arrays for new sites
    p_vec = Matrix{Unitful.Length}(undef, (3, n_sites))

    for i in 1:n_sites
        print("site $i/$n_sites \r")
        while true
            ref_vec = rand(Float64, 3)
            z_ref = ref_vec[1]*Δz + z_min
            x_ref = ref_vec[2]*Δx + x_min
            y_ref = ref_vec[3]*Δy + y_min

            # acceptance criterion, "reference"
            density_ref = trilinear(z_ref, x_ref, y_ref, z, x, y, quantity)
            # random sample, compare to reference
            density_ran = rand(Float64)*Δq + q_min
            if density_ref > density_ran
                # a point is accepted, store position and move on
                p_vec[:, i] .= (z_ref, x_ref, y_ref)
                # break to find next site
                break
            end
        end
    end
    print("                                                                 \r")
    return p_vec
end


"""
    find_sites_sorted(z_new::Array, x_new::Array, y_new::Array,
                      z_bounds::Tuple{Float64, Float64},
                      x_bounds::Tuple{Float64, Float64},
                      y_bounds::Tuple{Float64, Float64})

Identifies and counts number of sites in (z_new, x_new, y_new) that are inside
the cube with bounds in (z, x, y) defined by z_bounds, x_bounds and y_bounds.
Assumes the coordinate z_new to be in increasing order.
"""
function find_sites_sorted(p_new::AbstractArray,
                           z_bounds::Tuple,
                           x_bounds::Tuple,
                           y_bounds::Tuple)

    z_new = p_new[1, :]
    x_new = p_new[2, :]
    y_new = p_new[3, :]

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

"""
    find_sites(z_new::Array, x_new::Array, y_new::Array,
               z_bounds::Tuple{Float64, Float64},
               x_bounds::Tuple{Float64, Float64},
               y_bounds::Tuple{Float64, Float64})

Same as `find_sites_sorted`, but doesn't assume any array to be sorted
"""
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

@doc raw"""
    trilinear(x_mrk, y_mrk, z_mrk, hydrogen_populations)

Three-dimensional linear interpolation. Takes a function
$f: \mathbb{R}^3 -> \mathbb{R}$ and returns the trilinear interpolation in the
coordinates (x, y, z). Assumes an original cartesian grid (x, y, z), and
values for each grid point (hydrogen_populations) are defined
"""
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

function trilinear(z_mrk::AbstractFloat, x_mrk::AbstractFloat, y_mrk::AbstractFloat,
                   x::AbstractVector, y::AbstractVector, z::AbstractVector, vals::AbstractArray)
    # Returns the index of the first value in a greater than or equal to x
    # Subtract by 1 to get coordinate of lower corner
    idz = searchsortedfirst(z, z_mrk) - 1
    idx = searchsortedfirst(x, x_mrk) - 1
    idy = searchsortedfirst(y, y_mrk) - 1

    # bounding corner coordinates
    z0 = z[idz]; z1 = z[idz+1]
    x0 = x[idx]; x1 = x[idx+1]
    y0 = y[idy]; y1 = y[idy+1]

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

@doc raw"""
    bilinear(x_mrk, y_mrk, x_bounds, y_bounds, vals)

Two-dimensional linear interpolation. Takes a function
$f: \mathbb{R}^2 -> \mathbb{R}$ and returns the biilinear interpolation in the
coordinates (x_mrk, y_mrk). Assumes an underlying cartesian grid (x_i, y_i)
Values are defined on the corners of the rectangle spanned by x_bounds and
y_bounds. x_mrk and y_mrk have to lie inside this rectangle.
"""
function bilinear(x_mrk, y_mrk, x_bounds, y_bounds, vals)

    # corner coordinates
    x1 = x_bounds[1]; x2 = x_bounds[2]
    y1 = y_bounds[1]; y2 = y_bounds[2]

    # function value for each corner. 11 is lower left corner, 22 is upper right
    Q11 = vals[1,1]     # (x1, y1)
    Q12 = vals[1,2]     # (x1, y2)
    Q21 = vals[2,1]     # (x2, y1)
    Q22 = vals[2,2]     # (x2, y2)

    # Rectangle side length
    dx = x2 - x1
    dy = y2 - y1

    # Interpolate in x-direction
    f1 = ((x2 - x_mrk)*Q11 + (x_mrk - x1)*Q21)/dx
    f2 = ((x2 - x_mrk)*Q12 + (x_mrk - x1)*Q22)/dx

    # Interpolate in y-direction
    f = ((y2 - y_mrk)*f1 + (y_mrk - y1)*f2)/dy
    return f
end

"""
    mass_function(k::Int64, i::Int64, j::Int64)

Finds the mass of a cell by taking the average mass density over all cell
corners and multiplies with the volume of the cell. Assumes an original
cartesian grid (x, y, z) and a corresponding value (hydrogen_populations) for
each grid point
"""
function mass_function(k::Int64, i::Int64, j::Int64, atmos::Atmosphere)
    volume = (atmos.z[k+1] - atmos.z[k])*(atmos.x[i+1] - atmos.x[i])*(atmos.y[j+1] - atmos.y[j])
    avg_density = (atmos.hydrogen_populations[k, i, j] +
                   atmos.hydrogen_populations[k, i, j+1] +
                   atmos.hydrogen_populations[k, i+1, j+1] +
                   atmos.hydrogen_populations[k, i+1, j] +
                   atmos.hydrogen_populations[k+1, i, j] +
                   atmos.hydrogen_populations[k+1, i, j+1] +
                   atmos.hydrogen_populations[k+1, i+1, j+1] +
                   atmos.hydrogen_populations[k+1, i+1, j])/8*m_p
    mass = volume*avg_density
    return mass
end

"""
    function write_arrays(x::AbstractArray, y::AbstractArray, z::AbstractArray,
                          fname::String)

Writes the arrays z, x, and y to a file with filename fname.
Arrays are written in columns [ row number ] [ x ] [ y ] [ z ]
"""
function write_arrays(x::AbstractArray, y::AbstractArray, z::AbstractArray,
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

function write_boundaries(z_min, z_max, x_min, x_max, y_min, y_max, fname::String)
    open(fname, "w") do io
        println(io, "z_min = $z_min")
        println(io, "z_max = $z_max")
        println(io, "x_min = $x_min")
        println(io, "x_max = $x_max")
        println(io, "y_min = $y_min")
        println(io, "y_max = $y_max")
    end
end

"""
    function α_cont(λ::Unitful.Length, temperature::Unitful.Temperature,
               electron_density::NumberDensity, h_ground_density::NumberDensity,
               proton_density::NumberDensity)

Continuum extinction.
"""
function α_continuum(λ::Unitful.Length, temperature::Unitful.Temperature,
               electron_density::NumberDensity, h_ground_density::NumberDensity,
               proton_density::NumberDensity)

    α = max(0u"m^-1", Transparency.hminus_ff_stilley(λ, temperature, h_ground_density, electron_density))
    α += Transparency.hminus_bf_wbr(λ, temperature, h_ground_density, electron_density)
    α += hydrogenic_ff(c_0 / λ, temperature, electron_density, proton_density, 1)
    α += h2plus_ff(λ, temperature, h_ground_density, proton_density)
    α += h2plus_bf(λ, temperature, h_ground_density, proton_density)
    α += thomson(electron_density)
    α += rayleigh_h(λ, h_ground_density)
    return α
end

"""
    function α_scattering(λ::Unitful.Length, temperature::Unitful.Temperature,
               electron_density::NumberDensity, h_ground_density::NumberDensity,
               proton_density::NumberDensity)

Extinction from scattering processes.
"""
function α_scattering(λ::Unitful.Length, temperature::Unitful.Temperature,
               electron_density::NumberDensity, h_ground_density::NumberDensity,
               proton_density::NumberDensity)

   α = thomson(electron_density)
   α += rayleigh_h(λ, h_ground_density)
   return α
end

"""
    function α_absorption(λ::Unitful.Length, temperature::Unitful.Temperature,
               electron_density::NumberDensity, h_ground_density::NumberDensity,
               proton_density::NumberDensity)

Extinction from photon destruction processes.
"""
function α_absorption(λ::Unitful.Length, temperature::Unitful.Temperature,
               electron_density::NumberDensity, h_ground_density::NumberDensity,
               proton_density::NumberDensity)

    α = max(0u"m^-1", Transparency.hminus_ff_stilley(λ, temperature, h_ground_density, electron_density))
    α += Transparency.hminus_bf_wbr(λ, temperature, h_ground_density, electron_density)
    α += hydrogenic_ff(c_0 / λ, temperature, electron_density, proton_density, 1)
    α += h2plus_ff(λ, temperature, h_ground_density, proton_density)
    α += h2plus_bf(λ, temperature, h_ground_density, proton_density)
end

"""
    trapezoidal(Δx::AbstractFloat, a::AbstractFloat, b::AbstractFloat)

Trapezoidal rule.
"""
function trapezoidal(Δx, a, b)
    area = Δx*(a + b)/2
    return area
end

"""
    xy_intersect(ϕ::AbstractFloat)

Finds quadrant defined by the azimuthal angle ϕ
"""
function xy_intersect(ϕ::AbstractFloat)
    local sign_x, sign_y
    if ϕ < π/2
        # 1st quadrant. Positive x, positive y
        sign_x = 1
        sign_y = 1
    elseif π/2 < ϕ < π
        # 2nd quadrant. Negative x, positive y
        sign_x = -1
        sign_y = 1
    elseif π < ϕ < 3π/2
        # 3rd quadrant. Negative x, negative y
        sign_x = -1
        sign_y = -1
    elseif 3π/2 < ϕ < 2π
        # 4th quadrant. Positive x, negative y
        sign_x = 1
        sign_y = -1
    end
    return sign_x::Int, sign_y::Int
end

"""
    function read_quadrature(fname::String)

Read quarature weights and angles from file. Returns weights, horizontal angle,
azimuthal angle, and number of quadrature points. Quadratures found in
https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/645/A101#/browse
from Bestard & Bueno (2021)
"""
function read_quadrature(fname::String)
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

"""
    function range_bounds(sign::Int, bound::Int)

Given a quadrant from sign_x and sign_x from xy_intersect(), this function
determines the loop start end stop point for the short characteristics ray.
"""
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

"""
Computes weights for linear integration of source function,
approximating `exp(-Δτ)` for very small and very large values of `Δτ`.
"""
function weights2(Δτ::T) where T <: AbstractFloat
    if Δτ < 5e-4
        w1 = Δτ * (1 - Δτ / 2)
        w2 = Δτ^2 * (0.5f0 - Δτ / 3)
    elseif Δτ > 50
        w1 = w2 = one(T)
    else
        expΔτ = exp(-Δτ)
        w1 = 1 - expΔτ
        w2 = w1 - Δτ * expΔτ
    end
    return w1, w2
end

function weights(Δτ)
    e0 = 1 - exp(-Δτ)
    e1 = Δτ - e0
    return e0, e1
end


"""
    function coefficients(e0, e1, Δτ_upwind)

Coefficients for integrating the formal solution with a linear interpolation of
the source function. e0 and e1 are the weights.
"""
function coefficients(e0, e1, Δτ_upwind)
    if Δτ_upwind == 0
        a = 0
        b = 0
        c = 1
    else
        a = e0 - e1/Δτ_upwind
        b = e1/Δτ_upwind
        c = exp(-Δτ_upwind)
    end
    return a, b, c
end

function smallest_non_negative(arr::AbstractArray)
    minVal = Inf
    index = 0
    for i in 1:length(arr)
        if 0 < arr[i] < minVal
            minVal = arr[i]
            index = i
        end
    end
    return index::Int, minVal::Float64
end

"""
    function circle_shape(x, y, r)

Function for a circle with centre coordinates (x, y) and radius r. Useful for
plotting purposes.
"""
function circle_shape(x, y, r)
    θ = LinRange(0, 2π, 500)
    x .+ r*cos.(θ), y .+ r*sin.(θ)
end

function less_than(arr::AbstractArray, treshold)
    indices = Vector{Int}(undef, length(arr))
    j = 0
    for i in 1:length(arr)
        if arr[i] < treshold
            j += 1
            indices[j] = i
        end
    end
    return indices[1:j]
end

function arg_where(arr, num)
    indices = Vector{Int}(undef, length(arr))
    j = 0
    for i in eachindex(arr)
        if arr[i] == num
            j += 1
            indices[j] = i
        end
    end
    return indices[1:j]
end

function line_of_sight_velocity(atmos::Atmosphere, k::Vector)
    v_los = Array{Unitful.Velocity, 3}(undef, size(atmos.velocity_z))

    for kk in 1:length(atmos.z)
        for ii in 1:length(atmos.x)
            for jj in 1:length(atmos.y)
                velocity = [atmos.velocity_z[kk, ii, jj],
                            atmos.velocity_x[kk, ii, jj],
                            atmos.velocity_y[kk, ii, jj]]

                v_los[kk, ii, jj] = dot(velocity, k)
            end
        end
    end
    return v_los::Array{Unitful.Velocity, 3}
end
