using HDF5
using Random
using Unitful

include("io.jl")

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
    p_vec = Matrix{Float64}(undef, (3, n_sites))u"m"

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
    p_vec = Matrix{Float64}(undef, (3, n_sites))

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
    return p_vec*1u"m"
end

@doc raw"""
    trilinear(x_mrk, y_mrk, z_mrk, hydrogen_populations)

Three-dimensional linear interpolation. Takes a function
$f: \mathbb{R}^3 -> \mathbb{R}$ and returns the trilinear interpolation in the
coordinates (x, y, z). Assumes an original cartesian grid (x, y, z), and
values for each grid point (hydrogen_populations) are defined
"""
function trilinear(z_mrk::Unitful.Length, x_mrk::Unitful.Length, y_mrk::Unitful.Length,
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
    c10 = c010*(1 - x_d) + c110*x_d
    c11 = c011*(1 - x_d) + c111*x_d

    # interpolate in y direction
    c0 = c00*(1 - y_d) + c10*y_d
    c1 = c01*(1 - y_d) + c11*y_d

    # intepolate in z direction
    c = c0*(1 - z_d) + c1*z_d
    return c
end

function trilinear(z_mrk::Unitful.Length, x_mrk::Unitful.Length, y_mrk::Unitful.Length,
                   z::Vector{<:Unitful.Length}, x::Vector{<:Unitful.Length}, y::Vector{<:Unitful.Length},
                   vals::AbstractArray)
    # Returns the index of the first value in a greater than or equal to x
    # Subtract by 1 to get coordinate of lower corner
    idz = searchsortedfirst(z, z_mrk) - 1
    idx = searchsortedfirst(x, x_mrk) - 1
    idy = searchsortedfirst(y, y_mrk) - 1

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
    c10 = c010*(1 - x_d) + c110*x_d
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
function bilinear(x_mrk::Unitful.Length, y_mrk::Unitful.Length,
                  x_bounds::Tuple, y_bounds::Tuple,
                  vals::AbstractMatrix)

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

function bilinear(x_mrk::Float64, y_mrk::Float64,
                  x_bounds::Tuple, y_bounds::Tuple,
                  vals::AbstractMatrix)

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
    trapezoidal(Δx::AbstractFloat, a::AbstractFloat, b::AbstractFloat)

Trapezoidal rule.
"""
function trapezoidal(Δx, a, b)
    area = Δx*(a + b)/2
    return area
end

"""
    xy_intersect(ϕ::AbstractFloat)

Finds x and y direction ray is moving in defined by the azimuthal angle ϕ
"""
function xy_intersect(ϕ::Float64)
    local sign_x, sign_y
    @assert 0 <= ϕ < 2π "Bad angle. Expects ϕ ∈ [0, 2π), but ϕ = $ϕ"
    if 0 <= ϕ < π/2
        # 1st quadrant. Negative x, negative y
        sign_x = -1
        sign_y = -1
    elseif π/2 < ϕ < π
        # 2nd quadrant. Positive x, negative y
        sign_x = 1
        sign_y = -1
    elseif π < ϕ < 3π/2
        # 3rd quadrant. Positive x, positive y
        sign_x = 1
        sign_y = 1
    elseif 3π/2 < ϕ < 2π
        # 4th quadrant. Negative x, positive y
        sign_x = -1
        sign_y = 1
    end
    return sign_x::Int, sign_y::Int
end

"""
    xy_intersect(k::Vector{Float64})

Finds x and y direction ray is moving in defined by direction vector k
"""
function xy_intersect(k::Vector{Float64})
    local sign_x, sign_y
    @assert norm(k) ≈ 1 "Bad direction vector. Expects length to be 1"
    if k[2] > 0 && k[3] > 0
        # 1st quadrant. Negative x, negative y
        sign_x = -1
        sign_y = -1
    elseif k[2] < 0 && k[3] > 0
        # 2nd quadrant. Positive x, negative y
        sign_x = 1
        sign_y = -1
    elseif k[2] < 0 && k[3] < 0
        # 3rd quadrant. Positive x, positive y
        sign_x = 1
        sign_y = 1
    elseif k[2] > 0 && k[3] < 0
        # 4th quadrant. Negative x, positive y
        sign_x = -1
        sign_y = 1
    else
        # theta is 0 or 180, ray either straight up or down, doesn't really matter what
        # the sign is here, interpolation is exact anyways as the ray hits the
        # grid points...
        sign_x = 1
        sign_y = 1
    end
    return sign_x::Int, sign_y::Int
end


"""
    function range_bounds(sign::Int, bound::Int)

Given a quadrant from sign_x and sign_x from xy_intersect(), this function
determines the loop start and stop point for the short characteristics rays.
"""
function range_bounds(sign::Int, bound::Int)
    if sign == 1
        start = 2
        stop = bound-1
    elseif sign == -1
        start = bound-1
        stop = 2
    end
    return start::Int, stop::Int
end

"""
Computes weights for linear integration of source function,
approximating `exp(-Δτ)` for very small and very large values of `Δτ`.
"""
function weights(Δτ::T) where T <: AbstractFloat
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

"""
    coefficients(w1, w2, Δτ_upwind)

Coefficients for integrating the formal solution with a linear interpolation of
the source function. w1 and w2 are the weights.
"""
function coefficients(w1::Float64, w2::Float64, Δτ_upwind::Float64)
    if Δτ_upwind == 0
        a = 0
        b = 0
        c = 1
    else
        a = w2/Δτ_upwind
        b = w1 - w2/Δτ_upwind
        c = exp(-Δτ_upwind)
    end
    return a, b, c
end

"""
    linear_weights(Δτ::AbstractFloat)

Compute weights for linear integration of Source function in formal solution,
also return weight for I_0.
"""
function linear_weights(Δτ::AbstractFloat)
    local α, β, expΔτ
    if Δτ < 5e-4
        expΔτ = 1 - Δτ + 0.5*Δτ^2
        α = Δτ*(1/2 - Δτ/3)
        β = Δτ*(1/2 - Δτ/6)
    elseif Δτ > 50
        expΔτ = 0.0
        α = 1/Δτ
        β = 1.0 - α
    else
        expΔτ = exp(-Δτ)
        α = (1 - expΔτ)/Δτ - expΔτ
        β = 1 - α - expΔτ
    end
    return α, β, expΔτ
end

"""
    sample_from_extinction(atmos::Atmosphere,
                                λ0::Unitful.Length,
                                n_sites::Int)

Sample Voronoi sites by using continuum extinction as probalility density.
"""
function sample_from_extinction(atmos::Atmosphere,
                                λ0::Unitful.Length,
                                n_sites::Int)

    # Find continuum extinction and absorption extinction (without Thomson and Rayleigh)
    α_cont = α_continuum.(λ0,
                          atmos.temperature*1.0,
                          atmos.electron_density*1.0,
                          atmos.hydrogen_populations*1.0,
                          atmos.hydrogen_populations*1.0)

    positions = rejection_sampling(n_sites, atmos, ustrip.(α_cont))
    return positions
end
