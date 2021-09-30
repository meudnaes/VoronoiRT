include("functions.jl")

function short_characteristic_ray(θ, φ, idz, idx, idy, atmos, α, I_0, S_0)
    # This function is not general for all angels. Only for angles hitting
    # the box in top or bottom wall.
    if θ < 45
        idz_u = idz + 1
    elseif 135 < θ < 180
        idz_u = idz - 1
    else
        println("θ angle is not valid: $θ")
    end

    # centre point (c)
    x_c = atmos.x[idx]
    y_c = atmos.y[idy]
    z_c = atmos.z[idz]

    # fix angle stuff. And keep track of corners. I always use the "lower left"
    # corner for interpolation. Therefore I keep track of the indices for the
    # lower left corner: iiz, idx_u, and idy_u...
    if φ < 90
        sign_x = 1
        sign_y = 1
        idx_u = idx
        idy_u = idy
    elseif 90 < φ > 180
        sign_x = -1
        sign_y = 1
        idx_u = idx-1
        idy_u = idy
    elseif 180 < φ > 270
        sign_x = -1
        sign_y = -1
        idx_u = idx-1
        idy_u = idy-1
    elseif 270 < φ > 360
        sign_x = 1
        sign_y = -1
        idx_u = idx
        idy_u = idy-1
    else
        println("φ angle is not valid: $φ")
    end

    # convert to radians
    θ = θ*pi/180
    φ = φ*pi/180

    # find upwind (u) points
    z_u = atmos.z[idz_u]
    Δz = z_c - z_u
    # in the Cartesian grid
    r = Δz/cos(θ)
    Δx = r*cos(φ)*sin(θ)
    Δy = r*sin(φ)*sin(θ)
    x_u = atmos.x[idx] + sign_x*Δx
    y_u = atmos.y[idy] + sign_y*Δy

    # How to find S ? We know the initial value on the whole grid.
    # At first iteration, do linear interpolation on initial values.
    # Same for τ

    # exinction at centre
    α_p = α[idz, idx, idy]

    # exinction at upwind
    α_u = bilinear(z_u, x_u, y_u, idz_u, idx_u, idy_u, atmos, α)

    # Find the Δτ optical path from upwind to grid point
    Δτ_u = trapezoidal(sqrt(Δz^2 + Δx^2 + Δy^2), α_p, α_u)

    # Find source function at grid point and upwind point (intepolate)
    S_p = S_0[idz, idx, idy]
    S_u = bilinear(z_u, x_u, y_u, idz_u, idx_u, idy_u, atmos, S_0)

    # From Hennicker et. al. (2020) equation 13
    e_0 = 1 - exp(-Δτ_u)
    e_1 = Δτ_u - e_0

    # From Hennicker et. al. (2020) equation 12, with ω = 1
    a_ijk = e_0 - e_1/Δτ_u
    b_ijk = e_1/Δτ_u
    # c_ijk = 0
    d_ijk = exp(-Δτ_u)

    # Interpolate to find intensity at upwind point
    I_u = bilinear(z_u, x_u, y_u, idz_u, idx_u, idy_u, atmos, I_0)

    # Integrate intensity from two-point quadrature
    I_ijk = a_ijk*S_u + b_ijk*S_p + d_ijk*I_u
    return I_ijk
end

function J_ν(weights, intensities)
    # Gaussian quadrature to calculate mean intensity
    J = 0u"kW*m^-2*nm^-1"
    for (weight, intensity) in zip(weights, intensitites)
        J = J + weight*intensity
    end
    return J::Float64
end
