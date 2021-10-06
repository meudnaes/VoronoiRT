include("functions.jl")

#=
    function z_ray()

Ray moving upwards, upwind point intersecting with lower z plane. Assumes angle
in radians, θ ∈ [0, 2π) and ϕ ∈ [0, π].
=#
function z_up_ray(θ, ϕ, idz, sign_x, sign_y, I_0, S_0, α, atmos)
    # Allocate space for intensity
    I = zero(I_0)

    # Z center and upwind position
    z_c = atmos.z[idz]
    idz_u = idz - 1
    z_u = atmos.z[idz_u]

    for idx in 1:length(atmos.x)
        for idy in 1:length(atmos.y)
            # Centre point (c)
            x_c = atmos.x[idx]
            y_c = atmos.y[idy]

            # Lower corner to interpolate from
            idx_i = idx + sign_x
            idy_i = idy + sign_y

            # Length to plane
            Δz = z_c - z_u

            # calculate the length to intersections with z, x, and y plane
            r = Δz/cos(ϕ)
            x_u = x_c + r*cos(θ)*sin(ϕ)
            y_u = y_c + r*sin(θ)*sin(ϕ)

            # Do linear interpolation to find τ and S

            # exinction at centre and upwind point
            α_c = α[idz, idx, idy]
            α_u = bilinear_xy(x_u, y_u, idx_i, idy_i, atmos, α[idz_u])

            # Find the Δτ optical path from upwind to grid point
            Δτ_u = trapezoidal(r, α_c, α_u)

            # Find source function at grid point and upwind point (intepolate)
            S_c = S_0[idz, idx, idy]
            S_u = bilinear_xy(x_u, y_u, idx_i, idy_i, atmos, S_0[idz_u])

            # From Hennicker et. al. (2020) equation 13
            e_0 = 1 - exp(-Δτ_u)
            e_1 = Δτ_u - e_0

            # From Hennicker et. al. (2020) equation 12, with ω = 1
            a_ijk = e_0 - e_1/Δτ_u
            b_ijk = e_1/Δτ_u
            # c_ijk = 0
            d_ijk = exp(-Δτ_u)

            # Interpolate to find intensity at upwind point
            I_u = bilinear_xy(x_u, y_u, idx_i, idy_i, atmos, I_0)

            # Integrate intensity from two-point quadrature
            I[idx, idy] = a_ijk*S_u + b_ijk*S_c + d_ijk*I_u
        end
    end
    return I
end

function z_down_ray(θ, ϕ, idz, sign_x, sign_y, I_0, S_0, α, atmos)
    # Allocate space for intensity
    I = zero(I_0)

    # Z center and upwind position
    z_c = atmos.z[idz]
    idz_u = idz + 1
    z_u = atmos.z[idz_u]

    for idx in 1:length(atmos.x)
        for idy in 1:length(atmos.y)
            # Centre point (c)
            x_c = atmos.x[idx]
            y_c = atmos.y[idy]

            # Lower corner to interpolate from
            idx_i = idx + sign_x
            idy_i = idy + sign_y

            # Length to plane
            Δz = z_u - z_c

            # calculate the length to intersections with z, x, and y plane
            r = Δz/cos(ϕ)
            x_u = x_c + r*cos(θ)*sin(ϕ)
            y_u = y_c + r*sin(θ)*sin(ϕ)

            # Do linear interpolation to find τ and S
            # exinction at centre and upwind point
            α_c = α[idz, idx, idy]
            α_u = bilinear_xy(x_u, y_u, idx_i, idy_i, atmos, α[idz_u, :, :])

            # Find the Δτ optical path from upwind to grid point
            Δτ_u = trapezoidal(r, α_c, α_u)

            # Find source function at grid point and upwind point (intepolate)
            S_c = S_0[idz, idx, idy]
            S_u = bilinear_xy(x_u, y_u, idx_i, idy_i, atmos, S_0[idz_u, :, :])

            # From Hennicker et. al. (2020) equation 13
            e_0 = 1 - exp(-Δτ_u)
            e_1 = Δτ_u - e_0

            # From Hennicker et. al. (2020) equation 12, with ω = 1
            a_ijk = e_0 - e_1/Δτ_u
            b_ijk = e_1/Δτ_u
            # c_ijk = 0
            d_ijk = exp(-Δτ_u)

            # Interpolate to find intensity at upwind point
            I_u = bilinear_xy(x_u, y_u, idx_i, idy_i, atmos, I_0)

            # Integrate intensity from two-point quadrature
            I[idx, idy] = a_ijk*S_u + b_ijk*S_c + d_ijk*I_u
        end
    end

    return I
end

#=
    function x_up_ray()

Ray moving upwards, upwind point intersecting with lower z plane. Assumes angle
in radians, θ ∈ [0, 2π) and ϕ ∈ [0, π].
=#
function x_up_ray(θ, ϕ, idz, sign_x, sign_y, I_0, S_0, α, atmos)

    nx = length(atmos.x)

    # Allocate space for intensity
    I = zero(I_0)
    I_upper = zero(I_0[end])

    # Z center and interpolation position
    z_c = atmos.z[idz]
    idz_i = idz - 1

    # Do the ghost cell thing...
    idx = nx
    for idy in 1:length(atmos.y)
        # Centre point (c)
        x_c = atmos.x[idx]
        y_c = atmos.y[idy]

        # Lower corner to interpolate from
        idy_i = idy + sign_y

        # X upwind position
        idx_u = idx - 2*sign_x + 1
        if idx_u > nx
            idx_u = 1
        end
        x_u = atmos.x[idx_u]

        # Length to plane
        Δx = abs(x_c - x_u)

        # calculate the length to upwind position
        r = Δx/(cos(θ)*sin(ϕ))
        z_u = z_c + r*cos(ϕ)
        y_u = y_c + r*sin(θ)*sin(ϕ)

        # Do linear interpolation to find τ and S

        # exinction at centre and upwind point
        α_c = α[idz, idx, idy]
        α_vals = [α[idz_i, idx_i, idy_i], α[idz_i, idx_i, idy],
                  α[idz, idx_i, idy_i], α[idz, idx_i, idy]]
        α_u = bilinear_yz(z_u, x_u, y_u, idz_u, idx_i, idy_i, atmos, α_vals)

        # Find the Δτ optical path from upwind to grid point
        Δτ_u = trapezoidal(r, α_c, α_u)

        # Find source function at grid point and upwind point (intepolate)
        S_c = S_0[idz, idx, idy]
        S_vals = [S_0[idz_i, idx_i, idy_i], S_0[idz_i, idx_i, idy],
                  S_0[idz, idx_i, idy_i], S_0[idz, idx_i, idy]]
        S_u = bilinear_yz(z_u, x_u, y_u, idz_u, idx_i, idy_i, atmos, S_0)

        # From Hennicker et. al. (2020) equation 13
        e_0 = 1 - exp(-Δτ_u)
        e_1 = Δτ_u - e_0

        # From Hennicker et. al. (2020) equation 12, with ω = 1
        a_ijk = e_0 - e_1/Δτ_u
        b_ijk = e_1/Δτ_u
        # c_ijk = 0
        d_ijk = exp(-Δτ_u)

        # Interpolate to find intensity at upwind point
        I_vals = [I_0[idx_i, idy_i], I_0[idx_i, idy],
                   S_0[idx_i, idy_i], S_0[idx_i, idy]]

        I_ghost = bilinear_yz(z_u, y_u, idz_u, idy_i, atmos, I_ghost)

        # Integrate intensity from two-point quadrature
        I_upper[idy_i] = a_ijk*S_u + b_ijk*S_c + d_ijk*I_ghost
    end

    for idx in 1:length(atmos.x)
        for idy in 1:length(atmos.y)
            # Centre point (c)
            x_c = atmos.x[idx]
            y_c = atmos.y[idy]

            # Lower corner to interpolate from
            idy_i = idy + sign_y

            # X upwind position
            idx_u = idx - 2*sign_x + 1
            if idx_u > nx
                idx_u = 1
            end
            x_u = atmos.x[idx_u]

            # Length to plane
            Δx = abs(x_c - x_u)

            # calculate the length to upwind position
            r = Δx/(cos(θ)*sin(ϕ))
            z_u = z_c + r*cos(ϕ)
            y_u = y_c + r*sin(θ)*sin(ϕ)

            # Do linear interpolation to find τ and S

            # exinction at centre and upwind point
            α_c = α[idz, idx, idy]
            α_vals = [α[idz_i, idx_i, idy_i], α[idz_i, idx_i, idy],
                      α[idz, idx_i, idy_i], α[idz, idx_i, idy]]
            α_u = bilinear_yz(z_u, x_u, y_u, idz_u, idx_i, idy_i, atmos, α_vals)

            # Find the Δτ optical path from upwind to grid point
            Δτ_u = trapezoidal(r, α_c, α_u)

            # Find source function at grid point and upwind point (intepolate)
            S_c = S_0[idz, idx, idy]
            S_vals = [S_0[idz_i, idx_i, idy_i], S_0[idz_i, idx_i, idy],
                      S_0[idz, idx_i, idy_i], S_0[idz, idx_i, idy]]
            S_u = bilinear_yz(z_u, x_u, y_u, idz_u, idx_i, idy_i, atmos, S_0)

            # From Hennicker et. al. (2020) equation 13
            e_0 = 1 - exp(-Δτ_u)
            e_1 = Δτ_u - e_0

            # From Hennicker et. al. (2020) equation 12, with ω = 1
            a_ijk = e_0 - e_1/Δτ_u
            b_ijk = e_1/Δτ_u
            # c_ijk = 0
            d_ijk = exp(-Δτ_u)

            # Interpolate to find intensity at upwind point
            I_vals = [I_0[idx_i, idy_i], I_0[idx_i, idy],
                      I_upper[idy_i], I_upper[idy]]

            I_u = bilinear_yz(z_u, y_u, idz_u, idy_i, atmos, I_upwind)

            # Integrate intensity from two-point quadrature
            I[idx, idy] = a_ijk*S_u + b_ijk*S_c + d_ijk*I_u
        end
        # Update upper interpolate value
        I_upper = I[idx]
    end
    return I
end

#=
    function z_ray()

Ray moving upwards, upwind point intersecting with lower z plane. Assumes angle
in radians, θ ∈ [0, 2π) and ϕ ∈ [0, π].
=#
function y_ray(θ, ϕ, idz, idx, idy, sign_z, sign_x, sign_y, I_0, S_0, α, atmos)

    ny = length(atmos.y)

    # Centre point (c)
    x_c = atmos.x[idx]
    y_c = atmos.y[idy]
    z_c = atmos.z[idz]

    # Lower corner to interpolate from
    idz_i = idz + sign_z
    idx_i = idx + sign_x
    idy_i = idy + sign_y

    # Cut point
    idy_u = idy + 2*sign_y + 1

    if idy_u > ny
        idy_u = 1
    end

    y_u = atmos.y[idy_u]

    # Length to plane
    Δy = sign_y*(atmos.y[idy_u] - y_c)

    # calculate the length to intersections with z and x plane
    r = Δy/(sin(θ)*sin(ϕ))
    z_u = z_c + r*cos(ϕ)
    x_u = x_c + r*cos(θ)*sin(ϕ)

    # Do linear interpolation to find τ and S

    # exinction at centre and upwind point
    α_c = α[idz, idx, idy]
    α_u = bilinear_xy(z_u, x_u, y_u, idz_i, idx_i, idy_i, atmos, α)

    # Find the Δτ optical path from upwind to grid point
    Δτ_u = trapezoidal(r, α_c, α_u)

    # Find source function at grid point and upwind point (intepolate)
    S_c = S_0[idz, idx, idy]
    S_u = bilinear_xy(z_u, x_u, y_u, idz_i, idx_i, idy_i, atmos, S_0)

    # From Hennicker et. al. (2020) equation 13
    e_0 = 1 - exp(-Δτ_u)
    e_1 = Δτ_u - e_0

    # From Hennicker et. al. (2020) equation 12, with ω = 1
    a_ijk = e_0 - e_1/Δτ_u
    b_ijk = e_1/Δτ_u
    # c_ijk = 0
    d_ijk = exp(-Δτ_u)

    # Interpolate to find intensity at upwind point
    I_u = bilinear_xz(z_u, x_u, y_u, idz_i, idx_i, idy_i, atmos, I_0)

    # Integrate intensity from two-point quadrature
    I_ijk = a_ijk*S_u + b_ijk*S_c + d_ijk*I_u

    return I_ijk
end
