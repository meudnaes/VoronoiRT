include("functions.jl")

#=
    function z_ray()

Ray moving upwards, upwind point intersecting with lower z plane. Assumes angle
in radians, θ ∈ [0, 2π) and ϕ ∈ [0, π].
=#
function z_ray(θ, ϕ, k, i, j, sign_z, sign_x, sign_y, I_0, S_0, α, atmos)

    # Centre point (c)
    x_c = atmos.x[idx]
    y_c = atmos.y[idy]
    z_c = atmos.z[idz]

    # Lower corner to interpolate from
    idz_i = idz + sign_z
    idx_i = idx + sign_x
    idy_i = idy + sign_y

    # Cut point
    idz_u = idz + 2*sign_z + 1

    # Length to plane
    Δz = sign_z*(atmos.z[idz_u] - z_c)

    # calculate the length to intersections with z, x, and y plane
    r = Δz/cos(ϕ)
    x_u = x_c + r*cos(θ)*sin(ϕ)
    y_u = y_c + r*sin(θ)*sin(ϕ)

    # Do linear interpolation to find τ and S

    # exinction at centre and upwind point
    α_c = α[idz, idx, idy]
    α_u = bilinear(z_u, x_u, y_u, idz_i, idx_i, idy_i, atmos, α)

    # Find the Δτ optical path from upwind to grid point
    Δτ_u = trapezoidal(r, α_x, α_u)

    # Find source function at grid point and upwind point (intepolate)
    S_c = S_0[idz, idx, idy]
    S_u = bilinear(z_u, x_u, y_u, idz_i, idx_i, idy_i, atmos, S_0)

    # From Hennicker et. al. (2020) equation 13
    e_0 = 1 - exp(-Δτ_u)
    e_1 = Δτ_u - e_0

    # From Hennicker et. al. (2020) equation 12, with ω = 1
    a_ijk = e_0 - e_1/Δτ_u
    b_ijk = e_1/Δτ_u
    # c_ijk = 0
    d_ijk = exp(-Δτ_u)

    # Interpolate to find intensity at upwind point
    I_u = bilinear(z_u, x_u, y_u, idz_i, idx_i, idy_i, atmos, I_0)

    # Integrate intensity from two-point quadrature
    I_ijk = a_ijk*S_u + b_ijk*S_c + d_ijk*I_u

    return I_ijk
end

#=
    function z_ray()

Ray moving upwards, upwind point intersecting with lower z plane. Assumes angle
in radians, θ ∈ [0, 2π) and ϕ ∈ [0, π].
=#
function x_ray(θ, ϕ, k, i, j, sign_z, sign_x, sign_y, I_0, S_0, α, atmos)

    # Centre point (c)
    x_c = atmos.x[idx]
    y_c = atmos.y[idy]
    z_c = atmos.z[idz]

    # Lower corner to interpolate from
    idz_i = idz + sign_z
    idx_i = idx + sign_x
    idy_i = idy + sign_y

    # Cut point
    idz_u = idz + 2*sign_z + 1

    # Length to plane
    Δz = sign_z*(atmos.z[idz_u] - z_c)

    # calculate the length to intersections with z, x, and y plane
    r = Δz/cos(ϕ)
    x_u = x_c + r*cos(θ)*sin(ϕ)
    y_u = y_c + r*sin(θ)*sin(ϕ)

    # Do linear interpolation to find τ and S

    # exinction at centre and upwind point
    α_c = α[idz, idx, idy]
    α_u = bilinear(z_u, x_u, y_u, idz_i, idx_i, idy_i, atmos, α)

    # Find the Δτ optical path from upwind to grid point
    Δτ_u = trapezoidal(r, α_x, α_u)

    # Find source function at grid point and upwind point (intepolate)
    S_c = S_0[idz, idx, idy]
    S_u = bilinear(z_u, x_u, y_u, idz_i, idx_i, idy_i, atmos, S_0)

    # From Hennicker et. al. (2020) equation 13
    e_0 = 1 - exp(-Δτ_u)
    e_1 = Δτ_u - e_0

    # From Hennicker et. al. (2020) equation 12, with ω = 1
    a_ijk = e_0 - e_1/Δτ_u
    b_ijk = e_1/Δτ_u
    # c_ijk = 0
    d_ijk = exp(-Δτ_u)

    # Interpolate to find intensity at upwind point
    I_u = bilinear(z_u, x_u, y_u, idz_i, idx_i, idy_i, atmos, I_0)

    # Integrate intensity from two-point quadrature
    I_ijk = a_ijk*S_u + b_ijk*S_c + d_ijk*I_u

    return I_ijk
end

#=
    function z_ray()

Ray moving upwards, upwind point intersecting with lower z plane. Assumes angle
in radians, θ ∈ [0, 2π) and ϕ ∈ [0, π].
=#
function y_ray(θ, ϕ, k, i, j, sign_z, sign_x, sign_y, I_0, S_0, α, atmos)

    # Centre point (c)
    x_c = atmos.x[idx]
    y_c = atmos.y[idy]
    z_c = atmos.z[idz]

    # Lower corner to interpolate from
    idz_i = idz + sign_z
    idx_i = idx + sign_x
    idy_i = idy + sign_y

    # Cut point
    idz_u = idz + 2*sign_z + 1

    # Length to plane
    Δz = sign_z*(atmos.z[idz_u] - z_c)

    # calculate the length to intersections with z, x, and y plane
    r = Δz/cos(ϕ)
    x_u = x_c + r*cos(θ)*sin(ϕ)
    y_u = y_c + r*sin(θ)*sin(ϕ)

    # Do linear interpolation to find τ and S

    # exinction at centre and upwind point
    α_c = α[idz, idx, idy]
    α_u = bilinear(z_u, x_u, y_u, idz_i, idx_i, idy_i, atmos, α)

    # Find the Δτ optical path from upwind to grid point
    Δτ_u = trapezoidal(r, α_x, α_u)

    # Find source function at grid point and upwind point (intepolate)
    S_c = S_0[idz, idx, idy]
    S_u = bilinear(z_u, x_u, y_u, idz_i, idx_i, idy_i, atmos, S_0)

    # From Hennicker et. al. (2020) equation 13
    e_0 = 1 - exp(-Δτ_u)
    e_1 = Δτ_u - e_0

    # From Hennicker et. al. (2020) equation 12, with ω = 1
    a_ijk = e_0 - e_1/Δτ_u
    b_ijk = e_1/Δτ_u
    # c_ijk = 0
    d_ijk = exp(-Δτ_u)

    # Interpolate to find intensity at upwind point
    I_u = bilinear(z_u, x_u, y_u, idz_i, idx_i, idy_i, atmos, I_0)

    # Integrate intensity from two-point quadrature
    I_ijk = a_ijk*S_u + b_ijk*S_c + d_ijk*I_u

    return I_ijk
end
