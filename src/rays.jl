include("functions.jl")

#=
    function z_up_ray(θ::AbstractFloat, ϕ::AbstractFloat, idz::Int, sign_x::Int,
                       sign_y::Int, I_0::AbstractArray, S_0::AbstractArray,
                       α::AbstractArray, atmos::Atmosphere)

Ray moving upwards, upwind point intersecting with lower xy plane. Assumes angle
in radians.
=#
function xy_up_ray(θ::AbstractFloat, ϕ::AbstractFloat, idz::Int, sign_x::Int,
                   sign_y::Int, I_0::AbstractArray, S_0::AbstractArray,
                   α::AbstractArray, atmos::Atmosphere)

    nx = length(atmos.x)
    ny = length(atmos.y)

    # Allocate space for intensity
    I = zero(I_0)
    # Z center and upwind position
    z_c = atmos.z[idz]
    idz_u = idz - 1
    z_u = atmos.z[idz_u]

    # Length to plane
    Δz = z_u - z_c
    r = Δz/cos(θ)

    x_increment = r*cos(ϕ)*sin(θ)
    y_increment = r*sin(ϕ)*sin(θ)

    α_centre = α[idz,:,:]
    α_upwind = α[idz_u,:,:]

    S_centre = S_0[idz,:,:]
    S_upwind = S_0[idz_u,:,:]

    # Loop direction
    start_x, stop_x = range_bounds(-sign_x, nx)
    start_y, stop_y = range_bounds(-sign_y, ny)

    for idx in start_x:-sign_x:stop_x
        for idy in start_y:-sign_y:stop_y
            # Centre point (c)
            x_c = atmos.x[idx]
            y_c = atmos.y[idy]

            # Lower corner to interpolate from
            idx_i = Int(idx + (sign_x-1)/2)
            idy_i = Int(idy + (sign_y-1)/2)

            # Length to plane
            Δz = z_c - z_u

            # calculate the length to intersections with z, x, and y plane
            r = Δz/cos(θ)
            x_u = x_c + x_increment
            y_u = y_c + y_increment

            # Do linear interpolation to find τ and S

            # exinction at centre and upwind point
            α_c = α_centre[idx, idy]
            α_u = bilinear_xy(x_u, y_u, idx_i, idy_i, atmos, α_upwind)

            # Find the Δτ optical path from upwind to grid point
            Δτ_u = trapezoidal(r, α_c, α_u)

            # Find source function at grid point and upwind point (intepolate)
            S_c = S_centre[idx, idy]
            S_u = bilinear_xy(x_u, y_u, idx_i, idy_i, atmos, S_upwind)

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
    function xy_down_ray(θ::AbstractFloat, ϕ::AbstractFloat, idz::Int, sign_x::Int,
                         sign_y::Int, I_0::AbstractArray, S_0::AbstractArray,
                         α::AbstractArray, atmos::Atmosphere)

Ray moving downwards, upwind point intersecting with upper xy plane. Assumes
angle in radians.
=#
function xy_down_ray(θ::AbstractFloat, ϕ::AbstractFloat, idz::Int, sign_x::Int,
                     sign_y::Int, I_0::AbstractArray, S_0::AbstractArray,
                     α::AbstractArray, atmos::Atmosphere)

    nx = length(atmos.x)
    ny = length(atmos.y)

    # Allocate space for intensity
    I = zero(I_0)

    # Z center and upwind position
    z_c = atmos.z[idz]
    idz_u = idz + 1
    z_u = atmos.z[idz_u]

    # Length to plane
    Δz = z_u - z_c
    r = Δz/cos(θ)

    x_increment = r*cos(ϕ)*sin(θ)
    y_increment = r*sin(ϕ)*sin(θ)

    α_centre = α[idz,:,:]
    α_upwind = α[idz_u,:,:]

    S_centre = S_0[idz,:,:]
    S_upwind = S_0[idz_u,:,:]

    # Loop direction
    start_x, stop_x = range_bounds(-sign_x, nx)
    start_y, stop_y = range_bounds(-sign_y, ny)

    for idx in start_x:-sign_x:stop_x
        for idy in start_y:-sign_y:stop_y
            # Centre point (c)
            x_c = atmos.x[idx]
            y_c = atmos.y[idy]

            # Lower corner to interpolate from
            idx_i = idx + Int((sign_x-1)/2)
            idy_i = idy + Int((sign_y-1)/2)

            # calculate the length to intersections with z, x, and y plane
            x_u = x_c + x_increment
            y_u = y_c + y_increment

            # Do linear interpolation to find τ and S
            # exinction at centre and upwind point
            α_c = α_centre[idx, idy]
            α_u = bilinear_xy(x_u, y_u, idx_i, idy_i, atmos, α_upwind)

            # Find the Δτ optical path from upwind to grid point
            Δτ_u = trapezoidal(r, α_c, α_u)

            # Find source function at grid point and upwind point (intepolate)
            S_c = S_centre[idx, idy]
            S_u = bilinear_xy(x_u, y_u, idx_i, idy_i, atmos, S_upwind)

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
    function yz_up_ray(θ::AbstractFloat, ϕ::AbstractFloat, idz::Int, sign_x::Int,
                       sign_y::Int, I_0::AbstractArray, S_0::AbstractArray,
                       α::AbstractArray, atmos::Atmosphere)

Ray moving upwards, upwind point intersecting with yz plane. Assumes angle
in radians.
=#
function yz_up_ray(θ::AbstractFloat, ϕ::AbstractFloat, idz::Int, sign_x::Int,
                   sign_y::Int, I_0::AbstractArray, S_0::AbstractArray,
                   α::AbstractArray, atmos::Atmosphere)

    nx = length(atmos.x)
    ny = length(atmos.y)

    # Loop direction
    start_x, stop_x = range_bounds(-sign_x, nx)
    start_y, stop_y = range_bounds(-sign_y, ny)

    # Allocate space for intensity
    I = zero(I_0)
    I_upper = zero(I_0[end,:])

    # Z center and interpolation position
    z_c = atmos.z[idz]
    idz_i = idz - 1

    # Length to plane
    Δx = atmos.x[2] - atmos.x[1]

    # calculate the length to upwind position
    r = Δx/(cos(ϕ)*sin(θ))
    z_increment = r*cos(θ)
    y_increment = r*sin(ϕ)*sin(θ)

    z_u = z_c + z_increment

    α_centre = α[idz,:,:]
    α_upwind = α[idz_i,:,:]

    S_centre = S_0[idz,:,:]
    S_upwind = S_0[idz_i,:,:]

    # Do the ghost cell thing...
    idx = stop_x
    x_c = atmos.x[nx]

    # X upwind position
    idx_u = roll(idx+sign_x, nx)
    x_u = atmos.x[idx_u]


    for idy in start_y:-sign_y:stop_y
        # Centre point (c)
        y_c = atmos.y[idy]

        # Lower corner to interpolate from
        idy_i = roll(idy+Int((sign_y-1)/2), ny)

        y_u = y_c + y_increment

        # Do linear interpolation to find τ and S

        # exinction at centre and upwind point
        α_c = α_centre[idx, idy]
        α_vals = [α_upwind[idx_u, idy_i]    α_upwind[idx_u, idy]
                  α_centre[idx_u, idy_i]    α_centre[idx_u, idy]]
        α_u = bilinear_yz(z_u, y_u, idz_i, idy_i, atmos, α_vals)

        # Find the Δτ optical path from upwind to grid point
        Δτ_u = trapezoidal(r, α_c, α_u)

        # Find source function at grid point and upwind point (intepolate)
        S_c = S_centre[idx, idy]
        S_vals = [S_upwind[idx_u, idy_i]    S_upwind[idx_u, idy]
                  S_centre[idx_u, idy_i]    S_centre[idx_u, idy]]
        S_u = bilinear_yz(z_u, y_u, idz_i, idy_i, atmos, S_vals)

        # From Hennicker et. al. (2020) equation 13
        e_0 = 1 - exp(-Δτ_u)
        e_1 = Δτ_u - e_0

        # From Hennicker et. al. (2020) equation 12, with ω = 1
        a_ijk = e_0 - e_1/Δτ_u
        b_ijk = e_1/Δτ_u
        # c_ijk = 0
        d_ijk = exp(-Δτ_u)

        # Interpolate to find intensity at upwind point
        I_vals = [I_0[idx_u, idy_i]         I_0[idx_u, idy]
                  S_centre[idx_u, idy_i]    S_centre[idx_u, idy]]

        I_ghost = bilinear_yz(z_u, y_u, idz_i, idy_i, atmos, I_vals)

        # Integrate intensity from two-point quadrature
        I_upper[idy_i] = a_ijk*S_u + b_ijk*S_c + d_ijk*I_ghost
    end

    for idx in start_x:-sign_x:stop_x
        x_c = atmos.x[idx]
        # x upwind position
        idx_u = roll(idx+sign_x, nx)
        x_u = atmos.x[idx_u]

        for idy in start_y:-sign_y:stop_y
            # Centre point (c)
            y_c = atmos.y[idy]

            # Lower corner to interpolate from
            idy_i = roll(idy+sign_y, ny)

            y_u = y_c + y_increment

            # Do linear interpolation to find τ and S

            # exinction at centre and upwind point
            α_c = α_centre[idx, idy]
            α_vals = [α_upwind[idx_u, idy_i]   α_upwind[idx_u, idy]
                      α_centre[idx_u, idy_i]   α_centre[idx_u, idy]]
            α_u = bilinear_yz(z_u, y_u, idz_i, idy_i, atmos, α_vals)

            # Find the Δτ optical path from upwind to grid point
            Δτ_u = trapezoidal(r, α_c, α_u)

            # Find source function at grid point and upwind point (intepolate)
            S_c = S_centre[idx, idy]
            S_vals = [S_upwind[idx_u, idy_i]    S_upwind[idx_u, idy]
                      S_centre[idx_u, idy_i]    S_centre[idx_u, idy]]
            S_u = bilinear_yz(z_u, y_u, idz_i, idy_i, atmos, S_vals)

            # From Hennicker et. al. (2020) equation 13
            e_0 = 1 - exp(-Δτ_u)
            e_1 = Δτ_u - e_0

            # From Hennicker et. al. (2020) equation 12, with ω = 1
            a_ijk = e_0 - e_1/Δτ_u
            b_ijk = e_1/Δτ_u
            # c_ijk = 0
            d_ijk = exp(-Δτ_u)

            # Interpolate to find intensity at upwind point
            I_vals = [I_0[idx_u, idy_i]     I_0[idx_u, idy]
                      I_upper[idy_i]        I_upper[idy]   ]
            I_u = bilinear_yz(z_u, y_u, idz_i, idy_i, atmos, I_vals)

            # Integrate intensity from two-point quadrature
            I[idx, idy] = a_ijk*S_u + b_ijk*S_c + d_ijk*I_u
        end
        # Update upper interpolate value
        I_upper = I[idx, :]
    end
    return I
end

#=
    function yz_down_ray(θ::AbstractFloat, ϕ::AbstractFloat, idz::Int, sign_x::Int,
                         sign_y::Int, I_0::AbstractArray, S_0::AbstractArray,
                         α::AbstractArray, atmos::Atmosphere)

Ray moving downwards, upwind point intersecting with yz plane. Assumes angle
in radians
=#
function yz_down_ray(θ::AbstractFloat, ϕ::AbstractFloat, idz::Int, sign_x::Int,
                       sign_y::Int, I_0::AbstractArray, S_0::AbstractArray,
                       α::AbstractArray, atmos::Atmosphere)
    nx = length(atmos.x)
    ny = length(atmos.y)

    # Loop direction
    start_x, stop_x = range_bounds(-sign_x, nx)
    start_y, stop_y = range_bounds(-sign_y, ny)

    # Allocate space for intensity
    I = zero(I_0)
    I_lower = zero(I_0[end,:])

    # Z center and interpolation position
    z_c = atmos.z[idz]
    idz_i = idz + 1

    # Length to plane
    Δx = atmos.x[2] - atmos.x[1]

    # calculate the length to upwind position
    r = Δx/(cos(ϕ)*sin(θ))
    z_increment = r*cos(θ)
    y_increment = r*sin(ϕ)*sin(θ)

    z_u = z_c + z_increment

    α_centre = α[idz,:,:]
    α_upwind = α[idz_i,:,:]

    S_centre = S_0[idz,:,:]
    S_upwind = S_0[idz_i,:,:]

    # Do the ghost cell thing...
    idx = stop_x
    x_c = atmos.x[nx]

    # X upwind position
    idx_u = roll(idx+sign_x, nx)
    x_u = atmos.x[idx_u]


    for idy in start_y:-sign_y:stop_y
        # Centre point (c)
        y_c = atmos.y[idy]

        # Lower corner to interpolate from
        idy_i = roll(idy+Int((sign_y-1)/2), ny)

        y_u = y_c + y_increment

        # Do linear interpolation to find τ and S

        # exinction at centre and upwind point
        α_c = α_centre[idx, idy]
        α_vals = [α_centre[idx_u, idy_i]    α_centre[idx_u, idy]
                  α_upwind[idx_u, idy_i]    α_upwind[idx_u, idy]]
        α_u = bilinear_yz(z_u, y_u, idz_i, idy_i, atmos, α_vals)

        # Find the Δτ optical path from upwind to grid point
        Δτ_u = trapezoidal(r, α_c, α_u)

        # Find source function at grid point and upwind point (intepolate)
        S_c = S_centre[idx, idy]
        S_vals = [S_centre[idx_u, idy_i]    S_centre[idx_u, idy]
                  S_upwind[idx_u, idy_i]    S_upwind[idx_u, idy]]
        S_u = bilinear_yz(z_u, y_u, idz_i, idy_i, atmos, S_vals)

        # From Hennicker et. al. (2020) equation 13
        e_0 = 1 - exp(-Δτ_u)
        e_1 = Δτ_u - e_0

        # From Hennicker et. al. (2020) equation 12, with ω = 1
        a_ijk = e_0 - e_1/Δτ_u
        b_ijk = e_1/Δτ_u
        # c_ijk = 0
        d_ijk = exp(-Δτ_u)

        # Interpolate to find intensity at upwind point
        I_vals = [S_centre[idx_u, idy_i]    S_centre[idx_u, idy]
                  I_0[idx_u, idy_i]         I_0[idx_u, idy]]

        I_ghost = bilinear_yz(z_u, y_u, idz_i, idy_i, atmos, I_vals)

        # Integrate intensity from two-point quadrature
        I_lower[idy_i] = a_ijk*S_u + b_ijk*S_c + d_ijk*I_ghost
    end

    for idx in start_x:-sign_x:stop_x
        x_c = atmos.x[idx]
        # X upwind position
        idx_u = roll(idx+sign_x, nx)
        x_u = atmos.x[idx_u]

        for idy in start_y:-sign_y:stop_y
            # Centre point (c)
            y_c = atmos.y[idy]

            # Lower corner to interpolate from
            idy_i = roll(idy+sign_y, ny)

            y_u = y_c + r*sin(ϕ)*sin(θ)

            # Do linear interpolation to find τ and S

            # exinction at centre and upwind point
            α_c = α_centre[idx, idy]
            α_vals = [α_centre[idx_u, idy_i]   α_centre[idx_u, idy]
                      α_upwind[idx_u, idy_i]   α_upwind[idx_u, idy]]
            α_u = bilinear_yz(z_u, y_u, idz_i, idy_i, atmos, α_vals)

            # Find the Δτ optical path from upwind to grid point
            Δτ_u = trapezoidal(r, α_c, α_u)

            # Find source function at grid point and upwind point (intepolate)
            S_c = S_centre[idx, idy]
            S_vals = [S_centre[idx_u, idy_i]    S_centre[idx_u, idy]
                      S_upwind[idx_u, idy_i]    S_upwind[idx_u, idy]]
            S_u = bilinear_yz(z_u, y_u, idz_i, idy_i, atmos, S_vals)

            # From Hennicker et. al. (2020) equation 13
            e_0 = 1 - exp(-Δτ_u)
            e_1 = Δτ_u - e_0

            # From Hennicker et. al. (2020) equation 12, with ω = 1
            a_ijk = e_0 - e_1/Δτ_u
            b_ijk = e_1/Δτ_u
            # c_ijk = 0
            d_ijk = exp(-Δτ_u)

            # Interpolate to find intensity at upwind point
            I_vals = [I_0[idx_u, idy_i]     I_0[idx_u, idy]
                      I_lower[idy_i]        I_lower[idy]   ]
            I_u = bilinear_yz(z_u, y_u, idz_i, idy_i, atmos, I_vals)

            # Integrate intensity from two-point quadrature
            I[idx, idy] = a_ijk*S_u + b_ijk*S_c + d_ijk*I_u
        end
        # Update upper interpolate value
        I_lower = I[idx, :]
    end
    return I
end

#=
    function xz_up_ray(θ::AbstractFloat, ϕ::AbstractFloat, idz::Int, sign_x::Int,
                       sign_y::Int, I_0::AbstractArray, S_0::AbstractArray,
                       α::AbstractArray, atmos::Atmosphere)

Ray moving upwards, upwind point intersecting with xz plane. Assumes angle
in radians.
=#
function xz_up_ray(θ::AbstractFloat, ϕ::AbstractFloat, idz::Int, sign_x::Int,
                   sign_y::Int, I_0::AbstractArray, S_0::AbstractArray,
                   α::AbstractArray, atmos::Atmosphere)

    nx = length(atmos.x)
    ny = length(atmos.y)

    # Loop direction
    start_x, stop_x = range_bounds(-sign_x, nx)
    start_y, stop_y = range_bounds(-sign_y, ny)

    # Allocate space for intensity
    I = zero(I_0)
    I_upper = zero(I_0[:,end])

    # Z center and interpolation position
    z_c = atmos.z[idz]
    idz_i = idz - 1

    # Length to plane
    Δy = atmos.y[2] - atmos.y[1]

    # calculate the length to upwind position
    r = Δy/(sin(ϕ)*sin(θ))
    z_increment = r*cos(θ)
    x_increment = r*cos(ϕ)*sin(θ)

    z_u = z_c + z_increment

    α_centre = α[idz,:,:]
    α_upwind = α[idz_i,:,:]

    S_centre = S_0[idz,:,:]
    S_upwind = S_0[idz_i,:,:]

    # Do the ghost cell thing...
    idy = stop_y
    y_c = atmos.y[ny]

    # X upwind position
    idy_u = roll(idy+sign_y, ny)
    y_u = atmos.y[idy_u]


    for idx in start_x:-sign_x:stop_x
        # Centre point (c)
        x_c = atmos.x[idx]

        # Lower corner to interpolate from
        idx_i = roll(idx+Int((sign_x-1)/2), nx)

        x_u = x_c + x_increment

        # Do linear interpolation to find τ and S

        # exinction at centre and upwind point
        α_c = α_centre[idx, idy]
        α_vals = [α_upwind[idx_i, idy_u]    α_upwind[idx, idy_u]
                  α_centre[idx_i, idy_u]    α_centre[idx, idy_u]]
        α_u = bilinear_xz(z_u, x_u, idz_i, idx_i, atmos, α_vals)

        # Find the Δτ optical path from upwind to grid point
        Δτ_u = trapezoidal(r, α_c, α_u)

        # Find source function at grid point and upwind point (intepolate)
        S_c = S_centre[idx, idy]
        S_vals = [S_upwind[idx_i, idy_u]    S_upwind[idx, idy_u]
                  S_centre[idx_i, idy_u]    S_centre[idx, idy_u]]
        S_u = bilinear_xz(z_u, x_u, idz_i, idx_i, atmos, S_vals)

        # From Hennicker et. al. (2020) equation 13
        e_0 = 1 - exp(-Δτ_u)
        e_1 = Δτ_u - e_0

        # From Hennicker et. al. (2020) equation 12, with ω = 1
        a_ijk = e_0 - e_1/Δτ_u
        b_ijk = e_1/Δτ_u
        # c_ijk = 0
        d_ijk = exp(-Δτ_u)

        # Interpolate to find intensity at upwind point
        I_vals = [I_0[idx_i, idy_u]         I_0[idx, idy_u]
                  S_centre[idx_i, idy_u]    S_centre[idx, idy_u]]

        I_ghost = bilinear_xz(z_u, x_u, idz_i, idx_i, atmos, I_vals)

        # Integrate intensity from two-point quadrature
        I_upper[idx_i] = a_ijk*S_u + b_ijk*S_c + d_ijk*I_ghost
    end

    for idy in start_y:-sign_y:stop_y
        y_c = atmos.y[idy]
        # x upwind position
        idy_u = roll(idy+sign_y, ny)
        y_u = atmos.y[idy_u]

        for idx in start_x:-sign_x:stop_x
            # Centre point (c)
            x_c = atmos.x[idx]

            # Lower corner to interpolate from
            idx_i = roll(idx+sign_x, nx)

            x_u = x_c + x_increment

            # Do linear interpolation to find τ and S

            # exinction at centre and upwind point
            α_c = α_centre[idx, idy]
            α_vals = [α_upwind[idx_i, idy_u]   α_upwind[idx, idy_u]
                      α_centre[idx_i, idy_u]   α_centre[idx, idy_u]]
            α_u = bilinear_xz(z_u, x_u, idz_i, idx_i, atmos, α_vals)

            # Find the Δτ optical path from upwind to grid point
            Δτ_u = trapezoidal(r, α_c, α_u)

            # Find source function at grid point and upwind point (intepolate)
            S_c = S_centre[idx, idy]
            S_vals = [S_upwind[idx_i, idy_u]    S_upwind[idx, idy_u]
                      S_centre[idx_i, idy_u]    S_centre[idx, idy_u]]
            S_u = bilinear_xz(z_u, x_u, idz_i, idx_i, atmos, S_vals)

            # From Hennicker et. al. (2020) equation 13
            e_0 = 1 - exp(-Δτ_u)
            e_1 = Δτ_u - e_0

            # From Hennicker et. al. (2020) equation 12, with ω = 1
            a_ijk = e_0 - e_1/Δτ_u
            b_ijk = e_1/Δτ_u
            # c_ijk = 0
            d_ijk = exp(-Δτ_u)

            # Interpolate to find intensity at upwind point
            I_vals = [I_0[idx_i, idy_u]     I_0[idx, idy_u]
                      I_upper[idx_i]        I_upper[idx]   ]
            I_u = bilinear_xz(z_u, x_u, idz_i, idx_i, atmos, I_vals)

            # Integrate intensity from two-point quadrature
            I[idx, idy] = a_ijk*S_u + b_ijk*S_c + d_ijk*I_u
        end
        # Update upper interpolate value
        I_upper = I[:, idy]
    end
    return I
end

#=
    function xz_down_ray(θ::AbstractFloat, ϕ::AbstractFloat, idz::Int, sign_x::Int,
                           sign_y::Int, I_0::AbstractArray, S_0::AbstractArray,
                           α::AbstractArray, atmos::Atmosphere)

Ray moving downwards, upwind point intersecting with xz plane. Assumes angle
in radians.
=#
function xz_down_ray(θ::AbstractFloat, ϕ::AbstractFloat, idz::Int, sign_x::Int,
                       sign_y::Int, I_0::AbstractArray, S_0::AbstractArray,
                       α::AbstractArray, atmos::Atmosphere)
    nx = length(atmos.x)
    ny = length(atmos.y)

    # Loop direction
    start_x, stop_x = range_bounds(-sign_x, nx)
    start_y, stop_y = range_bounds(-sign_y, ny)

    # Allocate space for intensity
    I = zero(I_0)
    I_lower = zero(I_0[end,:])

    # Z center and interpolation position
    z_c = atmos.z[idz]
    idz_i = idz + 1

    # Length to plane
    Δy = atmos.y[2] - atmos.y[1]

    # calculate the length to upwind position
    r = Δy/(sin(ϕ)*sin(θ))
    z_increment = r*cos(θ)
    x_increment = r*cos(ϕ)*sin(θ)

    z_u = z_c + z_increment

    α_centre = α[idz,:,:]
    α_upwind = α[idz_i,:,:]

    S_centre = S_0[idz,:,:]
    S_upwind = S_0[idz_i,:,:]

    # Do the ghost cell thing...
    idy = stop_y
    y_c = atmos.y[ny]

    # X upwind position
    idy_u = roll(idy+sign_y, ny)
    y_u = atmos.y[idy_u]


    for idx in start_x:-sign_x:stop_x
        # Centre point (c)
        x_c = atmos.x[idx]

        # Lower corner to interpolate from
        idx_i = roll(idx+Int((sign_x-1)/2), nx)

        x_u = x_c + x_increment

        # Do linear interpolation to find τ and S

        # exinction at centre and upwind point
        α_c = α_centre[idx, idy]
        α_vals = [α_centre[idx_i, idy_u]    α_centre[idx, idy_u]
                  α_upwind[idx_i, idy_u]    α_upwind[idx, idy_u]]
        α_u = bilinear_xz(z_u, x_u, idz_i, idx_i, atmos, α_vals)

        # Find the Δτ optical path from upwind to grid point
        Δτ_u = trapezoidal(r, α_c, α_u)

        # Find source function at grid point and upwind point (intepolate)
        S_c = S_centre[idx, idy]
        S_vals = [S_centre[idx_i, idy_u]    S_centre[idx, idy_u]
                  S_upwind[idx_i, idy_u]    S_upwind[idx, idy_u]]
        S_u = bilinear_xz(z_u, x_u, idz_i, idx_i, atmos, S_vals)

        # From Hennicker et. al. (2020) equation 13
        e_0 = 1 - exp(-Δτ_u)
        e_1 = Δτ_u - e_0

        # From Hennicker et. al. (2020) equation 12, with ω = 1
        a_ijk = e_0 - e_1/Δτ_u
        b_ijk = e_1/Δτ_u
        # c_ijk = 0
        d_ijk = exp(-Δτ_u)

        # Interpolate to find intensity at upwind point
        I_vals = [S_centre[idx_i, idy_u]    S_centre[idx, idy_u]
                  I_0[idx_i, idy_u]         I_0[idx, idy_u]]

        I_ghost = bilinear_xz(z_u, x_u, idz_i, idx_i, atmos, I_vals)

        # Integrate intensity from two-point quadrature
        I_lower[idx_i] = a_ijk*S_u + b_ijk*S_c + d_ijk*I_ghost
    end

    for idy in start_y:-sign_y:stop_y
        y_c = atmos.y[idy]
        # X upwind position
        idy_u = roll(idy+sign_y, ny)
        y_u = atmos.y[idy_u]

        for idx in start_x:-sign_x:stop_x
            # Centre point (c)
            y_c = atmos.y[idy]

            # Lower corner to interpolate from
            idy_i = roll(idy+sign_y, ny)

            y_u = y_c + r*sin(ϕ)*sin(θ)

            # Do linear interpolation to find τ and S

            # exinction at centre and upwind point
            α_c = α_centre[idx, idy]
            α_vals = [α_centre[idx_i, idy_u]   α_centre[idx, idy_u]
                      α_upwind[idx_i, idy_u]   α_upwind[idx, idy_u]]
            α_u = bilinear_xz(z_u, x_u, idz_i, idx_i, atmos, α_vals)

            # Find the Δτ optical path from upwind to grid point
            Δτ_u = trapezoidal(r, α_c, α_u)

            # Find source function at grid point and upwind point (intepolate)
            S_c = S_centre[idx, idy]
            S_vals = [S_centre[idx_i, idy_u]    S_centre[idx, idy_u]
                      S_upwind[idx_i, idy_u]    S_upwind[idx, idy_u]]
            S_u = bilinear_xz(z_u, x_u, idz_i, idx_i, atmos, S_vals)

            # From Hennicker et. al. (2020) equation 13
            e_0 = 1 - exp(-Δτ_u)
            e_1 = Δτ_u - e_0

            # From Hennicker et. al. (2020) equation 12, with ω = 1
            a_ijk = e_0 - e_1/Δτ_u
            b_ijk = e_1/Δτ_u
            # c_ijk = 0
            d_ijk = exp(-Δτ_u)

            # Interpolate to find intensity at upwind point
            I_vals = [I_0[idx_i, idy_u]     I_0[idx, idy_u]
                      I_lower[idx_i]        I_lower[idx]   ]
            I_u = bilinear_xz(z_u, x_u, idz_i, idx_i, atmos, I_vals)

            # Integrate intensity from two-point quadrature
            I[idx, idy] = a_ijk*S_u + b_ijk*S_c + d_ijk*I_u
        end
        # Update upper interpolate value
        I_lower = I[:, idy]
    end
    return I
end
