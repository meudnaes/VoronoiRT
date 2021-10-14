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

    # Allocate space for intensity
    I = zero(I_0)

    # Z center and upwind position
    z_centre = atmos.z[idz]
    idz_upwind = idz - 1
    z_upwind = atmos.z[idz_upwind]

    # Length to plane
    Δz = z_upwind-z_centre
    r = abs(Δz/cos(θ))

    x_increment = r*cos(ϕ)*sin(θ)
    y_increment = r*sin(ϕ)*sin(θ)

    α_lower = α[idz_upwind,:,:]

    S_lower = S_0[idz_upwind,:,:]

    # Loop direction
    start_x, stop_x = range_bounds(sign_x, nx)
    start_y, stop_y = range_bounds(sign_y, ny)

    for idx in start_x:-sign_x:stop_x
        idx_lower = idx + Int((sign_x-1)/2)
        idx_upper = idx_lower + 1

        x_centre = atmos.x[idx]
        x_upwind = x_centre + x_increment

        x_bounds = [atmos.x[idx_lower], atmos.x[idx_upper]]

        for idy in start_y:-sign_y:stop_y
            # Centre point (c)
            y_centre = atmos.y[idy]

            # upwind coordinate
            y_upwind = y_centre + y_increment

            # Lower corner to interpolate from
            idy_lower = idy + Int((sign_y-1)/2)
            idy_upper = idy+1

            y_bounds = [atmos.y[idy_lower], atmos.y[idy_upper]]

            # Do linear interpolation to find τ and S

            # exinction at centre and upwind point
            α_vals = [α_lower[idx_lower, idy_lower] α_lower[idx_lower, idy_upper]
                      α_lower[idx_upper, idy_lower] α_lower[idx_upper, idy_upper]]

            α_centre = α[idz, idx, idy]
            α_upwind = bilinear(x_upwind, y_upwind, x_bounds, y_bounds, α_vals)

            # Find the Δτ optical path from upwind to grid point
            Δτ_upwind = trapezoidal(r, α_centre, α_upwind)

            # Find source function at grid point and upwind point (intepolate)
            S_vals = [S_lower[idx_lower, idy_lower] S_lower[idx_lower, idy_upper]
                      S_lower[idx_upper, idy_lower] S_lower[idx_upper, idy_upper]]

            S_centre = S_0[idz, idx, idy]
            S_upwind = bilinear(x_upwind, y_upwind, x_bounds, y_bounds, S_vals)

            # From Hennicker et. al. (2020) equation 13
            e_0 = 1 - exp(-Δτ_upwind)
            e_1 = Δτ_upwind - e_0

            # From Hennicker et. al. (2020) equation 12, with ω = 1
            a_ijk = e_0 - e_1/Δτ_upwind
            b_ijk = e_1/Δτ_upwind
            # c_ijk = 0
            d_ijk = exp(-Δτ_upwind)

            # Interpolate to find intensity at upwind point
            I_vals = [I_0[idx_lower, idy_lower] I_0[idx_lower, idy_upper]
                      I_0[idx_upper, idy_lower] I_0[idx_upper, idy_upper]]

            I_upwind = bilinear(x_upwind, y_upwind, x_bounds, y_bounds, I_vals)

            # Integrate intensity from two-point quadrature
            I[idx, idy] = a_ijk*S_upwind + b_ijk*S_centre + d_ijk*I_upwind
        end
        # Update ghost zones
        I[idx, 1] = I[idx, end-1]
        I[idx, end] = I[idx, 2]
    end
    # Update ghost zones
    I[1,:] = I[end-1,:]
    I[end,:] = I[2,:]
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

    # nx = length(atmos.x)
    # ny = length(atmos.y)

    # Allocate space for intensity
    I = zero(I_0)
    # Z center and upwind position
    z_c = atmos.z[idz]
    idz_u = idz + 1
    z_u = atmos.z[idz_u]

    # Length to plane
    Δz = abs(z_u - z_c)
    r = Δz/cos(θ)

    x_increment = r*cos(ϕ)*sin(θ)
    y_increment = r*sin(ϕ)*sin(θ)

    α_centre = α[idz,:,:]
    α_upwind = α[idz_u,:,:]

    S_centre = S_0[idz,:,:]
    S_upwind = S_0[idz_u,:,:]

    # Loop direction
    start_x, stop_x = range_bounds(sign_x, nx)
    start_y, stop_y = range_bounds(sign_y, ny)

    for idx in start_x:-sign_x:stop_x
        idx_i = idx + Int((sign_x-1)/2)
        x_c = atmos.x[idx]
        x_u = x_c + x_increment

        x_bounds = [x_c, atmos.x[idx_i]]
        for idy in start_y:-sign_y:stop_y
            # Centre point (c)
            y_c = atmos.y[idy]

            # Lower corner to interpolate from
            idy_i = idy + Int((sign_y-1)/2)

            # calculate the length to intersections with z, x, and y plane
            y_u = y_c + y_increment

            y_bounds = [y_c, atmos.y[idy_i]]

            # Do linear interpolation to find τ and S

            # exinction at centre and upwind point
            α_vals = [α_upwind[idx, idy]   α_upwind[idx, idy_i]
                      α_upwind[idx_i, idy] α_upwind[idx_i, idy_i]]

            α_c = α_centre[idx, idy]
            α_u = bilinear(x_u, y_u, x_bounds, y_bounds, α_vals)

            # Find the Δτ optical path from upwind to grid point
            Δτ_upwind = trapezoidal(abs(r), α_c, α_u)

            # Find source function at grid point and upwind point (intepolate)
            S_vals = [S_upwind[idx, idy]   S_upwind[idx, idy_i]
                      S_upwind[idx_i, idy] S_upwind[idx_i, idy_i]]

            S_c = S_centre[idx, idy]
            S_u = bilinear(x_u, y_u, x_bounds, y_bounds, S_vals)

            # From Hennicker et. al. (2020) equation 13
            e_0 = 1 - exp(-Δτ_upwind)
            e_1 = Δτ_upwind - e_0

            # From Hennicker et. al. (2020) equation 12, with ω = 1
            a_ijk = e_0 - e_1/Δτ_upwind
            b_ijk = e_1/Δτ_upwind
            # c_ijk = 0
            d_ijk = exp(-Δτ_upwind)

            # Interpolate to find intensity at upwind point
            I_vals = [I_0[idx, idy]   I_0[idx, idy_i]
                      I_0[idx_i, idy] I_0[idx_i, idy_i]]

            I_u = bilinear(x_u, y_u, x_bounds, y_bounds, I_vals)

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

    # Loop direction
    start_x, stop_x = range_bounds(sign_x, nx)
    start_y, stop_y = range_bounds(sign_y, ny)

    # Allocate space for intensity
    I = zero(I_0)
    I_upper = zero(I_0[end,:])

    # Z center and interpolation position
    z_centre = atmos.z[idz]
    idz_lower = idz-1

    # calculate the length to upwind position
    r = abs(Δx/(cos(ϕ)*sin(θ)))
    z_increment = r*cos(θ)
    y_increment = r*sin(ϕ)*sin(θ)

    z_upwind = z_centre + z_increment

    z_bounds = [atmos.z[idz_lower], z_centre]

    α_lower = α[idz_lower,:,:]
    α_upper = α[idz,:,:]

    S_lower = S_0[idz_lower,:,:]
    S_upper = S_0[idz,:,:]

    # Do the ghost cell thing...
    idx = stop_x
    x_centre = atmos.x[idx]

    # X upwind position
    idx_upwind = idx+sign_x
    x_upwind = atmos.x[idx_upwind]


    for idy in start_y:-sign_y:stop_y
        # Centre point (c)
        y_centre = atmos.y[idy]

        # Lower corner to interpolate from
        idy_lower = idy+Int((sign_y-1)/2)
        idy_upper = idy_lower+1

        y_upwind = y_centre + y_increment

        y_bounds = [atmos.y[idy_lower], atmos.y[idy_upper]]

        # Do linear interpolation to find τ and S

        # exinction at centre and upwind point
        # z is 1st coordinate, changes rows in vals. y is 2nd coordinate
        α_vals = [α_lower[idx_upwind, idy_lower]   α_lower[idx_upwind, idy_upper]
                  α_upper[idx_upwind, idy_lower]   α_upper[idx_upwind, idy_upper]]

        α_centre = α_upper[idx, idy]
        α_upwind = bilinear(z_upwind, y_upwind, z_bounds, y_bounds, α_vals)

        # Find the Δτ optical path from upwind to grid point
        Δτ_upwind = trapezoidal(r, α_centre, α_upwind)

        # Find source function at grid point and upwind point (intepolate)
        S_vals = [S_lower[idx_upwind, idy_lower]    S_lower[idx_upwind, idy_upper]
                  S_upper[idx_upwind, idy_lower]    S_upper[idx_upwind, idy_upper]]

        S_centre = S_upper[idx, idy]
        S_upwind = bilinear(z_upwind, y_upwind, z_bounds, y_bounds, S_vals)

        # From Hennicker et. al. (2020) equation 13, also see Hauschildt et. al. p. 276
        e_0 = 1 - exp(-Δτ_upwind)
        e_1 = Δτ_upwind - e_0

        # From Hennicker et. al. (2020) equation 12, with ω = 1
        a_ijk = e_0 - e_1/Δτ_upwind
        b_ijk = e_1/Δτ_upwind
        # c_ijk = 0
        d_ijk = exp(-Δτ_upwind)

        # Interpolate to find intensity at upwind point
        # I_0 is from the top
        I_vals = [I_0[idx_upwind, idy_lower]     I_0[idx_upwind, idy_upper]
                  S_upper[idx_upwind, idy_lower] S_upper[idx_upwind, idy_upper]]


        I_upwind = bilinear(z_upwind, y_upwind, z_bounds, y_bounds, I_vals)

        # Integrate intensity from two-point quadrature
        I_upper[idy_lower] = a_ijk*S_upwind + b_ijk*S_centre + d_ijk*I_upwind
    end

    I_upper[1] = I_upper[end-1]
    I_upper[end] = I_upper[2]

    for idx in start_x:-sign_x:stop_x
        x_centre = atmos.x[idx]
        # x upwind position
        idx_upwind = idx+sign_x
        x_upwind = atmos.x[idx_upwind]
        for idy in start_y:-sign_y:stop_y
            # Centre point (c)
            y_centre = atmos.y[idy]

            # Lower corner to interpolate from
            idy_lower = idy+Int((sign_y-1)/2)
            idy_upper = idy_lower+1

            y_upwind = y_centre + y_increment

            y_bounds = [atmos.y[idy_lower], atmos.y[idy_upper]]

            # Do linear interpolation to find τ and S

            # exinction at centre and upwind point
            # z is 1st coordinate, changes rows in vals. y is 2nd coordinate
            α_vals = [α_lower[idx_upwind, idy_lower]   α_lower[idx_upwind, idy_upper]
                      α_upper[idx_upwind, idy_lower]   α_upper[idx_upwind, idy_upper]]

            α_centre = α_upper[idx, idy]
            α_upwind = bilinear(z_upwind, y_upwind, z_bounds, y_bounds, α_vals)

            # Find the Δτ optical path from upwind to grid point
            Δτ_upwind = trapezoidal(r, α_centre, α_upwind)

            # Find source function at grid point and upwind point (intepolate)
            S_vals = [S_lower[idx_upwind, idy_lower]    S_lower[idx_upwind, idy_upper]
                      S_upper[idx_upwind, idy_lower]    S_upper[idx_upwind, idy_upper]]

            S_centre = S_upper[idx, idy]
            S_upwind = bilinear(z_upwind, y_upwind, z_bounds, y_bounds, S_vals)

            # From Hennicker et. al. (2020) equation 13, also see Hauschildt et. al. p. 276
            e_0 = 1 - exp(-Δτ_upwind)
            e_1 = Δτ_upwind - e_0

            # From Hennicker et. al. (2020) equation 12, with ω = 1
            a_ijk = e_0 - e_1/Δτ_upwind
            b_ijk = e_1/Δτ_upwind
            # c_ijk = 0
            d_ijk = exp(-Δτ_upwind)

            # Interpolate to find intensity at upwind point
            # I_0 is from the top
            I_vals = [I_0[idx_upwind, idy_lower]    I_0[idx_upwind, idy_upper]
                      I_upper[idy_lower]            I_upper[idy_upper]        ]


            I_upwind = bilinear(z_upwind, y_upwind, z_bounds, y_bounds, I_vals)

            # Integrate intensity from two-point quadrature
            I[idx, idy] = a_ijk*S_upwind + b_ijk*S_centre + d_ijk*I_upwind
        end
        # Update ghost zones
        I[idx, 1] = I[idx, end-1]
        I[idx, end] = I[idx, 2]

        # Update upper interpolate value
        I_upper = I[idx, :]
    end

    # Update ghost zones
    I[1,:] = I[end-1,:]
    I[end,:] = I[2,:]

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

    # Loop direction
    start_x, stop_x = range_bounds(sign_x, nx)
    start_y, stop_y = range_bounds(sign_y, ny)

    # Allocate space for intensity
    I = zero(I_0)
    I_lower = zero(I_0[end,:])

    # Z center and interpolation position
    z_centre = atmos.z[idz]
    idz_upper = idz+1

    # Length to plane
    # Δx = atmos.x[2] - atmos.x[1]

    # calculate the length to upwind position
    r = abs(Δx/(cos(ϕ)*sin(θ)))
    z_increment = r*cos(θ)
    y_increment = r*sin(ϕ)*sin(θ)

    z_upwind = z_centre + z_increment

    z_bounds = [z_centre, atmos.z[idz_upper]]

    α_lower = α[idz,:,:]
    α_upper = α[idz_upper,:,:]

    S_lower = S_0[idz,:,:]
    S_upper = S_0[idz_upper,:,:]

    # Do the ghost cell thing...
    idx = stop_x
    x_centre = atmos.x[idx]

    # X upwind position
    idx_upwind = idx+sign_x
    x_upwind = atmos.x[idx_upwind]

    for idy in start_y:-sign_y:stop_y
        # Centre point (c)
        y_centre = atmos.y[idy]

        # Lower corner to interpolate from
        idy_lower = idy+Int((sign_y-1)/2)
        idy_upper = idy_lower+1

        y_upwind = y_centre + y_increment

        y_bounds = [atmos.y[idy_lower], atmos.y[idy_upper]]

        # Do linear interpolation to find τ and S

        # exinction at centre and upwind point
        # z is 1st coordinate, changes rows in vals. y is 2nd coordinate
        α_vals = [α_lower[idx_upwind, idy_lower]   α_lower[idx_upwind, idy_upper]
                  α_upper[idx_upwind, idy_lower]   α_upper[idx_upwind, idy_upper]]

        α_centre = α_lower[idx, idy]
        α_upwind = bilinear(z_upwind, y_upwind, z_bounds, y_bounds, α_vals)

        # Find the Δτ optical path from upwind to grid point
        Δτ_upwind = trapezoidal(r, α_centre, α_upwind)

        # Find source function at grid point and upwind point (intepolate)
        S_vals = [S_lower[idx_upwind, idy_lower]    S_lower[idx_upwind, idy_upper]
                  S_upper[idx_upwind, idy_lower]    S_upper[idx_upwind, idy_upper]]

        S_centre = S_lower[idx, idy]
        S_upwind = bilinear(z_upwind, y_upwind, z_bounds, y_bounds, S_vals)

        # From Hennicker et. al. (2020) equation 13, also see Hauschildt et. al. p. 276
        e_0 = 1 - exp(-Δτ_upwind)
        e_1 = Δτ_upwind - e_0

        # From Hennicker et. al. (2020) equation 12, with ω = 1
        a_ijk = e_0 - e_1/Δτ_upwind
        b_ijk = e_1/Δτ_upwind
        # c_ijk = 0
        d_ijk = exp(-Δτ_upwind)

        # Interpolate to find intensity at upwind point
        # I_0 is from the top
        I_vals = [S_lower[idx_upwind, idy_lower]    S_lower[idx_upwind, idy_upper]
                  I_0[idx_upwind, idy_lower]        I_0[idx_upwind, idy_upper]    ]


        I_upwind = bilinear(z_upwind, y_upwind, z_bounds, y_bounds, I_vals)

        # Integrate intensity from two-point quadrature
        I_lower[idy_lower] = a_ijk*S_upwind + b_ijk*S_centre + d_ijk*I_upwind
    end

    I_lower[1] = I_lower[end-1]
    I_lower[end] = I_lower[2]

    for idx in start_x:-sign_x:stop_x
        x_centre = atmos.x[idx]
        # X upwind position
        idx_upwind = idx+sign_x
        x_upwind = atmos.x[idx_upwind]

        for idy in start_y:-sign_y:stop_y
            # Centre point (c)
            y_centre = atmos.y[idy]

            # Lower corner to interpolate from
            idy_lower = idy+Int((sign_y-1)/2)
            idy_upper = idy_lower+1

            y_upwind = y_centre + y_increment

            y_bounds = [atmos.y[idy_lower], atmos.y[idy_upper]]

            # Do linear interpolation to find τ and S

            # exinction at centre and upwind point
            # z is 1st coordinate, changes rows in vals. y is 2nd coordinate
            α_vals = [α_lower[idx_upwind, idy_lower]   α_lower[idx_upwind, idy_upper]
                      α_upper[idx_upwind, idy_lower]   α_upper[idx_upwind, idy_upper]]

            α_centre = α_lower[idx, idy]
            α_upwind = bilinear(z_upwind, y_upwind, z_bounds, y_bounds, α_vals)

            # Find the Δτ optical path from upwind to grid point
            Δτ_upwind = trapezoidal(r, α_centre, α_upwind)

            # Find source function at grid point and upwind point (intepolate)
            S_vals = [S_lower[idx_upwind, idy_lower]    S_lower[idx_upwind, idy_upper]
                      S_upper[idx_upwind, idy_lower]    S_upper[idx_upwind, idy_upper]]

            S_centre = S_lower[idx, idy]
            S_upwind = bilinear(z_upwind, y_upwind, z_bounds, y_bounds, S_vals)

            # From Hennicker et. al. (2020) equation 13
            e_0 = 1 - exp(-Δτ_upwind)
            e_1 = Δτ_upwind - e_0

            # From Hennicker et. al. (2020) equation 12, with ω = 1
            a_ijk = e_0 - e_1/Δτ_upwind
            b_ijk = e_1/Δτ_upwind
            # c_ijk = 0
            d_ijk = exp(-Δτ_upwind)

            # Interpolate to find intensity at upwind point
            # I_0 is from the top
            I_vals = [I_lower[idy_lower]         I_lower[idy_upper]
                      I_0[idx_upwind, idy_lower] I_0[idx_upwind, idy_upper]]

            I_upwind = bilinear(z_upwind, y_upwind, z_bounds, y_bounds, I_vals)

            # Integrate intensity from two-point quadrature
            I[idx, idy] = a_ijk*S_upwind + b_ijk*S_centre + d_ijk*I_upwind
        end

        # Update ghost zones
        I[idx, 1] = I[idx, end-1]
        I[idx, end] = I[idx, 2]

        # Update upper interpolate value
        I_lower = I[idx, :]
    end

    # Update ghost zones
    I[1,:] = I[end-1,:]
    I[end,:] = I[2,:]

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

    # Loop direction
    start_x, stop_x = range_bounds(sign_x, nx)
    start_y, stop_y = range_bounds(sign_y, ny)

    # Allocate space for intensity
    I = zero(I_0)
    I_upper = zero(I_0[end,:])

    # Z center and interpolation position
    z_centre = atmos.z[idz]
    idz_lower = idz-1

    # calculate the length to upwind position
    r = abs(Δy/(sin(ϕ)*sin(θ)))
    z_increment = r*cos(θ)
    x_increment = r*cos(ϕ)*sin(θ)

    z_upwind = z_centre + z_increment

    z_bounds = [atmos.z[idz_lower], z_centre]

    α_lower = α[idz_lower,:,:]
    α_upper = α[idz,:,:]

    S_lower = S_0[idz_lower,:,:]
    S_upper = S_0[idz,:,:]

    # Do the ghost cell thing...
    idy = stop_y
    y_centre = atmos.y[idy]

    # X upwind position
    idy_upwind = idy+sign_y
    y_upwind = atmos.y[idy_upwind]


    for idx in start_x:-sign_x:stop_x
        # Centre point (c)
        x_centre = atmos.x[idx]

        # Lower corner to interpolate from
        idx_lower = idx+Int((sign_x-1)/2)
        idx_upper = idx_lower+1

        x_upwind = x_centre + x_increment

        x_bounds = [atmos.x[idx_lower], atmos.x[idx_upper]]

        # Do linear interpolation to find τ and S

        # exinction at centre and upwind point
        # z is 1st coordinate, changes rows in vals. y is 2nd coordinate
        α_vals = [α_lower[idx_lower, idy_upwind]   α_lower[idx_upper, idy_upwind]
                  α_upper[idx_lower, idy_upwind]   α_upper[idx_upper, idy_upwind]]

        α_centre = α_upper[idx, idy]
        α_upwind = bilinear(z_upwind, x_upwind, z_bounds, x_bounds, α_vals)

        # Find the Δτ optical path from upwind to grid point
        Δτ_upwind = trapezoidal(r, α_centre, α_upwind)

        # Find source function at grid point and upwind point (intepolate)
        S_vals = [S_lower[idx_lower, idy_upwind]    S_lower[idx_upper, idy_upwind]
                  S_upper[idx_lower, idy_upwind]    S_upper[idx_upper, idy_upwind]]

        S_centre = S_upper[idx, idy]
        S_upwind = bilinear(z_upwind, x_upwind, z_bounds, x_bounds, S_vals)

        # From Hennicker et. al. (2020) equation 13, also see Hauschildt et. al. p. 276
        e_0 = 1 - exp(-Δτ_upwind)
        e_1 = Δτ_upwind - e_0

        # From Hennicker et. al. (2020) equation 12, with ω = 1
        a_ijk = e_0 - e_1/Δτ_upwind
        b_ijk = e_1/Δτ_upwind
        # c_ijk = 0
        d_ijk = exp(-Δτ_upwind)

        # Interpolate to find intensity at upwind point
        # I_0 is from the top
        I_vals = [I_0[idx_lower, idy_upwind]     I_0[idx_upper, idy_upwind]
                  S_upper[idx_lower, idy_upwind] S_upper[idx_upper, idy_upwind]]


        I_upwind = bilinear(z_upwind, x_upwind, z_bounds, x_bounds, I_vals)

        # Integrate intensity from two-point quadrature
        I_upper[idx_lower] = a_ijk*S_upwind + b_ijk*S_centre + d_ijk*I_upwind
    end

    I_upper[1] = I_upper[end-1]
    I_upper[end] = I_upper[2]

    for idy in start_y:-sign_y:stop_y
        y_centre = atmos.y[idy]
        # x upwind position
        idy_upwind = idy+sign_y
        y_upwind = atmos.y[idy_upwind]
        for idx in start_x:-sign_x:stop_x
            # Centre point (c)
            x_centre = atmos.x[idx]

            # Lower corner to interpolate from
            idx_lower = idx+Int((sign_x-1)/2)
            idx_upper = idx_lower+1

            x_upwind = x_centre + x_increment

            x_bounds = [atmos.x[idx_lower], atmos.x[idx_upper]]

            # Do linear interpolation to find τ and S

            # exinction at centre and upwind point
            # z is 1st coordinate, changes rows in vals. y is 2nd coordinate
            α_vals = [α_lower[idx_lower, idy_upwind]   α_lower[idx_upper, idy_upwind]
                      α_upper[idx_lower, idy_upwind]   α_upper[idx_upper, idy_upwind]]

            α_centre = α_upper[idx, idy]
            α_upwind = bilinear(z_upwind, x_upwind, z_bounds, x_bounds, α_vals)

            # Find the Δτ optical path from upwind to grid point
            Δτ_upwind = trapezoidal(r, α_centre, α_upwind)

            # Find source function at grid point and upwind point (intepolate)
            S_vals = [S_lower[idx_lower, idy_upwind]    S_lower[idx_upper, idy_upwind]
                      S_upper[idx_lower, idy_upwind]    S_upper[idx_upper, idy_upwind]]

            S_centre = S_upper[idx, idy]
            S_upwind = bilinear(z_upwind, x_upwind, z_bounds, x_bounds, S_vals)

            # From Hennicker et. al. (2020) equation 13, also see Hauschildt et. al. p. 276
            e_0 = 1 - exp(-Δτ_upwind)
            e_1 = Δτ_upwind - e_0

            # From Hennicker et. al. (2020) equation 12, with ω = 1
            a_ijk = e_0 - e_1/Δτ_upwind
            b_ijk = e_1/Δτ_upwind
            # c_ijk = 0
            d_ijk = exp(-Δτ_upwind)

            # Interpolate to find intensity at upwind point
            # I_0 is from the top
            I_vals = [I_0[idx_lower, idy_upwind]    I_0[idx_upper, idy_upwind]
                      I_upper[idx_lower]            I_upper[idx_upper]        ]


            I_upwind = bilinear(z_upwind, x_upwind, z_bounds, x_bounds, I_vals)

            # Integrate intensity from two-point quadrature
            I[idx, idy] = a_ijk*S_upwind + b_ijk*S_centre + d_ijk*I_upwind
        end
        # Update ghost zones
        I[1, idy] = I[end-1, idy]
        I[end, idy] = I[2, idy]

        # Update upper interpolate value
        I_upper = I[:, idy]
    end

    # Update ghost zones
    I[:, 1] = I[:, end-1]
    I[:, end] = I[:, 2]

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

    # Loop direction
    start_x, stop_x = range_bounds(sign_x, nx)
    start_y, stop_y = range_bounds(sign_y, ny)

    # Allocate space for intensity
    I = zero(I_0)
    I_lower = zero(I_0[end,:])

    # Z center and interpolation position
    z_centre = atmos.z[idz]
    idz_upper = idz+1

    # Length to plane
    # Δx = atmos.x[2] - atmos.x[1]

    # calculate the length to upwind position
    r = abs(Δx/(sin(ϕ)*sin(θ)))
    z_increment = r*cos(θ)
    x_increment = r*cos(ϕ)*sin(θ)

    z_upwind = z_centre + z_increment

    z_bounds = [z_centre, atmos.z[idz_upper]]

    α_lower = α[idz,:,:]
    α_upper = α[idz_upper,:,:]

    S_lower = S_0[idz,:,:]
    S_upper = S_0[idz_upper,:,:]

    # Do the ghost cell thing...
    idy = stop_y
    y_centre = atmos.y[idy]

    # X upwind position
    idy_upwind = idy+sign_y
    y_upwind = atmos.y[idy_upwind]

    for idx in start_x:-sign_x:stop_x
        # Centre point (c)
        x_centre = atmos.x[idx]

        # Lower corner to interpolate from
        idx_lower = idx+Int((sign_x-1)/2)
        idx_upper = idx_lower+1

        x_upwind = x_centre + x_increment

        x_bounds = [atmos.x[idx_lower], atmos.x[idx_upper]]

        # Do linear interpolation to find τ and S

        # exinction at centre and upwind point
        # z is 1st coordinate, changes rows in vals. y is 2nd coordinate
        α_vals = [α_lower[idx_lower, idy_upwind]   α_lower[idx_upper, idy_upwind]
                  α_upper[idx_lower, idy_upwind]   α_upper[idx_upper, idy_upwind]]

        α_centre = α_lower[idx, idy]
        α_upwind = bilinear(z_upwind, x_upwind, z_bounds, x_bounds, α_vals)

        # Find the Δτ optical path from upwind to grid point
        Δτ_upwind = trapezoidal(r, α_centre, α_upwind)

        # Find source function at grid point and upwind point (intepolate)
        S_vals = [S_lower[idx_lower, idy_upwind]    S_lower[idx_upper, idy_upwind]
                  S_upper[idx_lower, idy_upwind]    S_upper[idx_upper, idy_upwind]]

        S_centre = S_lower[idx, idy]
        S_upwind = bilinear(z_upwind, x_upwind, z_bounds, x_bounds, S_vals)

        # From Hennicker et. al. (2020) equation 13, also see Hauschildt et. al. p. 276
        e_0 = 1 - exp(-Δτ_upwind)
        e_1 = Δτ_upwind - e_0

        # From Hennicker et. al. (2020) equation 12, with ω = 1
        a_ijk = e_0 - e_1/Δτ_upwind
        b_ijk = e_1/Δτ_upwind
        # c_ijk = 0
        d_ijk = exp(-Δτ_upwind)

        # Interpolate to find intensity at upwind point
        # I_0 is from the top
        I_vals = [S_lower[idx_lower, idy_upwind]    S_lower[idx_upper, idy_upwind]
                  I_0[idx_lower, idy_upwind]        I_0[idx_upper, idy_upwind]    ]

        I_upwind = bilinear(z_upwind, x_upwind, z_bounds, x_bounds, I_vals)

        # Integrate intensity from two-point quadrature
        I_lower[idx_lower] = a_ijk*S_upwind + b_ijk*S_centre + d_ijk*I_upwind
    end

    I_lower[1] = I_lower[end-1]
    I_lower[end] = I_lower[2]

    for idy in start_y:-sign_y:stop_y
        y_centre = atmos.y[idy]
        # x upwind position
        idy_upwind = idy+sign_y
        y_upwind = atmos.y[idy_upwind]
        for idx in start_x:-sign_x:stop_x
            # Centre point (c)
            x_centre = atmos.x[idx]

            # Lower corner to interpolate from
            idx_lower = idx+Int((sign_x-1)/2)
            idx_upper = idx_lower+1

            x_upwind = x_centre + x_increment

            x_bounds = [atmos.x[idx_lower], atmos.x[idx_upper]]

            # Do linear interpolation to find τ and S

            # exinction at centre and upwind point
            # z is 1st coordinate, changes rows in vals. y is 2nd coordinate
            α_vals = [α_lower[idx_lower, idy_upwind]   α_lower[idx_upper, idy_upwind]
                      α_upper[idx_lower, idy_upwind]   α_upper[idx_upper, idy_upwind]]

            α_centre = α_upper[idx, idy]
            α_upwind = bilinear(z_upwind, x_upwind, z_bounds, x_bounds, α_vals)

            # Find the Δτ optical path from upwind to grid point
            Δτ_upwind = trapezoidal(r, α_centre, α_upwind)

            # Find source function at grid point and upwind point (intepolate)
            S_vals = [S_lower[idx_lower, idy_upwind]    S_lower[idx_upper, idy_upwind]
                      S_upper[idx_lower, idy_upwind]    S_upper[idx_upper, idy_upwind]]

            S_centre = S_upper[idx, idy]
            S_upwind = bilinear(z_upwind, x_upwind, z_bounds, x_bounds, S_vals)

            # From Hennicker et. al. (2020) equation 13, also see Hauschildt et. al. p. 276
            e_0 = 1 - exp(-Δτ_upwind)
            e_1 = Δτ_upwind - e_0

            # From Hennicker et. al. (2020) equation 12, with ω = 1
            a_ijk = e_0 - e_1/Δτ_upwind
            b_ijk = e_1/Δτ_upwind
            # c_ijk = 0
            d_ijk = exp(-Δτ_upwind)

            # Interpolate to find intensity at upwind point
            # I_0 is from the top
            I_vals = [I_lower[idx_lower]            I_lower[idx_upper]
                      I_0[idx_lower, idy_upwind]    I_0[idx_upper, idy_upwind]]

            I_upwind = bilinear(z_upwind, x_upwind, z_bounds, x_bounds, I_vals)

            # Integrate intensity from two-point quadrature
            I[idx, idy] = a_ijk*S_upwind + b_ijk*S_centre + d_ijk*I_upwind
        end

        # Update ghost zones
        I[1, idy] = I[end-1, idy]
        I[end, idy] = I[2, idy]

        # Update upper interpolate value
        I_lower = I[:, idy]
    end

    # Update ghost zones
    I[:, 1] = I[:, end-1]
    I[:, end] = I[:, 2]

    return I
end
