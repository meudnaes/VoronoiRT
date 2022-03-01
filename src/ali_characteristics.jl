#=
Short characteristic method to solve intensity along rays going through all grid
points.
=#

include("functions.jl")

"""
    function short_characteristics_up(k::Vector{Float64},
                                      α::Array{<:PerLength, 3},
                                      atmos::Atmosphere,
                                      Λstar::Matrix{Float64})

Computes intensity along rays traveling upwards from the bottom to the top
through all grid points in the atmosphere. Initial intensity I_0 = B_λ(T) at the
bottom of the domain.
"""
function short_characteristics_up(k::Vector{Float64},
                                  α::Array{<:PerLength, 3},
                                  atmos::Atmosphere,
                                  Λstar::Matrix{Float64})
    ############################################################################
    # | and - : Grid
    #     x   : Grid points
    #     *   : Ray
    #
    #                  x----------x----------x  z_i+1
    #                  |          |          |
    #                  |          |          |
    #                  |          |          |
    #                  x----------*----------x  z_i
    #                  |          |*         |
    #                  |          | *        |
    #                  |          |  *       |
    #                  x----------x---*------x  z_i-1
    #                                 :
    #                             upwind point
    ############################################################################

    nz = length(atmos.z)
    Δx = atmos.x[2] - atmos.x[1]
    Δy = atmos.y[2] - atmos.y[1]

    # I know that Δx = Δy = constant for all grid points
    # Δxy = atmos.x[2] - atmos.x[1]
    r_x = abs(Δx/k[2])
    r_y = abs(Δy/k[3])

    # find out direction in xy the ray moves (1 for positive, -1 for negative)
    sign_x, sign_y = xy_intersect(k)

    # loop upwards through atmosphere
    for idz in 2:nz
        # Grid diffference in z
        Δz = atmos.z[idz] - atmos.z[idz-1]

        # calculate length until ray hits z plane
        r_z = abs(Δz/k[1])

        # This finds which plane the ray intersects with
        plane_cut = argmin([r_z, r_x, r_y])
        if plane_cut == 1
            xy_up_ray(k, idz, sign_x, sign_y, α, atmos)
        elseif plane_cut==2
            yz_up_ray(k, idz, sign_x, sign_y, α, atmos)
        elseif plane_cut==3
            xz_up_ray(k, idz, sign_x, sign_y, α, atmos)
        end
    end

    return I
end

"""
    function short_characteristics_down(k::Vector{Float64},
                                        α::Array{<:PerLength, 3},
                                        atmos::Atmosphere,
                                        Λstar::Matrix{Float64})

Computes intensity along rays traveling downwards from the top to the bottom
through all grid points in the atmosphere. Initial intensity I_0 = 0 at the top
of the domain.
"""
function short_characteristics_down(k::Vector{Float64},
                                    α::Array{<:PerLength, 3},
                                    atmos::Atmosphere,
                                    Λstar::Matrix{Float64})
    ############################################################################
    # | and - : Grid
    #     x   : Grid points
    #     *   : Ray
    #
    #                  x----------x----------x  z_i+1
    #                  |          |          *  upwind point
    #                  |          |      *   |
    #                  |          |  *       |
    #                  x----------*----------x  z_i
    #                  |          |          |
    #                  |          |          |
    #                  |          |          |
    #                  x----------x----------x  z_i-1
    #
    ############################################################################

    nz = length(atmos.z)
    Δx = atmos.x[2] - atmos.x[1]
    Δy = atmos.y[2] - atmos.y[1]

    # I know that Δx = Δy = constant for all grid points
    r_x = abs(Δx/k[2])
    r_y = abs(Δy/k[3])

    # find out which plane the upwind part of the ray intersects
    sign_x, sign_y = xy_intersect(k)

    # loop downwards through atmosphere
    for idz in nz-1:-1:1
        # Grid diffference in z
        Δz = atmos.z[idz+1] - atmos.z[idz]

        # calculate length until ray hits z plane
        r_z = abs(Δz/k[1])

        # This finds which plane the ray intersects with
        plane_cut = argmin([r_z, r_x, r_y])
        if plane_cut == 1
            xy_down_ray(k, idz, sign_x, sign_y, α, atmos)
        elseif plane_cut==2
            yz_down_ray(k, idz, sign_x, sign_y, α, atmos)
        elseif plane_cut==3
            xz_down_ray(k, idz, sign_x, sign_y, α, atmos)
        end
    end

    return I
end

"""
    function xy_up_ray(θ::AbstractFloat, ϕ::AbstractFloat, idz::Int, sign_x::Int,
                       sign_y::Int, I_0::Matrix{<:UnitsIntensity_λ}, S_0::Array{<:UnitsIntensity_λ, 3},
                       α::Array{<:PerLength, 3}, atmos::Atmosphere)

Ray moving upwards, upwind point intersecting with lower xy plane. Assumes angle
in radians.
"""
function xy_up_ray(k::Vector{Float64}, idz::Int, sign_x::Int, sign_y::Int,
                   α::Array{<:PerLength, 3}, atmos::Atmosphere)

    nx = length(atmos.x)
    ny = length(atmos.y)

    # Z center and upwind position
    z_centre = atmos.z[idz]
    idz_upwind = idz - 1
    z_upwind = atmos.z[idz_upwind]

    # Length to plane
    Δz = z_upwind-z_centre
    r = abs(Δz/k[1])

    x_increment = r*k[2]
    y_increment = r*k[3]

    α_lower = α[idz_upwind,:,:]

    for idx in 2:nx-1
        idx_lower = idx - Int((sign_x+1)/2)
        idx_upper = idx_lower + 1

        x_centre = atmos.x[idx]
        x_upwind = x_centre + x_increment

        x_bounds = (atmos.x[idx_lower], atmos.x[idx_upper])

        for idy in 2:ny-1
            # Centre point (c)
            y_centre = atmos.y[idy]

            # upwind coordinate
            y_upwind = y_centre + y_increment

            # Lower corner to interpolate from
            idy_lower = idy - Int((sign_y+1)/2)
            idy_upper = idy_lower + 1

            y_bounds = (atmos.y[idy_lower], atmos.y[idy_upper])

            # Do linear interpolation to find τ and S

            # exinction at centre and upwind point
            α_vals = [α_lower[idx_lower, idy_lower] α_lower[idx_lower, idy_upper]
                      α_lower[idx_upper, idy_lower] α_lower[idx_upper, idy_upper]]

            α_centre = α[idz, idx, idy]
            α_upwind = bilinear(x_upwind, y_upwind, x_bounds, y_bounds, α_vals)

            # Find the Δτ optical path from upwind to grid point
            Δτ_upwind = trapezoidal(r, α_centre, α_upwind)

            w1, w2 =  weights(Δτ_upwind)
            # Fill approximate Λ operator
        end

        # Update ghost zones
        # I[idx, 1] = I[idx, end-1]
        # I[idx, end] = I[idx, 2]
    end
    # Update ghost zones
    # I[1,:] = I[end-1,:]
    # I[end,:] = I[2,:]
end

"""
    function xy_down_ray(θ::AbstractFloat, ϕ::AbstractFloat, idz::Int, sign_x::Int,
                         sign_y::Int, I_0::Matrix{<:UnitsIntensity_λ}, S_0::Array{<:UnitsIntensity_λ, 3},
                         α::Array{<:PerLength, 3}, atmos::Atmosphere)

Ray moving downwards, upwind point intersecting with upper xy plane. Assumes
angle in radians.
"""
function xy_down_ray(k::Vector{Float64}, idz::Int, sign_x::Int, sign_y::Int,
                     α::Array{<:PerLength, 3}, atmos::Atmosphere)

    nx = length(atmos.x)
    ny = length(atmos.y)

    # Z center and upwind position
    z_centre = atmos.z[idz]
    idz_upwind = idz + 1
    z_upwind = atmos.z[idz_upwind]

    # Length to plane
    Δz = z_upwind - z_centre
    r = abs(Δz/k[1])

    x_increment = r*k[2]
    y_increment = r*k[3]

    α_upper = α[idz_upwind,:,:]

    for idx in 2:nx-1
        idx_lower = idx - Int((sign_x+1)/2)
        idx_upper = idx_lower+1

        x_centre = atmos.x[idx]
        x_upwind = x_centre + x_increment

        x_bounds = (atmos.x[idx_lower], atmos.x[idx_upper])

        for idy in 2:ny-1
            idy_lower = idy - Int((sign_y+1)/2)
            idy_upper = idy_lower+1

            y_centre = atmos.y[idy]
            y_upwind = y_centre + y_increment

            y_bounds = (atmos.y[idy_lower], atmos.y[idy_upper])
            # Do linear interpolation to find τ and S

            # exinction at centre and upwind point
            α_vals = [α_upper[idx_lower, idy_lower] α_upper[idx_lower, idy_upper]
                      α_upper[idx_upper, idy_lower] α_upper[idx_upper, idy_upper]]

            α_centre = α[idz, idx, idy]
            α_upwind = bilinear(x_upwind, y_upwind, x_bounds, y_bounds, α_vals)

            # Find the Δτ optical path from upwind to grid point
            Δτ_upwind = trapezoidal(r, α_centre, α_upwind)

            w1, w2 =  weights(Δτ_upwind)
        end
        # Update ghost zones
        # I[idx, 1] = I[idx, end-1]
        # I[idx, end] = I[idx, 2]
    end
    # Update ghost zones
    # I[1,:] = I[end-1,:]
    # I[end,:] = I[2,:]
end

"""
    function yz_up_ray(θ::AbstractFloat, ϕ::AbstractFloat, idz::Int, sign_x::Int,
                       sign_y::Int, I_0::Matrix{<:UnitsIntensity_λ}, S_0::Array{<:UnitsIntensity_λ, 3},
                       α::Array{<:PerLength, 3}, atmos::Atmosphere)

Ray moving upwards, upwind point intersecting with yz plane. Assumes angle
in radians.
"""
function yz_up_ray(k::Vector{Float64}, idz::Int, sign_x::Int, sign_y::Int,
                   α::Array{<:PerLength, 3}, atmos::Atmosphere)

    nx = length(atmos.x)
    ny = length(atmos.y)

    Δx = atmos.x[2] - atmos.x[1]

    # Loop direction
    start_x, stop_x = range_bounds(sign_x, nx)
    start_y, stop_y = range_bounds(sign_y, ny)

    # Z center and interpolation position
    z_centre = atmos.z[idz]
    idz_lower = idz-1

    # calculate the length to upwind position
    r = abs(Δx/k[2])
    z_increment = r*k[1]
    y_increment = r*k[3]

    z_upwind = z_centre + z_increment

    z_bounds = (atmos.z[idz_lower], z_centre)

    α_lower = α[idz_lower,:,:]
    α_upper = α[idz,:,:]

    for idx in start_x:sign_x:stop_x
        x_centre = atmos.x[idx]
        # x upwind position
        idx_upwind = idx+sign_x
        x_upwind = atmos.x[idx_upwind]
        for idy in start_y:sign_y:stop_y
            # Centre point (c)
            y_centre = atmos.y[idy]

            # Lower corner to interpolate from
            idy_lower = idy-Int((sign_y+1)/2)
            idy_upper = idy_lower+1

            y_upwind = y_centre + y_increment

            y_bounds = (atmos.y[idy_lower], atmos.y[idy_upper])

            # Do linear interpolation to find τ and S

            # exinction at centre and upwind point
            # z is 1st coordinate, changes rows in vals. y is 2nd coordinate
            α_vals = [α_lower[idx_upwind, idy_lower]   α_lower[idx_upwind, idy_upper]
                      α_upper[idx_upwind, idy_lower]   α_upper[idx_upwind, idy_upper]]

            α_centre = α_upper[idx, idy]
            α_upwind = bilinear(z_upwind, y_upwind, z_bounds, y_bounds, α_vals)

            # Find the Δτ optical path from upwind to grid point
            Δτ_upwind = trapezoidal(r, α_centre, α_upwind)

            w1, w2 =  weights(Δτ_upwind)
        end
        # Update ghost zones
        # I[idx, 1] = I[idx, end-1]
        # I[idx, end] = I[idx, 2]

    end

    # Update ghost zones
    # I[1,:] = I[end-1,:]
    # I[end,:] = I[2,:]
end

"""
    function yz_down_ray(k::Vector{Float64}, idz::Int, sign_x::Int, sign_y::Int,
                       α::Array{<:PerLength, 3}, atmos::Atmosphere)

Ray moving downwards, upwind point intersecting with yz plane. Assumes angle
in radians
"""
function yz_down_ray(k::Vector{Float64}, idz::Int, sign_x::Int, sign_y::Int,
                     α::Array{<:PerLength, 3}, atmos::Atmosphere)

    nx = length(atmos.x)
    ny = length(atmos.y)

    Δx = atmos.x[2] - atmos.x[1]

    # Loop direction
    start_x, stop_x = range_bounds(sign_x, nx)
    start_y, stop_y = range_bounds(sign_y, ny)

    # Z center and interpolation position
    z_centre = atmos.z[idz]
    idz_upper = idz+1

    # Length to plane
    # Δx = atmos.x[2] - atmos.x[1]

    # calculate the length to upwind position
    r = abs(Δx/k[2])
    z_increment = r*k[1]
    y_increment = r*k[3]

    z_upwind = z_centre + z_increment

    z_bounds = (z_centre, atmos.z[idz_upper])

    α_lower = α[idz,:,:]
    α_upper = α[idz_upper,:,:]

    for idx in start_x:sign_x:stop_x
        x_centre = atmos.x[idx]
        # X upwind position
        idx_upwind = idx+sign_x
        x_upwind = atmos.x[idx_upwind]

        for idy in start_y:sign_y:stop_y
            # Centre point (c)
            y_centre = atmos.y[idy]

            # Lower corner to interpolate from
            idy_lower = idy-Int((sign_y+1)/2)
            idy_upper = idy_lower+1

            y_upwind = y_centre + y_increment

            y_bounds = (atmos.y[idy_lower], atmos.y[idy_upper])

            # Do linear interpolation to find τ and S

            # exinction at centre and upwind point
            # z is 1st coordinate, changes rows in vals. y is 2nd coordinate
            α_vals = [α_lower[idx_upwind, idy_lower]   α_lower[idx_upwind, idy_upper]
                      α_upper[idx_upwind, idy_lower]   α_upper[idx_upwind, idy_upper]]

            α_centre = α_lower[idx, idy]
            α_upwind = bilinear(z_upwind, y_upwind, z_bounds, y_bounds, α_vals)

            # Find the Δτ optical path from upwind to grid point
            Δτ_upwind = trapezoidal(r, α_centre, α_upwind)

            w1, w2 =  weights(Δτ_upwind)
        end

        # Update ghost zones
        # I[idx, 1] = I[idx, end-1]
        # I[idx, end] = I[idx, 2]
    end

    # Update ghost zones
    # I[1,:] = I[end-1,:]
    # I[end,:] = I[2,:]
end

"""
    function xz_up_ray(k::Vector{Float64}, idz::Int, sign_x::Int, sign_y::Int,
                       α::Array{<:PerLength, 3}, atmos::Atmosphere)

Ray moving upwards, upwind point intersecting with xz plane. Assumes angle
in radians.
"""
function xz_up_ray(k::Vector{Float64}, idz::Int, sign_x::Int, sign_y::Int,
                   α::Array{<:PerLength, 3}, atmos::Atmosphere)

    nx = length(atmos.x)
    ny = length(atmos.y)

    Δy = atmos.y[2] - atmos.y[1]

    # Loop direction
    start_x, stop_x = range_bounds(sign_x, nx)
    start_y, stop_y = range_bounds(sign_y, ny)

    # Z center and interpolation position
    z_centre = atmos.z[idz]
    idz_lower = idz-1

    # calculate the length to upwind position
    r = abs(Δy/k[3])
    z_increment = r*k[1]
    x_increment = r*k[2]

    z_upwind = z_centre + z_increment

    z_bounds = (atmos.z[idz_lower], z_centre)

    α_lower = α[idz_lower,:,:]
    α_upper = α[idz,:,:]

    for idy in start_y:sign_y:stop_y
        y_centre = atmos.y[idy]
        # x upwind position
        idy_upwind = idy+sign_y
        y_upwind = atmos.y[idy_upwind]
        for idx in start_x:sign_x:stop_x
            # Centre point (c)
            x_centre = atmos.x[idx]

            # Lower corner to interpolate from
            idx_lower = idx-Int((sign_x+1)/2)
            idx_upper = idx_lower+1

            x_upwind = x_centre + x_increment

            x_bounds = (atmos.x[idx_lower], atmos.x[idx_upper])

            # Do linear interpolation to find τ and S

            # exinction at centre and upwind point
            # z is 1st coordinate, changes rows in vals. y is 2nd coordinate
            α_vals = [α_lower[idx_lower, idy_upwind]   α_lower[idx_upper, idy_upwind]
                      α_upper[idx_lower, idy_upwind]   α_upper[idx_upper, idy_upwind]]

            α_centre = α_upper[idx, idy]
            α_upwind = bilinear(z_upwind, x_upwind, z_bounds, x_bounds, α_vals)

            # Find the Δτ optical path from upwind to grid point
            Δτ_upwind = trapezoidal(r, α_centre, α_upwind)

            w1, w2 =  weights(Δτ_upwind)
        end
        # Update ghost zones
        # I[1, idy] = I[end-1, idy]
        # I[end, idy] = I[2, idy]
    end

    # Update ghost zones
    # I[:, 1] = I[:, end-1]
    # I[:, end] = I[:, 2]
end

"""
    function xz_down_ray(θ::AbstractFloat, ϕ::AbstractFloat, idz::Int, sign_x::Int,
                           sign_y::Int, I_0::Matrix{<:UnitsIntensity_λ}, S_0::Array{<:UnitsIntensity_λ, 3},
                           α::Array{<:PerLength, 3}, atmos::Atmosphere)

Ray moving downwards, upwind point intersecting with xz plane. Assumes angle
in radians.
"""
function xz_down_ray(k::Vector{Float64}, idz::Int, sign_x::Int, sign_y::Int,
                     α::Array{<:PerLength, 3}, atmos::Atmosphere)

    nx = length(atmos.x)
    ny = length(atmos.y)

    Δy = atmos.y[2] - atmos.y[1]

    # Loop direction
    start_x, stop_x = range_bounds(sign_x, nx)
    start_y, stop_y = range_bounds(sign_y, ny)


    # Z center and interpolation position
    z_centre = atmos.z[idz]
    idz_upper = idz+1

    # Length to plane
    # Δx = atmos.x[2] - atmos.x[1]

    # calculate the length to upwind position
    r = abs(Δy/k[3])
    z_increment = r*k[1]
    x_increment = r*k[2]

    z_upwind = z_centre + z_increment

    z_bounds = (z_centre, atmos.z[idz_upper])

    α_lower = α[idz,:,:]
    α_upper = α[idz_upper,:,:]

    for idy in start_y:sign_y:stop_y
        y_centre = atmos.y[idy]
        # x upwind position
        idy_upwind = idy+sign_y
        y_upwind = atmos.y[idy_upwind]
        for idx in start_x:sign_x:stop_x
            # Centre point (c)
            x_centre = atmos.x[idx]

            # Lower corner to interpolate from
            idx_lower = idx-Int((sign_x+1)/2)
            idx_upper = idx_lower+1

            x_upwind = x_centre + x_increment

            x_bounds = (atmos.x[idx_lower], atmos.x[idx_upper])

            # Do linear interpolation to find τ and S

            # exinction at centre and upwind point
            # z is 1st coordinate, changes rows in vals. y is 2nd coordinate
            α_vals = [α_lower[idx_lower, idy_upwind]   α_lower[idx_upper, idy_upwind]
                      α_upper[idx_lower, idy_upwind]   α_upper[idx_upper, idy_upwind]]

            α_centre = α_upper[idx, idy]
            α_upwind = bilinear(z_upwind, x_upwind, z_bounds, x_bounds, α_vals)

            # Find the Δτ optical path from upwind to grid point
            Δτ_upwind = trapezoidal(r, α_centre, α_upwind)

            w1, w2 =  weights(Δτ_upwind)
        end

        # Update ghost zones
        # I[1, idy] = I[end-1, idy]
        # I[end, idy] = I[2, idy]
    end

    # Update ghost zones
    # I[:, 1] = I[:, end-1]
    # I[:, end] = I[:, 2]
end
