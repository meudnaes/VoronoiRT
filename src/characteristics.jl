include("functions.jl")
include("rays.jl")

function short_characteristic_ray(θ, ϕ, I_0, S_0, α, atmos)
    # First find out which wall the ray hits, to determine which direction to
    # interpolate. Exploit the fact that xy grid is equidistand partitioned for
    # all z, however z is in general not equipartitioned. Find out which Δ - r
    # is largest or smallest???

    # Convert to radians if input angle is in degrees...
    if degrees == true
        θ = θ*pi/180
        ϕ = ϕ*pi/180
    end

    # Allocate array for new intensity
    I = zero(I_0)

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

    # I know that Δx = Δy = constant for all grid points
    Δxy = atmos.x[2] - atmos.x[1]
    r_x = Δxy/(cos(ϕ)*sin(θ))
    r_y = Δxy/(sin(ϕ)*sin(θ))

    # find out which plane the upwind part of the ray intersects
    sign_z = z_intersect(θ)
    sign_x, sign_y = xy_intersect(ϕ)

    # Cut point
    idz_u = idz + 2*sign_z + 1

    # loop through atmosphere
    for k in 2:length(atmos.z) - 1
        # calculate length until ray hits z plane

        # Length to plane
        Δz = sign_z*(atmos.z[idz_u] - atmos.z[k])
        r_z = Δz/cos(ϕ)

        # This finds which plane the ray intersects with
        plane_cut = argmin([r_z, r_x, r_y])
        for i in 1:length(atmos.x)
            for j in 1:length(atmos.y)
                if plane_cut == 1
                    ΔI = z_ray(θ, ϕ, k, i, j, sign_z, sign_x, sign_y,
                                                                 I_0, S_0, α, atmos)

                elseif plane_cut==2
                    ΔI = x_ray(θ, ϕ, k, i, j, sign_z, sign_x, sign_y,
                                                                 I_0, S_0, α, atmos)

                elseif plane_cut==3
                    ΔI = y_ray(θ, ϕ, k, i, j, sign_z, sign_x, sign_y,
                                                                 I_0, S_0, α, atmos)
                end

                I[k, i, j] = ΔI
            end
        end
        print("$k \r")
    end
    return I
end
