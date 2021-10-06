include("functions.jl")
include("rays.jl")


"""
TODO: fix boundary at the top. Add "ghost layer" for rays moving from xz or yz
planes... how to? Start with last I_0, do extra iteration:

for i in 1:nx
    calculate SC
    update I_0
i = 1
    calculate SC from last I_0
"""

function short_characteristic_ray(θ, ϕ, S_0, α, atmos; degrees)
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
    I = zero(S_0)

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

    if ϕ > π/2
        # Boundary condition
        I_0 = S_0[1,:,:]

        # loop upwards through atmosphere
        for idz in 2:length(atmos.z)-1
            # Grid diffference in z
            Δz = atmos.z[idz] - atmos.z[idz-1]

            # calculate length until ray hits z plane
            r_z = Δz/cos(ϕ)

            # This finds which plane the ray intersects with
            plane_cut = argmin([r_z, r_x, r_y])
            if plane_cut == 1
                ΔI = z_up_ray(θ, ϕ, idz, sign_x, sign_y, I_0, S_0, α, atmos)

            elseif plane_cut==2
                ΔI = x_up_ray(θ, ϕ, idz, sign_x, sign_y, I_0, S_0, α, atmos)

            elseif plane_cut==3
                ΔI = y_up_ray(θ, ϕ, idz, sign_x, sign_y, I_0, S_0, α, atmos)
            end

            I_0 = I[idz]
            percent = trunc(Int, 100*idz/(length(atmos.z)-1))
            print("\t\t$percent% \r")
        end
    elseif ϕ < π/2
        # Boundary condition
        I_0 = zero(S_0[end,:,:])
        # loop downwards through atmosphere
        for idz in length(atmos.z)-1:-1:2
            # Grid diffference in z
            Δz = atmos.z[idz+1] - atmos.z[idz]

            # calculate length until ray hits z plane
            r_z = Δz/cos(ϕ)

            # This finds which plane the ray intersects with
            plane_cut = argmin([r_z, r_x, r_y])
            if plane_cut == 1
                I[idz,:,:] = z_down_ray(θ, ϕ, idz, sign_x, sign_y, I_0, S_0, α, atmos)

            elseif plane_cut==2
                I[idz,:,:] = x_down_ray(θ, ϕ, idz, sign_x, sign_y, I_0, S_0, α, atmos)

            elseif plane_cut==3
                I[idz,:,:] = y_down_ray(θ, ϕ, idz, sign_x, sign_y, I_0, S_0, α, atmos)
            end

            I_0 = I[idz,:,:]
            percent = trunc(Int, 100*(length(atmos.z)-idz+1)/(length(atmos.z)-1))
            print("\t\t$percent% \r")
        end
    end

    return I
end
