include("functions.jl")
include("rays.jl")


function short_characteristics_up(θ, ϕ, S_0, α, atmos; degrees)
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
    # Δxy = atmos.x[2] - atmos.x[1]
    r_x = Δx/(cos(ϕ)*sin(θ))
    r_y = Δy/(sin(ϕ)*sin(θ))

    # find out direction in xy the ray moves (1 for positive, -1 for negative)
    sign_x, sign_y = xy_intersect(ϕ)

    # Boundary condition
    I_0 = S_0[1,:,:]

    # loop upwards through atmosphere
    for idz in 2:length(atmos.z)-1
        # Grid diffference in z
        Δz = atmos.z[idz] - atmos.z[idz-1]

        # calculate length until ray hits z plane
        r_z = abs(Δz/cos(θ))

        # This finds which plane the ray intersects with
        plane_cut = argmin([r_z, r_x, r_y])
        if plane_cut == 1
            if idz == 2
                println("xy_up")
            end
            I[idz,:,:] = xy_up_ray(θ, ϕ, idz, sign_x, sign_y, I_0, S_0, α, atmos)
        elseif plane_cut==2
            if idz == 2
                println("yz_up")
            end
            I[idz,:,:] = yz_up_ray(θ, ϕ, idz, sign_x, sign_y, I_0, S_0, α, atmos)
        elseif plane_cut==3
            if idz == 2
                println("xz_up")
            end
            I[idz,:,:] = xz_up_ray(θ, ϕ, idz, sign_x, sign_y, I_0, S_0, α, atmos)
        end

        I_0 = I[idz,:,:]
        percent = trunc(Int, 100*idz/(length(atmos.z)-1))
        # print("\t\t$percent% \r")
    end

    return I
end


function short_characteristics_down(θ, ϕ, S_0, α, atmos; degrees)
    # Convert to radians if input angle is in degrees...
    if degrees == true
        θ = θ*pi/180
        ϕ = ϕ*pi/180
    end

    # Allocate array for new intensity
    I = zero(S_0)

    # I know that Δx = Δy = constant for all grid points
    # Δxy = atmos.x[2] - atmos.x[1]
    r_x = Δx/(cos(ϕ)*sin(θ))
    r_y = Δy/(sin(ϕ)*sin(θ))

    # find out which plane the upwind part of the ray intersects
    sign_x, sign_y = xy_intersect(ϕ)

    # Boundary condition
    I_0 = zero(S_0[end,:,:])
    # loop downwards through atmosphere
    for idz in length(atmos.z)-1:-1:2
        # Grid diffference in z
        Δz = atmos.z[idz+1] - atmos.z[idz]

        # calculate length until ray hits z plane
        r_z = abs(Δz/cos(θ))

        # This finds which plane the ray intersects with
        plane_cut = argmin([r_z, r_x, r_y])
        if plane_cut == 1
            if idz == length(atmos.z)-1
                println("xy_down")
            end
            I[idz,:,:] = xy_down_ray(θ, ϕ, idz, sign_x, sign_y, I_0, S_0, α, atmos)
        elseif plane_cut==2
            if idz == length(atmos.z)-1
                println("yz_down")
            end
            I[idz,:,:] = yz_down_ray(θ, ϕ, idz, sign_x, sign_y, I_0, S_0, α, atmos)
        elseif plane_cut==3
            if idz == length(atmos.z)-1
                println("xz_down")
            end
            I[idz,:,:] = xz_down_ray(θ, ϕ, idz, sign_x, sign_y, I_0, S_0, α, atmos)
        end

        I_0 = I[idz,:,:]
        percent = trunc(Int, 100*(length(atmos.z)-idz+1)/(length(atmos.z)-1))
        # print("\t\t$percent% \r")
    end

    return I
end
