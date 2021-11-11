include("voronoi_utils.jl")
include("functions.jl")

using Distances
using LinearAlgebra

function irregular_SC_up(sites::VoronoiSites, cells::AbstractArray{VoronoiCell}, layers::Vector,
                         I_0::AbstractArray, S_0::AbstractArray, α_cont::AbstractArray)

    # Inverse distance power law parameter
    p = 3.0

    # Traces rays through an irregular grid
    θ = 10
    ϕ = 170

    # start at the bottom
    # shoot rays through every grid cell

    # precalculate trigonometric functions
    cosθ = cos(θ)
    sinθ = sin(θ)

    cosϕ = cos(ϕ)
    sinϕ = sin(ϕ)

    # Unit vector towards upwind direction of the ray
    k = -[cosθ, cosϕ*sinθ, sinϕ*sinθ]

    # Allocate space for intensity
    I = zero(S_0)

    x_min = sites.x_min
    x_max = sites.x_max
    y_min = sites.y_min
    y_max = sites.y_max

    nx = length(I_0[:,1])
    ny = length(I_0[1,:])

    Δx = (x_max - x_min)/(nx - 1)
    Δy = (y_max - y_min)/(ny - 1)

    max_layer = maximum(layers)

    for layer in 1:max_layer
        lower_idx = searchsortedfirst(layers, layer)
        upper_idx = searchsortedfirst(layers, layer+1)
        if layer == 1
            for i in lower_idx:upper_idx-1
                cell = cells[i]
                # coordinate
                position = cell.position

                # 1st layer
                z_upwind = sites.z_min

                Δz = position[1] - z_upwind
                r = abs(Δz/cosθ)

                x_upwind = position[2] + r*cosϕ*sinθ
                if x_upwind > sites.x_max
                    x_upwind = sites.x_min + (x_upwind - sites.x_max)
                elseif x_upwind < sites.x_min
                    x_upwind = sites.x_max + (x_upwind - sites.x_min)
                end

                y_upwind = position[3] + r*sinϕ*sinθ
                if y_upwind > sites.y_max
                    y_upwind = sites.y_min + (y_upwind - sites.y_max)
                elseif y_upwind < sites.y_min
                    y_upwind = sites.y_max + (y_upwind - sites.y_min)
                end

                upwind_position = [z_upwind, x_upwind, y_upwind]

                distances = Vector{Float64}(undef, cell.n)
                # This has length = cell.n - 1
                neighbours = cell.neighbours[cell.neighbours .> 0]

                for (j, neighbour) in enumerate(neighbours)
                    distances[j] = euclidean(upwind_position, sites.positions[:,neighbour])
                end

                distances[cell.n] = euclidean(upwind_position, position)

                # Pass α as an array
                α_centre = α_cont[cell.ID]
                α_upwind = inv_dist_itp_test([neighbours..., cell.ID],
                                             distances, p, α_cont)

                # Find the Δτ optical path from upwind to grid point
                Δτ_upwind = trapezoidal(r, α_centre, α_upwind)

                S_centre = S_0[cell.ID]
                S_upwind = inv_dist_itp_test([neighbours..., cell.ID],
                                             distances, p, S_0)

                w1, w2 =  weights(Δτ_upwind)
                a_ijk, b_ijk, c_ijk = coefficients(w1, w2, Δτ_upwind)

                # Interpolate I_0
                # find x and y on the I_0 values...
                idx_0 = floor(Int, abs(x_upwind/Δx)) + 1
                idy_0 = floor(Int, abs(y_upwind/Δy)) + 1

                x_bounds = Δx .* [idx_0-1, idx_0]
                y_bounds = Δy .* [idy_0-1, idy_0]

                I_vals = [I_0[idx_0, idy_0]     I_0[idx_0, idy_0+1]
                          I_0[idx_0+1, idy_0]   I_0[idx_0+1, idy_0+1]]

                I_upwind = bilinear(x_upwind, y_upwind, x_bounds, y_bounds, I_vals)

                I[cell.ID] = a_ijk*S_upwind + b_ijk*S_centre + c_ijk*I_upwind
            end

        elseif 1 < layer <= max_layer
            for i in lower_idx:upper_idx-1
                cell = cells[i]

                # coordinate
                position = cell.position

                # Find the intersection and the neighbour the ray is coming from
                upwind_index, upwind_position = _rayIntersection(k, cell, sites)

                distances = Vector{Float64}(undef, cell.n+1)

                # Control for periodic boundaries, fix !!
                for (j, neighbour) in enumerate(cell.neighbours)
                    if neighbour > 0
                        distances[j] = euclidean(upwind_position, sites.positions[:,neighbour])
                    elseif neighbour == -6
                        # Boundary
                        z_upwind = sites.z_max
                        z = position[1]
                        distances[j] = abs(k[1]/(z_upwind - z))
                    end
                end

                distances[cell.n+1] = euclidean(upwind_position, position)

                # Pass α as an array
                α_centre = α_cont[cell.ID]
                α_upwind = inv_dist_itp_test([cell.neighbours..., cell.ID],
                                             distances, p, α_cont)

                # Find the Δτ optical path from upwind to grid point
                Δτ_upwind = trapezoidal(distances[cell.n], α_centre, α_upwind)

                S_centre = S_0[cell.ID]
                S_upwind = inv_dist_itp_test([cell.neighbours..., cell.ID],
                                             distances, p, S_0)

                w1, w2 =  weights(Δτ_upwind)
                a_ijk, b_ijk, c_ijk = coefficients(w1, w2, Δτ_upwind)

                indices = greater_than(I[cell.neighbours[cell.neighbours .> 0]], 0)
                I_vals = I[cell.neighbours[cell.neighbours .> 0][indices]]
                I_distances = distances[indices]

                I_upwind = inv_dist_itp_test(Vector(1:length(I_vals)), I_distances, p, I_vals)

                I[cell.ID] = a_ijk*S_upwind + b_ijk*S_centre + c_ijk*I_upwind
            end
        end
    end

    return I
end


function _rayIntersection(k::AbstractArray, cell::VoronoiCell, sites::VoronoiSites)
    # k is unit vector towards upwind direction of the ray

    p_r = cell.position

    x_r_r = sites.x_max - p_r[2]
    x_r_l = p_r[2] - sites.x_min

    y_r_r = sites.y_max - p_r[3]
    y_r_l = p_r[3] - sites.y_min

    s = Vector{Float64}(undef, cell.n)
    for (i, neighbour) in enumerate(cell.neighbours)
        if neighbour > 0
            # Natural neighbor position
            p_i = sites.positions[:,neighbour]

            x_i_r = sites.x_max - p_i[2]
            x_i_l = p_i[2] - sites.x_min

            # Test for periodic
            if x_r_r + x_i_l < p_r[2] - p_i[2]
                p_i[2] = sites.x_max + p_i[2] - sites.x_min
            elseif x_r_l + x_i_r < p_i[2] - p_r[2]
                p_i[2] = sites.x_min + sites.x_max - p_i[2]
            end

            y_i_r = sites.y_max - p_i[3]
            y_i_l = p_i[3] - sites.y_min

            # Test for periodic
            if y_r_r + y_i_l < p_r[3] - p_i[3]
                p_i[3] = sites.y_max + p_i[3] - sites.y_min
            elseif y_r_l + y_i_r < p_i[3] - p_r[3]
                p_i[3] = sites.y_min + sites.y_max - p_i[3]
            end

            # Calculate plane biseting site and neighbor
            # Normal vector
            n = p_i - p_r

            # Point on the plane
            p = (p_i + p_r)/2

            s[i] = dot(n, p - p_r)/dot(n, k)
        elseif neighbour == -6
            # Boundary
            z_upwind = sites.z_max
            z = p_r[1]
            s[i] = abs(k[1]/(z_upwind - z))
        elseif neighbour == -5
            # Boundary
            z_upwind = sites.z_min
            z = p_r[1]
            s[i] = abs(k[1]/(z_upwind - z))
        end
    end

    # Find the smallst non-negative s
    index, s_q = smallestNonNegative(s)

    q = p_r + s_q*k

    return index, q
end

function greater_than(arr::AbstractArray, treshold)
    indices = Vector{Int}(undef, length(arr))
    j = 0
    for i in 1:length(arr)
        if arr[i] > treshold
            j += 1
            indices[j] = i
        end
    end
    return indices[1:j]
end
