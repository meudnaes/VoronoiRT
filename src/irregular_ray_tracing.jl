include("voronoi_utils.jl")
include("functions.jl")

using Distances
using LinearAlgebra

function SC_NNintersection_up(sites::VoronoiSites,
                         I_0::AbstractMatrix, S_0::AbstractMatrix, α_0::AbstractMatrix,
                         S::AbstractVector, α_cont::AbstractVector)

    # Inverse distance power law parameter
    p = 3.0

    # Weights for interpolation
    ω = [1, 1]

    # Sweeeeeps
    n_sweeps = 2

    # Traces rays through an irregular grid
    θ = 10*π/180
    ϕ = 10*π/180

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
    I = zero(S)

    x_min = sites.x_min
    x_max = sites.x_max
    y_min = sites.y_min
    y_max = sites.y_max

    nx = length(I_0[:,1])
    ny = length(I_0[1,:])

    Δx = (x_max - x_min)/(nx - 1)
    Δy = (y_max - y_min)/(ny - 1)

    max_layer = maximum(sites.layers)

    for layer in 1:max_layer
        lower_idx = searchsortedfirst(sites.layers, layer)
        upper_idx = searchsortedfirst(sites.layers, layer+1)
        if layer == 1
            z_upwind = sites.z_min
            for i in lower_idx:upper_idx-1
                # coordinate
                position = sites.positions[:,i]

                # 1st layer
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

                # Find distance between border and cell centre
                r = euclidean(upwind_position, position)

                # Interpolate boundary values
                # find x and y on the I_0 values...
                idx_0 = floor(Int, abs(x_upwind/Δx)) + 1
                idy_0 = floor(Int, abs(y_upwind/Δy)) + 1

                x_bounds = Δx .* [idx_0-1, idx_0]
                y_bounds = Δy .* [idy_0-1, idy_0]

                # Pass α as an array
                α_centre = α_cont[i]
                α_vals = [α_0[idx_0, idy_0]     α_0[idx_0, idy_0+1]
                          α_0[idx_0+1, idy_0]   α_0[idx_0+1, idy_0+1]]
                α_upwind = bilinear(x_upwind, y_upwind, x_bounds, y_bounds, α_vals)

                # Find the Δτ optical path from upwind to grid point
                Δτ_upwind = trapezoidal(r, α_centre, α_upwind)

                S_centre = S[i]
                S_vals = [S_0[idx_0, idy_0]     S_0[idx_0, idy_0+1]
                          S_0[idx_0+1, idy_0]   S_0[idx_0+1, idy_0+1]]
                S_upwind = bilinear(x_upwind, y_upwind, x_bounds, y_bounds, S_vals)

                w1, w2 =  weights(Δτ_upwind)
                a_ijk, b_ijk, c_ijk = coefficients(w1, w2, Δτ_upwind)


                I_vals = [I_0[idx_0, idy_0]     I_0[idx_0, idy_0+1]
                          I_0[idx_0+1, idy_0]   I_0[idx_0+1, idy_0+1]]

                I_upwind = bilinear(x_upwind, y_upwind, x_bounds, y_bounds, I_vals)

                I[i] = a_ijk*S_upwind + b_ijk*S_centre + c_ijk*I_upwind
            end

        else
            for sweep in 1:n_sweeps
                for i in lower_idx:upper_idx-1
                    # coordinate
                    position = sites.positions[:,i]

                    # number of neighbours
                    n_neighbours = sites.neighbours[i,1]
                    neighbours = sites.neighbours[i, 2:n_neighbours+1]
                    cell_neighbours = neighbours[neighbours .> 0]
                    neighbour_height = sites.positions[1, cell_neighbours]

                    # Find the intersection and the neighbour the ray is coming from
                    upwind_position, upwind_index = ray_intersection(k, neighbours, position, sites, i)


                    distances = Vector{Float64}(undef, n_neighbours+1)
                    distances[1:n_neighbours] = upwind_distances(neighbours, n_neighbours, i,
                                                                 upwind_position, sites)

                    distances[n_neighbours+1] = euclidean(upwind_position, position)

                    # Pass α as an array
                    α_centre = α_cont[i]
                    α_upwind = inv_dist_itp([neighbours..., i],
                                            distances, p, α_cont)

                    # Find the Δτ optical path from upwind to grid point
                    Δτ_upwind = trapezoidal(distances[n_neighbours+1], α_centre, α_upwind)

                    S_centre = S[i]
                    S_upwind = inv_dist_itp([neighbours..., i],
                                            distances, p, S)

                    w1, w2 =  weights(Δτ_upwind)
                    a_ijk, b_ijk, c_ijk = coefficients(w1, w2, Δτ_upwind)

                    indices = less_than(neighbour_height, upwind_position[1])
                    I_vals = I[cell_neighbours[indices]]
                    I_distances = distances[indices]
                    I_neighbours = inv_dist_itp(collect(1:length(indices)), I_distances, p, I_vals)
                    I_voronoi = I[upwind_index]
                    I_upwind = 0.9*I_voronoi + 0.1*I_neighbours
                    if length(indices) < 1
                        I_upwind = 0
                    end
                    I[i] = a_ijk*S_upwind + b_ijk*S_centre + c_ijk*I_upwind
                end
            end
        end
    end
    return I
end

function SC_Delaunay_up(sites::VoronoiSites,
                    I_0::AbstractMatrix, S_0::AbstractMatrix, α_0::AbstractMatrix,
                    S::AbstractVector, α_cont::AbstractVector)

    #
    # Inverse distance power law parameter
    p = 3.0

    #
    Nran = 3

    # Sweeeeeps
    n_sweeps = 3

    # Traces rays through an irregular grid
    θ = 10*π/180
    ϕ = 10*π/180

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
    I = zero(S)

    x_min = sites.x_min
    x_max = sites.x_max
    y_min = sites.y_min
    y_max = sites.y_max

    nx = length(I_0[:,1])
    ny = length(I_0[1,:])

    Δx = (x_max - x_min)/(nx - 1)
    Δy = (y_max - y_min)/(ny - 1)

    p = sortperm(sites.layers_up)
    layers_sorted = sites.layers_up[p]
    max_layer = layers_sorted[end]

    for layer in 1:max_layer
        lower_idx = searchsortedfirst(layers_sorted, layer)
        upper_idx = searchsortedfirst(layers_sorted, layer+1)
        if layer == 1
            z_upwind = sites.z_min
            for i in lower_idx:upper_idx-1
                # coordinate
                idx = p[i]
                position = sites.positions[:,idx]

                # 1st layer
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

                # Find distance between border and cell centre
                r = euclidean(upwind_position, position)

                # Interpolate boundary values
                # find x and y on the I_0 values...
                idx_0 = floor(Int, abs(x_upwind/Δx)) + 1
                idy_0 = floor(Int, abs(y_upwind/Δy)) + 1

                x_bounds = Δx .* [idx_0-1, idx_0]
                y_bounds = Δy .* [idy_0-1, idy_0]

                # Pass α as an array
                α_centre = α_cont[idx]
                α_vals = [α_0[idx_0, idy_0]     α_0[idx_0, idy_0+1]
                          α_0[idx_0+1, idy_0]   α_0[idx_0+1, idy_0+1]]
                α_upwind = bilinear(x_upwind, y_upwind, x_bounds, y_bounds, α_vals)

                # Find the Δτ optical path from upwind to grid point
                Δτ_upwind = trapezoidal(r, α_centre, α_upwind)

                S_centre = S[idx]
                S_vals = [S_0[idx_0, idy_0]     S_0[idx_0, idy_0+1]
                          S_0[idx_0+1, idy_0]   S_0[idx_0+1, idy_0+1]]
                S_upwind = bilinear(x_upwind, y_upwind, x_bounds, y_bounds, S_vals)

                w1, w2 =  weights(Δτ_upwind)
                a_ijk, b_ijk, c_ijk = coefficients(w1, w2, Δτ_upwind)


                I_vals = [I_0[idx_0, idy_0]     I_0[idx_0, idy_0+1]
                          I_0[idx_0+1, idy_0]   I_0[idx_0+1, idy_0+1]]

                I_upwind = bilinear(x_upwind, y_upwind, x_bounds, y_bounds, I_vals)
                I[i] = a_ijk*S_upwind + b_ijk*S_centre + c_ijk*I_upwind
            end

        else
            for sweep in 1:n_sweeps
                for i in lower_idx:upper_idx-1
                    # coordinate
                    idx = p[i]
                    position = sites.positions[:,idx]

                    # number of neighbours
                    n_neighbours = sites.neighbours[idx, 1]
                    neighbours = sites.neighbours[idx, 2:n_neighbours+1]

                    smallest∠, ∠_indices = smallest_angle(position, neighbours, -k, sites, 2)
                    # ∠_index = smallest_angle(position, neighbours, -k, sites, 2)
                    I[idx] = 0
                    for rn in 1:Nran
                        ∠_index = choose_random(smallest∠, ∠_indices)
                        upwind_index = neighbours[∠_index]
                        upwind_position = sites.positions[:, upwind_index]


                        # Find the intersection and the neighbour the ray is coming from
                        # upwind_position, upwind_index = ray_intersection(k, neighbours, position, sites, i)

                        # Pass α as an array
                        α_centre = α_cont[idx]
                        α_upwind = α_cont[upwind_index]

                        r = euclidean(position, upwind_position)
                        # Find the Δτ optical path from upwind to grid point
                        Δτ_upwind = trapezoidal(r, α_centre, α_upwind)

                        S_centre = S[idx]
                        S_upwind = S[upwind_index]

                        w1, w2 =  weights(Δτ_upwind)
                        a_ijk, b_ijk, c_ijk = coefficients(w1, w2, Δτ_upwind)

                        I_upwind = I[upwind_index]
                        I[idx] += (a_ijk*S_upwind + b_ijk*S_centre + c_ijk*I_upwind)/Nran
                        # maybe try a mix of this and the ray intersection?
                    end
                end
            end
        end
    end
    return I
end

function SC_Delaunay_down(sites::VoronoiSites,
                    I_0::AbstractMatrix, S_0::AbstractMatrix, α_0::AbstractMatrix,
                    S::AbstractVector, α_cont::AbstractVector)

    #
    # Inverse distance power law parameter
    p = 3.0

    #
    Nran = 3

    # Sweeeeeps
    n_sweeps = 3

    # Traces rays through an irregular grid
    θ = 170*π/180
    ϕ = 10*π/180

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
    I = zero(S)

    x_min = sites.x_min
    x_max = sites.x_max
    y_min = sites.y_min
    y_max = sites.y_max

    nx = length(I_0[:,1])
    ny = length(I_0[1,:])

    Δx = (x_max - x_min)/(nx - 1)
    Δy = (y_max - y_min)/(ny - 1)

    p = sortperm(sites.layers_down, rev=true)
    layers_sorted = sites.layers_down[p]
    max_layer = layers_sorted[end]

    for layer in max_layer:-1:1
        lower_idx = searchsortedfirst(layers_sorted, layer)
        upper_idx = searchsortedfirst(layers_sorted, layer+1)
        if layer == max_layer
            z_upwind = sites.z_max
            for i in upper_idx-1:-1:lower_idx
                # coordinate
                idx=p[i]
                position = sites.positions[:,idx]

                # 1st layer
                Δz = z_upwind - position[1]
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

                # Find distance between border and cell centre
                r = euclidean(upwind_position, position)

                # Interpolate boundary values
                # find x and y on the I_0 values...
                idx_0 = floor(Int, abs(x_upwind/Δx)) + 1
                idy_0 = floor(Int, abs(y_upwind/Δy)) + 1

                x_bounds = Δx .* [idx_0-1, idx_0]
                y_bounds = Δy .* [idy_0-1, idy_0]

                # Pass α as an array
                α_centre = α_cont[idx]
                α_vals = [α_0[idx_0, idy_0]     α_0[idx_0, idy_0+1]
                          α_0[idx_0+1, idy_0]   α_0[idx_0+1, idy_0+1]]
                α_upwind = bilinear(x_upwind, y_upwind, x_bounds, y_bounds, α_vals)

                # Find the Δτ optical path from upwind to grid point
                Δτ_upwind = trapezoidal(r, α_centre, α_upwind)

                S_centre = S[idx]
                S_vals = [S_0[idx_0, idy_0]     S_0[idx_0, idy_0+1]
                          S_0[idx_0+1, idy_0]   S_0[idx_0+1, idy_0+1]]
                S_upwind = bilinear(x_upwind, y_upwind, x_bounds, y_bounds, S_vals)

                w1, w2 =  weights(Δτ_upwind)
                a_ijk, b_ijk, c_ijk = coefficients(w1, w2, Δτ_upwind)


                I_vals = [I_0[idx_0, idy_0]     I_0[idx_0, idy_0+1]
                          I_0[idx_0+1, idy_0]   I_0[idx_0+1, idy_0+1]]

                I_upwind = bilinear(x_upwind, y_upwind, x_bounds, y_bounds, I_vals)
                I[i] = a_ijk*S_upwind + b_ijk*S_centre + c_ijk*I_upwind
            end

        else
            for sweep in 1:n_sweeps
                for i in upper_idx-1:-1:lower_idx
                    # coordinate
                    position = sites.positions[:,idx]

                    # number of neighbours
                    n_neighbours = sites.neighbours[idx, 1]
                    neighbours = sites.neighbours[idx, 2:n_neighbours+1]

                    smallest∠, ∠_indices = smallest_angle(position, neighbours, -k, sites, 2)
                    # ∠_index = smallest_angle(position, neighbours, -k, sites, 2)
                    I[i] = 0
                    for rn in 1:Nran
                        ∠_index = choose_random(smallest∠, ∠_indices)
                        upwind_index = neighbours[∠_index]
                        upwind_position = sites.positions[:, upwind_index]


                        # Find the intersection and the neighbour the ray is coming from
                        # upwind_position, upwind_index = ray_intersection(k, neighbours, position, sites, i)

                        # Pass α as an array
                        α_centre = α_cont[idx]
                        α_upwind = α_cont[upwind_index]

                        r = euclidean(position, upwind_position)
                        # Find the Δτ optical path from upwind to grid point
                        Δτ_upwind = trapezoidal(r, α_centre, α_upwind)

                        S_centre = S[idx]
                        S_upwind = S[upwind_index]

                        w1, w2 =  weights(Δτ_upwind)
                        a_ijk, b_ijk, c_ijk = coefficients(w1, w2, Δτ_upwind)

                        I_upwind = I[upwind_index]
                        I[idx] += (a_ijk*S_upwind + b_ijk*S_centre + c_ijk*I_upwind)/Nran
                        # maybe try a mix of this and the ray intersection?
                        # Need to think more
                    end
                end
            end
        end
    end
    return I
end
