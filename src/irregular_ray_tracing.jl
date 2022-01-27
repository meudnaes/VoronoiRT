include("voronoi_utils.jl")
include("functions.jl")

using Distances
using LinearAlgebra

function SC_NNintersection_up(sites::VoronoiSites,
                         I_0::AbstractMatrix, S_0::AbstractMatrix, α_0::AbstractMatrix,
                         S::AbstractVector, α::AbstractVector)

    # Inverse distance power law parameter
    p = 3.0

    # Weights for interpolation
    ω = [1, 1]

    # sweeps
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

    perm = sortperm(sites.layers_up)
    layers_sorted = sites.layers_up[perm]
    max_layer = layers_sorted[end]

    lower_boundary(I, sites, I_0, S_0, α_0, S, α, k)

    for layer in 2:max_layer
        lower_idx = searchsortedfirst(layers_sorted, layer)
        upper_idx = searchsortedfirst(layers_sorted, layer+1)
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
                α_centre = α[i]
                α_upwind = inv_dist_itp([neighbours..., i],
                                        distances, p, α)

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
    return I
end

function SC_Delaunay_up(sites::VoronoiSites,
                    I_0::AbstractMatrix, S_0::AbstractMatrix, α_0::AbstractMatrix,
                    S::AbstractVector, α::AbstractVector, k::AbstractVector, n_sweeps::Int, Nran::Int)

    #
    # Inverse distance power law parameter
    p = 3.0

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

    perm = sortperm(sites.layers_up)
    layers_sorted = sites.layers_up[perm]
    max_layer = layers_sorted[end]

    lower_boundary(I, sites, I_0, S_0, α_0, S, α, k)

    for layer in 2:max_layer
        lower_idx = searchsortedfirst(layers_sorted, layer)
        upper_idx = searchsortedfirst(layers_sorted, layer+1)
        for sweep in 1:n_sweeps
            for i in lower_idx:upper_idx-1
                # coordinate
                idx = perm[i]
                position = sites.positions[:,idx]

                # number of neighbours
                n_neighbours = sites.neighbours[idx, 1]
                neighbours = sites.neighbours[idx, 2:n_neighbours+1]

                smallest∠, ∠_indices, upwind_positions = smallest_angle(position, neighbours, -k, sites, 2)
                # ∠_index = smallest_angle(position, neighbours, -k, sites, 2)
                I[idx] = 0
                for rn in 1:Nran
                    ∠_index = choose_random(smallest∠, ∠_indices)
                    upwind_index = neighbours[∠_indices[∠_index]]
                    upwind_position = upwind_positions[:, ∠_index]


                    # Find the intersection and the neighbour the ray is coming from
                    # upwind_position, upwind_index = ray_intersection(k, neighbours, position, sites, i)

                    # Pass α as an array
                    α_centre = α[idx]
                    α_upwind = α[upwind_index]

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
    return I
end

function SC_Delaunay_down(sites::VoronoiSites,
                        I_0::AbstractMatrix, S_0::AbstractMatrix, α_0::AbstractMatrix,
                        S::AbstractVector, α::AbstractVector, k::AbstractVector, n_sweeps::Int, Nran::Int)

    #
    # Inverse distance power law parameter
    p = 3.0

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

    perm = sortperm(sites.layers_down)
    layers_sorted = sites.layers_down[perm]
    max_layer = layers_sorted[end]

    upper_boundary(I, sites, I_0, S_0, α_0, S, α, k)

    for layer in 2:max_layer
        lower_idx = searchsortedfirst(layers_sorted, layer)
        upper_idx = searchsortedfirst(layers_sorted, layer+1)
        for sweep in 1:n_sweeps
            for i in upper_idx-1:-1:lower_idx
                # coordinate
                idx=perm[i]
                position = sites.positions[:,idx]

                # number of neighbours
                n_neighbours = sites.neighbours[idx, 1]
                neighbours = sites.neighbours[idx, 2:n_neighbours+1]

                smallest∠, ∠_indices, upwind_positions = smallest_angle(position, neighbours, -k, sites, 2)
                # ∠_index = smallest_angle(position, neighbours, -k, sites, 2)
                I[idx] = 0
                for rn in 1:Nran
                    ∠_index = choose_random(smallest∠, ∠_indices)
                    upwind_index = neighbours[∠_indices[∠_index]]
                    upwind_position = upwind_positions[:, ∠_index]


                    # Find the intersection and the neighbour the ray is coming from
                    # upwind_position, upwind_index = ray_intersection(k, neighbours, position, sites, i)

                    # Pass α as an array
                    α_centre = α[idx]
                    α_upwind = α[upwind_index]

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
    return I
end

function Delaunay_up(sites::VoronoiSites,
                        I_0::AbstractVector,
                        S::AbstractVector,
                        α::AbstractVector,
                        k::AbstractVector,
                        n_sweeps::Int)

    # Weighting parameter
    p = 10.0

    # Allocate space for intensity
    I = zero(S)

    x_min = sites.x_min
    x_max = sites.x_max
    y_min = sites.y_min
    y_max = sites.y_max

    perm = sortperm(sites.layers_up)
    layers_sorted = sites.layers_up[perm]
    max_layer = layers_sorted[end]

    lower_idx = searchsortedfirst(layers_sorted, 2)-1

    I[perm[1:lower_idx]] = I_0

    for layer in 2:max_layer
        lower_idx = searchsortedfirst(layers_sorted, layer)
        upper_idx = searchsortedfirst(layers_sorted, layer+1)
        for sweep in 1:n_sweeps
            for i in lower_idx:upper_idx-1
                # coordinate
                idx = perm[i]
                position = sites.positions[:,idx]

                # number of neighbours
                n_neighbours = sites.neighbours[idx, 1]
                neighbours = sites.neighbours[idx, 2:n_neighbours+1]

                dot_products, ∠_indices, upwind_positions = smallest_angle(position, neighbours, -k, sites, 2)
                dot_weights = [dot_products[1]^p, dot_products[2]^p]./sum(dot_products.^p)

                I[idx] = 0u"kW*nm^-1*m^-2"
                for rn in 1:2
                    ∠_index = rn
                    upwind_index = neighbours[∠_indices[∠_index]]
                    upwind_position = upwind_positions[:, ∠_index]


                    # Find the intersection and the neighbour the ray is coming from
                    # upwind_position, upwind_index = ray_intersection(k, neighbours, position, sites, i)

                    # Pass α as an array
                    α_centre = α[idx]
                    α_upwind = α[upwind_index]

                    r = euclidean(position*1.0, upwind_position)
                    # Find the Δτ optical path from upwind to grid point
                    Δτ_upwind = trapezoidal(r, α_centre, α_upwind)

                    S_centre = S[idx]
                    S_upwind = S[upwind_index]

                    w1, w2 =  weights(Δτ_upwind)
                    a_ijk, b_ijk, c_ijk = coefficients(w1, w2, Δτ_upwind)

                    I_upwind = I[upwind_index]
                    I[idx] += (a_ijk*S_upwind + b_ijk*S_centre + c_ijk*I_upwind)*dot_weights[rn]
                    # maybe try a mix of this and the ray intersection?
                end
            end
        end
    end
    return I
end

function Delaunay_down(sites::VoronoiSites,
                        I_0::AbstractVector,
                        S::AbstractVector,
                        α::AbstractVector,
                        k::AbstractVector,
                        n_sweeps::Int)

    # Weighting parameter
    p = 7.0

    # Allocate space for intensity
    I = zero(S)

    x_min = sites.x_min
    x_max = sites.x_max
    y_min = sites.y_min
    y_max = sites.y_max

    perm = sortperm(sites.layers_down)
    layers_sorted = sites.layers_down[perm]
    max_layer = layers_sorted[end]

    lower_idx = searchsortedfirst(layers_sorted, 2)-1

    I[perm[1:lower_idx]] = I_0

    for layer in 2:max_layer
        lower_idx = searchsortedfirst(layers_sorted, layer)
        upper_idx = searchsortedfirst(layers_sorted, layer+1)
        for sweep in 1:n_sweeps
            for i in upper_idx-1:-1:lower_idx
                # coordinate
                idx=perm[i]
                position = sites.positions[:,idx]

                # number of neighbours
                n_neighbours = sites.neighbours[idx, 1]
                neighbours = sites.neighbours[idx, 2:n_neighbours+1]

                dot_products, ∠_indices, upwind_positions = smallest_angle(position, neighbours, -k, sites, 2)

                dot_weights = [dot_products[1]^p, dot_products[2]^p]./sum(dot_products.^p)

                I[idx] = 0u"kW*nm^-1*m^-2"
                for rn in 1:2
                    #∠_index = choose_random(smallest∠, ∠_indices)
                    ∠_index = rn
                    upwind_index = neighbours[∠_indices[∠_index]]
                    upwind_position = upwind_positions[:, ∠_index]


                    # Find the intersection and the neighbour the ray is coming from
                    # upwind_position, upwind_index = ray_intersection(k, neighbours, position, sites, i)

                    # Pass α as an array
                    α_centre = α[idx]
                    α_upwind = α[upwind_index]

                    r = euclidean(position*1.0, upwind_position)
                    # Find the Δτ optical path from upwind to grid point
                    Δτ_upwind = trapezoidal(r, α_centre, α_upwind)

                    S_centre = S[idx]
                    S_upwind = S[upwind_index]

                    w1, w2 =  weights(Δτ_upwind)
                    a_ijk, b_ijk, c_ijk = coefficients(w1, w2, Δτ_upwind)

                    I_upwind = I[upwind_index]
                    I[idx] += (a_ijk*S_upwind + b_ijk*S_centre + c_ijk*I_upwind)*dot_weights[rn]
                    # maybe try a mix of this and the ray intersection?
                    # Need to think more
                end
            end
        end
    end
    return I
end

"""
Make a method that combines the ray intersection and the Delaunay lines. Use
Delaunay lines to find quantities, and ray intersection to find lenght the ray
travels. Mix with MC method (??).
"""

function Delaunay_ray_up(sites::VoronoiSites,
                    I_0::AbstractMatrix, S_0::AbstractMatrix, α_0::AbstractMatrix,
                    S::AbstractVector, α::AbstractVector, k::AbstractVector, n_sweeps::Int, Nran::Int)

    #
    # Inverse distance power law parameter
    p = 3.0

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

    perm = sortperm(sites.layers_up)
    layers_sorted = sites.layers_up[perm]
    max_layer = layers_sorted[end]

    # Calculate intensity in boundary
    lower_boundary(I, sites, I_0, S_0, α_0, S, α, k)

    # Find out which cells passes radiation to cells not in the same layer... (a good start)

    for layer in 2:max_layer
        # Calculate radiation from bottom layer upwards
        lower_idx = searchsortedfirst(layers_sorted, layer)
        upper_idx = searchsortedfirst(layers_sorted, layer+1)

        # Keep track of completed sites
        completed_cells = zeros(upper_idx - lower_idx - 1)

        for i in 1:upper_idx - lower_idx - 1
            idx = perm[i + lower_idx]

            position = sites.positions[:, idx]

            # number of neighbours
            n_neighbours = sites.neighbours[idx, 1]
            neighbours = sites.neighbours[idx, 2:n_neighbours+1]

            # Check if lower cells are solved
            dot_products, ∠_indices, upwind_positions = smallest_angle(position, neighbours, -k, sites, 2)

            dot_weights = [dot_products[1], dot_products[2]]./(sum(dot_products))

            if layers_sorted[n1] < layer && layers_sorted[n2] < layer
                # Calculate radiation propagating to neighbouring cells
                for rn in 1:2
                    ∠_index = rn
                    upwind_position = upwind_positions[:, ∠_index]
                    upwind_idx = neighbours[∠_indices[∠_index]]

                    # Find the intersection and the neighbour the ray is coming from
                    # upwind_position, upwind_index = ray_intersection(k, neighbours, position, sites, i)

                    # Pass α as an array
                    α_centre = α[idx]
                    α_upwind = α[upwind_idx]

                    r = euclidean(position, upwind_position)
                    # Find the Δτ optical path from upwind to grid point
                    Δτ_upwind = trapezoidal(r, α_centre, α_upwind)

                    S_centre = S[idx]
                    S_upwind = S[upwind_idx]

                    w1, w2 =  weights(Δτ_upwind)
                    a_ijk, b_ijk, c_ijk = coefficients(w1, w2, Δτ_upwind)

                    I_upwind = I[upwind_idx]
                    I[idx] += (a_ijk*S_upwind + b_ijk*S_centre + c_ijk*I_upwind)/dot_weights[rn]

                end
                completed_cells[i] = 1
            end
        end

        if sum(completed_cells) == upper_idx - lower_idx
            break
        end

        for i in arg_where(completed_cells, 0)
            idx = perm[i + lower_idx]

            position = sites.positions[:, idx]

            # number of neighbours
            n_neighbours = sites.neighbours[idx, 1]
            neighbours = sites.neighbours[idx, 2:n_neighbours+1]

            # Check if lower cells are solved
            smallest∠, ∠_indices, upwind_positions = smallest_angle(position, neighbours, -k, sites, 2)

            n1 = neighbours[∠_indices[1]]
            n2 = neighbours[∠_indices[2]]

            if layers_sorted[n2] < layer && layers_sorted[n1] >= layer
                # Calculate radiation propagating to neighbouring cells
                upwind_idx = n2
                upwind_position = upwind_positions[:, 2]

                # Find the intersection and the neighbour the ray is coming from
                # upwind_position, upwind_index = ray_intersection(k, neighbours, position, sites, i)

                # Pass α as an array
                α_centre = α[idx]
                α_upwind = α[upwind_idx]

                r = euclidean(position, upwind_position)
                # Find the Δτ optical path from upwind to grid point
                Δτ_upwind = trapezoidal(r, α_centre, α_upwind)

                S_centre = S[idx]
                S_upwind = S[upwind_idx]

                w1, w2 =  weights(Δτ_upwind)
                a_ijk, b_ijk, c_ijk = coefficients(w1, w2, Δτ_upwind)

                I_upwind = I[upwind_idx]
                I[idx] += (a_ijk*S_upwind + b_ijk*S_centre + c_ijk*I_upwind)

                completed_cells[i] = 1
            elseif layers_sorted[n1] < layer && layers_sorted[n2] >= layer
                # Calculate radiation propagating to neighbouring cells
                upwind_idx = n1
                upwind_position = upwind_positions[:, 1]

                # Find the intersection and the neighbour the ray is coming from
                # upwind_position, upwind_index = ray_intersection(k, neighbours, position, sites, i)

                # Pass α as an array
                α_centre = α[idx]
                α_upwind = α[upwind_idx]

                r = euclidean(position, upwind_position)
                # Find the Δτ optical path from upwind to grid point
                Δτ_upwind = trapezoidal(r, α_centre, α_upwind)

                S_centre = S[idx]
                S_upwind = S[upwind_idx]

                w1, w2 =  weights(Δτ_upwind)
                a_ijk, b_ijk, c_ijk = coefficients(w1, w2, Δτ_upwind)

                I_upwind = I[upwind_idx]
                I[idx] += (a_ijk*S_upwind + b_ijk*S_centre + c_ijk*I_upwind)

                completed_cells[i] = 1
            end
        end

        if sum(completed_cells) == upper_idx - lower_idx
            break
        end

        for sweep in n_sweeps
            for i in arg_where(completed_cells, 0)
                idx = perm[i + lower_idx]

                I[idx] = 0

                position = sites.positions[:, idx]

                # number of neighbours
                n_neighbours = sites.neighbours[idx, 1]
                neighbours = sites.neighbours[idx, 2:n_neighbours+1]

                # Check if lower cells are solved
                smallest∠, ∠_indices, upwind_positions = smallest_angle(position, neighbours, -k, sites, 2)

                for rn in 1:Nran
                    ∠_index = choose_random(smallest∠, ∠_indices)
                    upwind_position = upwind_positions[:, ∠_index]
                    upwind_idx = neighbours[∠_indices[∠_index]]


                    # Find the intersection and the neighbour the ray is coming from
                    # upwind_position, upwind_index = ray_intersection(k, neighbours, position, sites, i)

                    # Pass α as an array
                    α_centre = α[idx]
                    α_upwind = α[upwind_idx]

                    r = euclidean(position, upwind_position)
                    # Find the Δτ optical path from upwind to grid point
                    Δτ_upwind = trapezoidal(r, α_centre, α_upwind)

                    S_centre = S[idx]
                    S_upwind = S[upwind_idx]

                    w1, w2 =  weights(Δτ_upwind)
                    a_ijk, b_ijk, c_ijk = coefficients(w1, w2, Δτ_upwind)

                    I_upwind = I[upwind_idx]
                    I[idx] += (a_ijk*S_upwind + b_ijk*S_centre + c_ijk*I_upwind)/Nran
                end
            end
        end
    end

    return I
end

function lower_boundary(I, sites::VoronoiSites,
                    I_0::AbstractMatrix, S_0::AbstractMatrix, α_0::AbstractMatrix,
                    S::AbstractVector, α::AbstractVector, k::AbstractVector)

    x_min = sites.x_min
    x_max = sites.x_max
    y_min = sites.y_min
    y_max = sites.y_max

    nx = length(I_0[:,1])
    ny = length(I_0[1,:])

    Δx = (x_max - x_min)/(nx - 1)
    Δy = (y_max - y_min)/(ny - 1)

    perm = sortperm(sites.layers_up)
    layers_sorted = sites.layers_up[perm]
    max_layer = layers_sorted[end]

    z_upwind = sites.z_min

    lower_idx = searchsortedfirst(layers_sorted, 1)
    upper_idx = searchsortedfirst(layers_sorted, 2)
    for i in lower_idx:upper_idx-1
        # coordinate
        idx = perm[i]
        position = sites.positions[:,idx]

        # 1st layer
        Δz = position[1] - z_upwind
        r = abs(Δz/k[1])

        x_upwind = position[2] - r*k[2]
        if x_upwind > sites.x_max
            x_upwind = sites.x_min + (x_upwind - sites.x_max)
        elseif x_upwind < sites.x_min
            x_upwind = sites.x_max + (x_upwind - sites.x_min)
        end

        y_upwind = position[3] - r*k[3]
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
        α_centre = α[idx]
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
        I[idx] = a_ijk*S_upwind + b_ijk*S_centre + c_ijk*I_upwind
    end
end

function upper_boundary(I, sites::VoronoiSites,
                    I_0::AbstractMatrix, S_0::AbstractMatrix, α_0::AbstractMatrix,
                    S::AbstractVector, α::AbstractVector, k::AbstractVector)


    #
    x_min = sites.x_min
    x_max = sites.x_max
    y_min = sites.y_min
    y_max = sites.y_max

    nx = length(I_0[:,1])
    ny = length(I_0[1,:])

    Δx = (x_max - x_min)/(nx - 1)
    Δy = (y_max - y_min)/(ny - 1)

    perm = sortperm(sites.layers_down)
    layers_sorted = sites.layers_down[perm]
    max_layer = layers_sorted[end]

    z_upwind = sites.z_max

    lower_idx = searchsortedfirst(layers_sorted, 1)
    upper_idx = searchsortedfirst(layers_sorted, 2)
    for i in lower_idx:upper_idx-1
        # coordinate
        idx=perm[i]
        position = sites.positions[:,idx]

        # 1st layer
        Δz = z_upwind - position[1]
        r = abs(Δz/k[1])

        x_upwind = position[2] - r*k[2]
        if x_upwind > sites.x_max
            x_upwind = sites.x_min + (x_upwind - sites.x_max)
        elseif x_upwind < sites.x_min
            x_upwind = sites.x_max + (x_upwind - sites.x_min)
        end

        y_upwind = position[3] - r*k[3]
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
        α_centre = α[idx]
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
        I[idx] = a_ijk*S_upwind + b_ijk*S_centre + c_ijk*I_upwind
    end
end
