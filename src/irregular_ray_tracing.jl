include("voronoi_utils.jl")
include("functions.jl")

using Distances
using LinearAlgebra

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
