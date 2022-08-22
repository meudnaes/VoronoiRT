include("voronoi_utils.jl")
include("functions.jl")

using Distances
using LinearAlgebra

const p = 7.0           # Weighting upwind rays

"""
    Delaunay_upII(k::Vector{Float64},
                       S::Vector{<:UnitsIntensity_λ},
                       I_0::Vector{<:UnitsIntensity_λ},
                       α::Vector{<:PerLength},
                       sites::VoronoiSites,
                       n_sweeps::Int)

Computes intensity along rays traveling upwards from the bottom to the top
through all grid points in the atmosphere, in an irregular grid. Initial
intensity I_0 at the bottom of the domain, usually B_λ.
"""
function Delaunay_upII(k::Vector{Float64},
                       S::Vector{<:UnitsIntensity_λ},
                       I_0::Vector{<:UnitsIntensity_λ},
                       α::Vector{<:PerLength},
                       sites::VoronoiSites,
                       n_sweeps::Int)

    # Allocate space for intensity
    I = zero(S)

    x_min = sites.x_min
    x_max = sites.x_max
    y_min = sites.y_min
    y_max = sites.y_max

    perm = sites.perm_up
    max_layer = length(sites.layers_up)

    lower_idx = sites.layers_up[2]-1

    I[perm[1:lower_idx]] = I_0

    for layer in 2:max_layer-1
        lower_idx = sites.layers_up[layer]
        upper_idx = sites.layers_up[layer+1]
        for sweep in 1:n_sweeps
            for i in lower_idx:upper_idx-1
                # coordinate
                idx = perm[i]
                position = sites.positions[:,idx]

                # number of neighbours
                n_neighbours = sites.neighbours[idx, 1]
                neighbours = sites.neighbours[idx, 2:n_neighbours+1]

                dot_products, upwind_indices = smallest_angle(idx, neighbours, k, sites)
                dot_weights = [dot_products[1]^p, dot_products[2]^p]./sum(dot_products.^p)

                I[idx] = 0u"kW*nm^-1*m^-2"
                for rn in 1:2
                    upwind_index = upwind_indices[rn]
                    upwind_position = sites.positions[:, upwind_index]


                    # Find the intersection and the neighbour the ray is coming from
                    # upwind_position, upwind_index = ray_intersection(k, neighbours, position, sites, i)

                    # Pass α as an array
                    α_centre = α[idx]
                    α_upwind = α[upwind_index]

                    r = euclidean(position*1.0, upwind_position)#*dot_products[rn]
                    # Find the Δτ optical path from upwind to grid point
                    Δτ_upwind = trapezoidal(r, α_centre, α_upwind)

                    S_centre = S[idx]
                    S_upwind = S[upwind_index]

                    α_ijk, β_ijk, expΔτ = linear_weights(Δτ_upwind)

                    I_upwind = I[upwind_index]
                    I[idx] += (expΔτ*I_upwind + α_ijk*S_upwind + β_ijk*S_centre)*dot_weights[rn]
                end
            end
        end
    end
    return I
end

"""
    Delaunay_downII(k::Vector{Float64},
                       S::Vector{<:UnitsIntensity_λ},
                       I_0::Vector{<:UnitsIntensity_λ},
                       α::Vector{<:PerLength},
                       sites::VoronoiSites,
                       n_sweeps::Int)

Computes intensity along rays traveling downwards from the top to the bottom
through all grid points in the atmosphere, in an irregular grid. Initial
intensity I_0 at the bottom of the domain, usually 0.
"""
function Delaunay_downII(k::Vector{Float64},
                         S::Vector{<:UnitsIntensity_λ},
                         I_0::Vector{<:UnitsIntensity_λ},
                         α::Vector{<:PerLength},
                         sites::VoronoiSites,
                         n_sweeps::Int)

    # Allocate space for intensity
    I = zero(S)

    x_min = sites.x_min
    x_max = sites.x_max
    y_min = sites.y_min
    y_max = sites.y_max

    perm = sites.perm_down
    max_layer = length(sites.layers_down)

    lower_idx = sites.layers_down[2] - 1

    I[perm[1:lower_idx]] = I_0

    for layer in 2:max_layer-1
        lower_idx = sites.layers_down[layer]
        upper_idx = sites.layers_down[layer+1]
        for sweep in 1:n_sweeps
            for i in upper_idx-1:-1:lower_idx
                # coordinate
                idx=perm[i]
                position = sites.positions[:,idx]

                # number of neighbours
                n_neighbours = sites.neighbours[idx, 1]
                neighbours = sites.neighbours[idx, 2:n_neighbours+1]

                dot_products, upwind_indices = smallest_angle(idx, neighbours, k, sites)
                dot_weights = [dot_products[1]^p, dot_products[2]^p]./sum(dot_products.^p)

                I[idx] = 0u"kW*nm^-1*m^-2"
                for rn in 1:2
                    upwind_index = upwind_indices[rn]
                    upwind_position = sites.positions[:, upwind_index]


                    # Find the intersection and the neighbour the ray is coming from
                    # upwind_position, upwind_index = ray_intersection(k, neighbours, position, sites, i)

                    # Pass α as an array
                    α_centre = α[idx]
                    α_upwind = α[upwind_index]

                    r = euclidean(position*1.0, upwind_position)#*dot_products[rn]
                    # Find the Δτ optical path from upwind to grid point
                    Δτ_upwind = trapezoidal(r, α_centre, α_upwind)

                    S_centre = S[idx]
                    S_upwind = S[upwind_index]

                    α_ijk, β_ijk, expΔτ = linear_weights(Δτ_upwind)

                    I_upwind = I[upwind_index]
                    I[idx] += (expΔτ*I_upwind + α_ijk*S_upwind + β_ijk*S_centre)*dot_weights[rn]
                end
            end
        end
    end
    return I
end
