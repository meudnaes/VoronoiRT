include("voronoi_utils.jl")

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

    # 1st layer
    z_upwind = sites.z_min

    # precalculate trigonometric functions
    cosθ = cos(θ)
    sinθ = sin(θ)

    cosϕ = cos(ϕ)
    sinϕ = sin(ϕ)

    # Unit vector towards upwind direction of the ray
    k = -[cosθ, cosϕ*sinθ, sinϕ*sinθ]

    # Allocate space for intensity
    I = zero(S_0)

    for layer in 1:maximum(layers)
        lower_idx = searchsortedfirst(layers, layer)
        upper_idx = searchsortedfirst(layers, layer+1)
        for cell in cells
            # coordinate
            position = cell.position

            # Find the intersection and the neighbour the ray is coming from
            upwind_index, upwind_position = _rayIntersection(k, cell, position)

            distances = Vector{Float64}(undef, cell.n+1)

            for (i, neighbour) in enumerate(cell.neighbours)
                distances[i] = euclidean(upwind_position, sites[neighbour])
            end
            distances[cell.n+1] = euclidean(upwind_position, position)

            # Pass α as an array
            α_centre = α_cont[cell.ID]
            α_upwind = inv_dist_itp([cell.neighbours..., cell.ID],
                                    distances, p, α_cont)

            # Find the Δτ optical path from upwind to grid point
            Δτ_upwind = trapezoidal(distances[cell.n+1], α_centre, α_upwind)

            S_centre = S_0[cell.ID]
            S_upwind = inv_dist_itp([cell.neighbours..., cell.ID],
                                    distances, p, S_0)

            w1, w2 =  weights(Δτ_upwind)
            a_ijk, b_ijk, c_ijk = coefficients(w1, w2, Δτ_upwind)

            # FIX!
            I_upwind = I[upwind_index]

            I[cell.ID] = a_ijk*S_upwind + b_ijk*S_centre + c_ijk*I_upwind
        end
    end
end


function _rayIntersection(k::AbstractArray, cell::VoronoiCell, sites::AbstractArray)
    # k is unit vector towards upwind direction of the ray

    r = [cell.z, cell.x, cell.y]
    s = Vector{Float64}(undef, cell.n)
    for i in 1:cell.n
        # Natural neighbor position
        p_i = sites[:,i]

        # Calculate plane biseting site and neighbor
        # Normal vector
        n = p_i - p_r

        # Point on the plane
        p = (p_i + p_r)/2

        s[i] = dot(n, p - r)/dot(n, k)
    end

    # Find the smallst non-negative s
    index, s_q = smallestNonNegative(s)

    q = r + s_q*k

    return index, q
end
