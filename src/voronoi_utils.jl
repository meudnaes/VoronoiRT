using NearestNeighbors

include("functions.jl")

struct VoronoiSites
    positions::Matrix{Float64}
    neighbours::Matrix{Int}
    layers::Vector{Int}
    temperature::Vector{Float64}
    electron_density::Vector{Float64}
    hydrogen_populations::Vector{Float64}
    z_min::Float64
    z_max::Float64
    x_min::Float64
    x_max::Float64
    y_min::Float64
    y_max::Float64
    n::Int
end

struct VoronoiCell
    ID::Int
    position::Vector{Float64}
    volume::Float64 # Volume of the cell
    neighbours::Vector{Int}
    faces::Vector{Float64}  # Area of the faces
    n::Int # Number of neighbours
end

"""
    read_neighbours(fname::String, n_sites::Int)

Reads a file containing neighbouring cells for each grid point in the voronoi
tesselation.
"""
function read_cell(fname::String, n_sites::Int, positions::AbstractMatrix)
    # Guess (overshoot), maybe do exact later?
    max_guess = 100
    NeighbourMatrix = zeros(Int, n_sites, max_guess+1)
    open(fname, "r") do io
        for (i, l) in enumerate(eachline(io))
            lineLength = length(split(l))

            N = Int(lineLength - 1)

            # Which cell
            ID = parse(Int64, split(l)[1])

            # Neighbouring cells
            neighbours = Vector{Int}(undef, N)
            for j in 1:N
                neighbours[j] = parse(Int64, split(l)[j+1])
            end

            # Store all neighbour information for current cell
            NeighbourMatrix[ID, 1] = N
            NeighbourMatrix[ID, 2:N+1] = neighbours
        end
    end
    max_neighbours = maximum(NeighbourMatrix[:,1])
    if max_neighbours == max_guess
        println("Guess too low!")
    end
    NeighbourMatrix = NeighbourMatrix[:,1:max_neighbours+1]
    layers = _sort_by_layer(NeighbourMatrix, n_sites)

    # Sorting works for layers and sites, but loses neighbour information!
    p = sortperm(layers)
    positions = positions[:,p]
    layers = layers[p]
    NeighbourMatrix = _sort_neighbours(NeighbourMatrix, p)
    return positions, NeighbourMatrix, layers
end

function _sort_by_layer(neighbours::AbstractMatrix, n_sites::Int)

    layers = zeros(Int, n_sites)

    lower_boundary = -5
    for i in 1:n_sites
        n_neighbours = neighbours[i,1]
        for j in 1:n_neighbours
            if neighbours[i,j+1] == lower_boundary
                layers[i] = 1
            end
        end
    end

    lower_layer=1
    while true
        for i in 1:n_sites
            if layers[i] == 0
                n_neighbours = neighbours[i,1]
                for j in 1:n_neighbours
                    neighbour = neighbours[i,j+1]
                    if neighbour > 0 && layers[neighbour] == lower_layer
                        layers[i] = lower_layer+1
                        break
                    end
                end
            end
        end

        if !any(i -> i==0, layers)
            break
        end

        lower_layer += 1
    end
    return layers
end

function _sort_neighbours(NeighbourMatrix::AbstractMatrix, permVec::AbstractVector)

    sortedNeighbours = zero(NeighbourMatrix)
    sortedIdx = sortperm(permVec)
    for (i, idp) in enumerate(permVec)
        n = NeighbourMatrix[idp, 1]
        neighbours = NeighbourMatrix[idp, 2:n+1]
        for j in 1:n
            neighbour = neighbours[j]
            if neighbour > 0
                neighbours[j] = sortedIdx[neighbour]
            end
        end
        sortedNeighbours[i, 1] = n
        sortedNeighbours[i, 2:n+1] = neighbours
    end
    return sortedNeighbours
end

function inv_dist_itp_test(idxs, dists, p, sites::VoronoiSites)
    avg_inv_dist = 0
    f = 0
    for i in 1:length(idxs)
        idx = idxs[i]
        inv_dist = 1/dists[i]^p
        avg_inv_dist += inv_dist
        f += values[idx]*inv_dist
    end
    f = f/avg_inv_dist
end

function inv_dist_itp(idxs::AbstractVector, dists::AbstractVector,
                      p::AbstractFloat, values::AbstractVector)
    avg_inv_dist = 0
    f = 0
    for i in 1:length(idxs)
        idx = idxs[i]
        if idx > 0
            inv_dist = 1/dists[i]^p
            avg_inv_dist += inv_dist
            f += values[idx]*inv_dist
        elseif idx == -5
            # lower boundary
            f+=0
        elseif idx == -6
            # upper boundary
            f+=0
        end
    end
    f = f/avg_inv_dist
end
