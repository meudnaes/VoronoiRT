using NearestNeighbors

include("functions.jl")

struct VoronoiSites
    positions::Matrix{Float64}
    temperature::Vector{Float64}
    electron_density::Vector{Float64}
    hydrogen_populations::Vector{Float64}
    z_min::Float64
    z_max::Float64
    x_min::Float64
    x_max::Float64
    y_min::Float64
    y_max::Float64
end

struct VoronoiCell
    ID::Int
    position::Vector{Float64}
    volume::Float64 # Volume of the cell
    neighbours::Vector{Int}
    faces::Vector{Float64}  # Area of the faces
    n::Int # Number of neighbours
end

struct RasterDomain
    z::Vector{<:Unitful.Length}
    x::Vector{<:Unitful.Length}
    y::Vector{<:Unitful.Length}
end

"""
    read_neighbours(fname::String, n_sites::Int)

Reads a file containing neighbouring cells for each grid point in the voronoi
tesselation.
"""
function read_cell(fname::String, n_sites::Int, sites::VoronoiSites)
    ID = Vector{Int64}(undef, n_sites)
    cells = Vector{VoronoiCell}(undef, n_sites)
    open(fname, "r") do io
        for (i, l) in enumerate(eachline(io))
            lineLength = length(split(l))

            N = Int((lineLength-2)/2)

            # Which cell
            ID[i] = parse(Int64, split(l)[1])

            # Coordinates
            z = sites.positions[1, ID[i]]
            x = sites.positions[2, ID[i]]
            y = sites.positions[3, ID[i]]

            # Volume of cell
            volume = parse(Float64, split(l)[2])

            # Neighbouring cells
            neighbours = Vector{Int}(undef, N)
            for j in 3:N+2
                neighbours[j-2] = parse(Int64, split(l)[j])
            end

            # Area of faces
            area = Vector{Float64}(undef, N)
            for j in N+3:lineLength
                area[j-(N+2)] = parse(Float64, split(l)[j])
            end

            # Store all neighbour information for current cell
            cells[i] = VoronoiCell(ID[i], [z, x, y], volume, neighbours, area, N)
        end
    end

    layers = _sort_by_layer(n_sites, cells)


    # Sort the arrays
    p = sortperm(layers)
    layers = layers[p]
    cells = cells[p, :]

    return cells, layers

end

function _sort_by_layer(n_sites::Int, cells::AbstractArray{VoronoiCell})
    layers = zeros(Int, n_sites)

    lower_boundary = -5
    for cell in cells
        if any(i -> i==lower_boundary, cell.neighbours)
            layers[cell.ID] = 1
        end
    end

    lower_layer=1
    while true
        for cell in cells
            if layers[cell.ID] == 0
                for neighbor_ID in cell.neighbours
                    if neighbor_ID > 0 && layers[neighbor_ID] == lower_layer
                        layers[cell.ID] = lower_layer+1
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

function inv_dist_itp(idxs, dists, p, sites::VoronoiSites)
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

function inv_dist_itp_test(idxs, dists, p, values)
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
