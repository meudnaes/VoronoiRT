using NearestNeighbors

struct VoronoiSites
    z::Vector{Float64}
    x::Vector{Float64}
    y::Vector{Float64}
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
    z::Float64
    x::Float64
    y::Float64
    volume::Float64 # Volume of the cell
    neighbours::Vector{Int}
    faces::Vector{Float64}  # area of the faces
    n::Int
end

struct RasterDomain
    z::Vector{<:Unitful.Length}
    x::Vector{<:Unitful.Length}
    y::Vector{<:Unitful.Length}
end

function sort_array(array::AbstractArray; axis=1)::AbstractArray
    ix = sortperm(array[axis, :])
    array_sorted = array[:, ix]
end

"""
    read_neighbours(fname::String, n_sites::Int)

Reads a file containing neighbouring cells for each grid point in the voronoi
tesselation.
"""
function read_cell(fname::String, n_sites::Int, sites::VoronoiSites)
    ID = Vector{Int64}(undef, n_sites)
    cell = Vector{VoronoiCell}(undef, n_sites)
    open(fname, "r") do io
        for (i, l) in enumerate(eachline(io))
            lineLength = length(split(l))

            N = Int((lineLength-2)/2)

            # Which cell
            ID[i] = parse(Int64, split(l)[1])
            z = sites.z[ID[i]]
            x = sites.x[ID[i]]
            y = sites.y[ID[i]]
            # Volume of cell
            volume = parse(Float64, split(l)[2])

            # neighbouring cells
            neighbours = Vector{Int}(undef, N)
            for j in 3:N+2
                neighbours[j-2] = parse(Int64, split(l)[j])
            end

            area = Vector{Float64}(undef, N)
            # area of faces
            for j in N+3:lineLength
                area[j-(N+2)] = parse(Float64, split(l)[j])
            end

            # store all neighbour information
            cell[i] = VoronoiCell(ID[i], z, x, y, volume, neighbours, area, N)
        end
    end
    return cell[sortperm(ID), :]
end

function inv_dist_itp(idxs, dists, p, sites::VoronoiSites)
    avg_inv_dist = 0
    f = 0
    for i in 1:length(idxs)
        idx = idxs[i]
        inv_dist = 1/dists[i]^p
        avg_inv_dist += inv_dist
        f += sites.hydrogen_populations[idx]*inv_dist
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
