struct VoronoiSites
    z::Vector{<:Unitful.Length}
    x::Vector{<:Unitful.Length}
    y::Vector{<:Unitful.Length}
    temperature::Vector{<:Unitful.Temperature}
    electron_density::Vector{<:NumberDensity}
    hydrogen_populations::Vector{<:NumberDensity}
end

struct VoronoiCell
    ID::Int
    z::Unitful.Length
    x::Unitful.Length
    y::Unitful.Length
    volume::Volume
    neighbours::Vector{Int}
    n::Int
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
function read_neighbours(fname::String, n_sites::Int, sites::VoronoiSites)
    ID = Vector{Int64}(undef, n_sites)
    neighbours = Vector{VoronoiCell}(undef, n_sites)
    open(fname, "r") do io
        for (i, l) in enumerate(eachline(io))
            lineLength = length(split(l))

            # Which cell
            ID[i] = parse(Int64, split(l)[1])
            z = sites.z[ID[i]]
            x = sites.x[ID[i]]
            y = sites.y[ID[i]]
            # Volume of cell
            volume = parse(Float64, split(l)[2])*u"m^3"

            # neighbouring cells
            line = Vector{Int}(undef, lineLength-2)
            for j in 3:lineLength
                line[j-2] = parse(Int64, split(l)[j])
            end
            # store all neighbour information
            neighbours[i] = VoronoiCell(ID[i], z, x, y, volume, line, lineLength-2)
        end
    end
    return neighbours[sortperm(ID), :]
end
