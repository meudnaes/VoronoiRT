using NearestNeighbors

include("functions.jl")

struct VoronoiSites
    positions::Matrix{Float64}
    neighbours::Matrix{Int}
    layers_up::Vector{Int}
    layers_down::Vector{Int}
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

    global NeighbourMatrix
    NeighbourMatrix = NeighbourMatrix[:,1:max_neighbours+1]
    layers_up = _sort_by_layer_up(NeighbourMatrix, n_sites)
    layers_down = _sort_by_layer_down(NeighbourMatrix, n_sites)

    # Sorting works for layers and sites, but loses neighbour information!
    return positions, NeighbourMatrix, layers_up, layers_down
end

function _sort_by_layer_up(neighbours::AbstractMatrix, n_sites::Int)

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

function _sort_by_layer_down(neighbours::AbstractMatrix, n_sites::Int)
    layers = zeros(Int, n_sites)

    upper_boundary = -6
    for i in 1:n_sites
        n_neighbours = neighbours[i,1]
        for j in 1:n_neighbours
            if neighbours[i,j+1] == upper_boundary
                layers[i] = 1
            end
        end
    end

    upper_layer=1
    while true
        for i in 1:n_sites
            if layers[i] == 0
                n_neighbours = neighbours[i,1]
                for j in 1:n_neighbours
                    neighbour = neighbours[i,j+1]
                    if neighbour > 0 && layers[neighbour] == upper_layer
                        layers[i] = upper_layer+1
                        break
                    end
                end
            end
        end

        if !any(i -> i==0, layers)
            break
        end

        upper_layer += 1
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

function beam(v::AbstractVector, v0::AbstractVector, R0::Float64, k::AbstractVector)
    r = abs(v[1] - v0[1])/k[1]

    vh = v0 .+ r .* k

    if sqrt((v[2] - vh[2])^2 + (v[3] - vh[3])^2) < R0
        return 1
    else
        return 0.1
    end
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

function sample_beam(n_sites::Int, boundaries::Matrix, func, v0::AbstractVector, R0::Float64, k::AbstractVector)
    # Find max and min to convert random number between 0 and 1 to coordinate
    println("---Sampling new sites---")
    z_min = boundaries[1,1]; z_max = boundaries[1,2]
    x_min = boundaries[2,1]; x_max = boundaries[2,2]
    y_min = boundaries[3,1]; y_max = boundaries[3,2]

    Δz = z_max - z_min
    Δx = x_max - x_min
    Δy = y_max - y_min

    # allocate arrays for new sites
    p_vec = Matrix{Float64}(undef, (3, n_sites))

    for i in 1:n_sites
        print("site $i/$n_sites \r")
        while true
            ref_vec = rand(Float64, 3)
            z_ref = ref_vec[1]*Δz + z_min
            x_ref = ref_vec[2]*Δx + x_min
            y_ref = ref_vec[3]*Δy + y_min

            # acceptance criterion, "reference"
            density_ref = beam(ref_vec, v0, R0, k)
            # random sample, compare to reference
            density_ran = rand(Float64)
            if density_ref > density_ran
                # a point is accepted, store position and move on
                p_vec[:, i] .= (z_ref, x_ref, y_ref)
                # break to find next site
                break
            end
        end
    end
    print("                                                                 \r")
    return p_vec
end

function upwind_distances(neighbours::Vector{Int}, n_neighbours::Integer, id::Integer,
                upwind_position::Vector{Float64}, sites::VoronoiSites)

    p_r = sites.positions[:, id]

    distances = Vector{Float64}(undef, n_neighbours)

    # This has length = cell.n - 1
    cell_neighbours = neighbours[neighbours .> 0]

    x_r_r = sites.x_max - p_r[2]
    x_r_l = p_r[2] - sites.x_min

    y_r_r = sites.y_max - p_r[3]
    y_r_l = p_r[3] - sites.y_min

    for (j, neighbour) in enumerate(cell_neighbours)
        if neighbour > 0
            p_n = sites.positions[:, neighbour]

            x_n_r = sites.x_max - p_n[2]
            x_n_l = p_n[2] - sites.x_min

            # Test for periodic
            if x_r_r + x_n_l < p_r[2] - p_n[2]
                p_n[2] = sites.x_max + p_n[2] - sites.x_min
            elseif x_r_l + x_n_r < p_n[2] - p_r[2]
                p_n[2] = sites.x_min + sites.x_max - p_n[2]
            end

            y_n_r = sites.y_max - p_n[3]
            y_n_l = p_n[3] - sites.y_min

            # Test for periodic
            if y_r_r + y_n_l < p_r[3] - p_n[3]
                p_n[3] = sites.y_max + p_n[3] - sites.y_min
            elseif y_r_l + y_n_r < p_n[3] - p_r[3]
                p_n[3] = sites.y_min + sites.y_max - p_n[3]
            end
            distances[j] = euclidean(upwind_position, p_n)
        elseif neighbour == -5
            continue
        elseif neighbour == -6
            continue
        end
    end

    distances[n_neighbours] = euclidean(upwind_position, p_r)
    return distances
end

function ray_intersection(k::AbstractArray, neighbours::Vector{Int}, position::Vector{Float64},
                          sites::VoronoiSites, ID)
    # k is unit vector towards upwind direction of the ray

    p_r = position

    n_neighbours = length(neighbours)

    x_r_r = abs(sites.x_max - p_r[2])
    x_r_l = abs(p_r[2] - sites.x_min)

    y_r_r = abs(sites.y_max - p_r[3])
    y_r_l = abs(p_r[3] - sites.y_min)

    s = Vector{Float64}(undef, n_neighbours)
    for (i, neighbour) in enumerate(neighbours)
        if neighbour > 0
            # Natural neighbor position
            p_i = sites.positions[:, neighbour]

            x_i_r = abs(sites.x_max - p_i[2])
            x_i_l = abs(p_i[2] - sites.x_min)

            # Test for periodic
            if x_r_r + x_i_l < p_r[2] - p_i[2]
                p_i[2] = sites.x_max + p_i[2] - sites.x_min
            elseif x_r_l + x_i_r < p_i[2] - p_r[2]
                p_i[2] = sites.x_min + sites.x_max - p_i[2]
            end

            y_i_r = abs(sites.y_max - p_i[3])
            y_i_l = abs(p_i[3] - sites.y_min)

            # Test for periodic
            if y_r_r + y_i_l < p_r[3] - p_i[3]
                p_i[3] = sites.y_max + p_i[3] - sites.y_min
            elseif y_r_l + y_i_r < p_i[3] - p_r[3]
                p_i[3] = sites.y_min + sites.y_max - p_i[3]
            end

            # Calculate plane bisecting site and neighbor
            # Normal vector
            n = p_i .- p_r

            # Point on the plane
            p = (p_i .+ p_r) ./ 2
            s[i] = dot(n, p .- p_r)/dot(n, k)
        elseif neighbour == -6
            # Boundary
            z_upwind = sites.z_max
            z = p_r[1]
            s[i] = (z_upwind - z)/k[1]
        elseif neighbour == -5
            # Boundary
            z_upwind = sites.z_min
            z = p_r[1]
            s[i] = (z_upwind - z)/k[1]
        end
    end



    # Find the smallst non-negative s
    index, s_q = smallest_non_negative(s)
    q = p_r .+ s_q .* k

    return q, neighbours[index]
end

function smallest_angle(position::AbstractVector, neighbours::AbstractVector, k::AbstractVector, sites::VoronoiSites, n::Int)

    dots = Vector{Float64}(undef, length(neighbours))
    upwind_positions = Matrix{Float64}(undef, (3, length(neighbours)))

    x_r_r = sites.x_max - position[2]
    x_r_l = position[2] - sites.x_min

    y_r_r = sites.y_max - position[3]
    y_r_l = position[3] - sites.y_min

    for i in 1:length(neighbours)
        neighbour = neighbours[i]
        if neighbour > 0
            p_n = sites.positions[:, neighbour]

            x_i_r = abs(sites.x_max - p_n[2])
            x_i_l = abs(p_n[2] - sites.x_min)

            # Test for periodic
            if x_r_r + x_i_l < position[2] - p_n[2]
                p_n[2] = sites.x_max + p_n[2] - sites.x_min
            elseif x_r_l + x_i_r < p_n[2] - position[2]
                p_n[2] = sites.x_min + sites.x_max - p_n[2]
            end

            y_i_r = abs(sites.y_max - p_n[3])
            y_i_l = abs(p_n[3] - sites.y_min)

            # Test for periodic
            if y_r_r + y_i_l < position[3] - p_n[3]
                p_n[3] = sites.y_max + p_n[3] - sites.y_min
            elseif y_r_l + y_i_r < p_n[3] - position[3]
                p_n[3] = sites.y_min + sites.y_max - p_n[3]
            end

            direction = position .- p_n
            norm_dir = direction/(norm(direction))
            # Two normalized direction vectors, denominator is 1
            # Save the dot product, don't bother calculating the angle
            dots[i] = dot(k, norm_dir)

            # Save the upwind positions
            upwind_positions[:, i] = p_n
        else
            dots[i] = -1
        end
    end

    p=sortperm(dots)
    dots=dots[p]
    upwind_positions=upwind_positions[:, p]

    return dots[end-(n-1):end], p[end-(n-1):end], upwind_positions[:, end-(n-1):end]
end

function choose_random(angles::AbstractVector, indices::AbstractVector)

    p2 = rand()

    r1 = acos(angles[2])/acos(angles[1])
    p1 = r1*rand()

    return argmax([p1, p2])
end
