using NearestNeighbors

include("functions.jl")

struct VoronoiSites
    positions::Matrix{<:Unitful.Length}
    neighbours::Matrix{Int}
    layers_up::Vector{Int}
    layers_down::Vector{Int}
    temperature::Vector{<:Unitful.Temperature}
    electron_density::Vector{<:NumberDensity}
    hydrogen_populations::Vector{<:NumberDensity}
    velocity_z::Vector{<:Unitful.Velocity}
    velocity_x::Vector{<:Unitful.Velocity}
    velocity_y::Vector{<:Unitful.Velocity}
    z_min::Unitful.Length
    z_max::Unitful.Length
    x_min::Unitful.Length
    x_max::Unitful.Length
    y_min::Unitful.Length
    y_max::Unitful.Length
    n::Int
end

"""
    read_neighbours(fname::String, n_sites::Int)

Reads a file containing neighbouring cells for each grid point in the voronoi
tesselation.
"""
function read_cell(fname::String, n_sites::Int, positions::AbstractMatrix)
    println("---Reading neighbour information---")
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
        println("===================Guess was too low!=======================")
    end

    global NeighbourMatrix
    NeighbourMatrix = NeighbourMatrix[:,1:max_neighbours+1]

    layers_up = _sort_by_layer_up(NeighbourMatrix, n_sites)
    ######################
    # Preparing to make code more efficient
    # Remember to make searchlight test work first!!!
    # I need the permutation vectors, and the sorted vectors.
    # Layer sorting is unneccessary after this function
    #######################

    # perm_up = sortperm(layers_up)
    # layers_sorted_up = layers_up[perm_up]

    layers_down = _sort_by_layer_down(NeighbourMatrix, n_sites)
    # perm_down = sortperm(layers_down)
    # layers_sorted = layers_down[perm_down]

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

function _sort_neighbours(NeighbourMatrix::AbstractMatrix,
                          permVec::AbstractVector)

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

function beam(v::AbstractVector,
              v0::AbstractVector,
              R0::Float64,
              k::AbstractVector)
    r = abs(v[1] - v0[1])/k[1]

    vh = v0 .+ r .* k

    if sqrt((v[2] - vh[2])^2 + (v[3] - vh[3])^2) < R0
        return 1
    else
        return 0.1
    end
end

function sample_beam(n_sites::Int,
                     boundaries::Matrix,
                     func,
                     v0::AbstractVector,
                     R0::Float64,
                     k::AbstractVector)
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

function upwind_distances(neighbours::Vector{Int},
                          n_neighbours::Integer,
                          id::Integer,
                          upwind_position::Vector{Float64},
                          sites::VoronoiSites)

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

function smallest_angle(position::Vector{<:Unitful.Length},
                        neighbours::Vector{Int},
                        k::Vector{Float64},
                        sites::VoronoiSites,
                        n::Int)

    dots = Vector{Float64}(undef, length(neighbours))
    upwind_positions = Matrix{Float64}(undef, (3, length(neighbours)))u"m"

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
    upwind_positions = upwind_positions[:, p]

    return dots[end-(n-1):end], p[end-(n-1):end], upwind_positions[:, end-(n-1):end]
end

function choose_random(angles::AbstractVector, indices::AbstractVector)

    p2 = rand()

    r1 = acos(angles[2])/acos(angles[1])
    p1 = r1*rand()

    return argmax([p1, p2])
end

function Voronoi_to_Raster(sites::VoronoiSites,
                           atmos::Atmosphere,
                           S_λ::Array{<:UnitsIntensity_λ},
                           α_tot::Array{<:PerLength},
                           populations::Array{<:NumberDensity},
                           r_factor)

    z = collect(LinRange(sites.z_min, sites.z_max,
                         floor(Int, r_factor*length(atmos.z))))
    x = collect(LinRange(sites.x_min, sites.x_max,
                         floor(Int, r_factor*length(atmos.x))))
    y = collect(LinRange(sites.y_min, sites.y_max,
                         floor(Int, r_factor*length(atmos.y))))

    nz = length(z)
    nx = length(x)
    ny = length(y)
    nλ = size(S_λ)[1]

    temperature = Array{Float64, 3}(undef, (nz, ny, nx))u"K"
    electron_density = Array{Float64, 3}(undef, (nz, ny, nx))u"m^-3"
    hydrogen_populations = Array{Float64, 3}(undef, (nz, ny, nx))u"m^-3"
    velocity_z = Array{Float64, 3}(undef, (nz, ny, nx))u"m*s^-1"
    velocity_x = copy(velocity_z)
    velocity_y = copy(velocity_z)
    S_λ_grid = Array{Float64, 4}(undef, (nλ, nz, nx, ny))u"kW*nm^-1*m^-2"
    α_grid = Array{Float64, 4}(undef, (nλ, nz, nx, ny))u"m^-1"
    populations_grid = Array{Float64, 4}(undef, (nz, nx, ny, 3))u"m^-3"

    tree = KDTree(ustrip(sites.positions))
    for k in 1:length(z)
        for i in 1:length(x)
            for j in 1:length(y)
                grid_point = [z[k], x[i], y[j]]
                idx, dist = nn(tree, ustrip(grid_point))
                temperature[k, i, j] = sites.temperature[idx]
                electron_density[k, i, j] = sites.electron_density[idx]
                hydrogen_populations[k, i, j] = sites.hydrogen_populations[idx]
                velocity_z[k, i, j] = sites.velocity_z[idx]
                velocity_x[k, i, j] = sites.velocity_x[idx]
                velocity_y[k, i, j] = sites.velocity_y[idx]
                S_λ_grid[:, k, i, j] = S_λ[:, idx]
                α_grid[:, k, i, j] = α_tot[:, idx]
                populations_grid[k, i, j, :] = populations[idx, :]
            end
        end
    end

    voronoi_atmos = Atmosphere(z, x, y, temperature, electron_density,
                           hydrogen_populations, velocity_z, velocity_x, velocity_y)

    return voronoi_atmos, S_λ_grid, α_grid, populations_grid
end

function _initialise(p_vec, atmos::Atmosphere)
    println("---Interpolating quantities to new grid---")
    n_sites = length(p_vec[1,:])

    temperature_new = Vector{Float64}(undef, n_sites)u"K"
    N_e_new = Vector{Float64}(undef, n_sites)u"m^-3"
    N_H_new = Vector{Float64}(undef, n_sites)u"m^-3"
    velocity_z_new = Vector{Float64}(undef, n_sites)u"m*s^-1"
    velocity_x_new = Vector{Float64}(undef, n_sites)u"m*s^-1"
    velocity_y_new = Vector{Float64}(undef, n_sites)u"m*s^-1"

    for k in 1:n_sites
        zk, xk, yk = p_vec[:, k]
        temperature_new[k] = trilinear(zk, xk, yk, atmos, atmos.temperature)
        N_e_new[k] = trilinear(zk, xk, yk, atmos, atmos.electron_density)
        N_H_new[k] = trilinear(zk, xk, yk, atmos, atmos.hydrogen_populations)
        velocity_z_new[k] = trilinear(zk, xk, yk, atmos, atmos.velocity_z)
        velocity_x_new[k] = trilinear(zk, xk, yk, atmos, atmos.velocity_x)
        velocity_y_new[k] = trilinear(zk, xk, yk, atmos, atmos.velocity_y)
    end
    return temperature_new, N_e_new, N_H_new, velocity_z_new, velocity_x_new, velocity_y_new
end

function line_of_sight_velocity(sites::VoronoiSites, k::Vector)
    v_los = Vector{Unitful.Velocity}(undef, sites.n)

    for ii in 1:sites.n
        velocity = [sites.velocity_z[ii],
                    sites.velocity_x[ii],
                    sites.velocity_y[ii]]
        v_los[ii] = dot(velocity, k)
    end
    return v_los::Vector{Unitful.Velocity}
end
