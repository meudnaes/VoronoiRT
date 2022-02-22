using NearestNeighbors

include("functions.jl")

"""
    VoronoiSites

Structure containing quantities on the irregular grid, and information relevant
for ray-tracing.
"""
struct VoronoiSites
    positions::Matrix{typeof(1.0u"m")}
    neighbours::Matrix{Int}
    layers_up::Vector{Int}
    layers_down::Vector{Int}
    perm_up::Vector{Int}
    perm_down::Vector{Int}
    temperature::Vector{typeof(1.0u"K")}
    electron_density::Vector{typeof(1.0u"m^-3")}
    hydrogen_populations::Vector{typeof(1.0u"m^-3")}
    velocity_z::Vector{typeof(1.0u"m*s^-1")}
    velocity_x::Vector{typeof(1.0u"m*s^-1")}
    velocity_y::Vector{typeof(1.0u"m*s^-1")}
    z_min::typeof(1.0u"m")
    z_max::typeof(1.0u"m")
    x_min::typeof(1.0u"m")
    x_max::typeof(1.0u"m")
    y_min::typeof(1.0u"m")
    y_max::typeof(1.0u"m")
    n::Int
end

"""
    read_neighbours(fname::String, n_sites::Int, positions::Matrix{<:Unitful.Length})

Reads a file containing neighbouring cells for each grid point in the voronoi
tesselation.
"""
function read_cell(fname::String, n_sites::Int, positions::Matrix{<:Unitful.Length})
    println("---Reading neighbour information---")
    # Guess (overshoot), maybe do exact later?
    max_guess = 70
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

    NeighbourMatrix = NeighbourMatrix[:,1:max_neighbours+1]
    layers_up = _sort_by_layer_up(NeighbourMatrix, n_sites)
    perm_up = sortperm(layers_up)
    layers_up = layers_up[perm_up]
    layers_up = reduce_layers(layers_up)

    layers_down = _sort_by_layer_down(NeighbourMatrix, n_sites)
    perm_down = sortperm(layers_down)
    layers_down = layers_down[perm_down]
    layers_down = reduce_layers(layers_down)

    return positions, NeighbourMatrix, layers_up, layers_down, perm_up, perm_down
end

"""
    _sort_by_layer_up(neighbours::Matrix{Int}, n_sites::Int)

Sort Voronoi sites by layer, from all cells neighbouring the bottom boundary and
upwards. Second layer consists of sites neighbouring first layer, and so on.
"""
function _sort_by_layer_up(neighbours::Matrix{Int}, n_sites::Int)

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

"""
    _sort_by_layer_down(neighbours::Matrix{Int}, n_sites::Int)

Sort Voronoi sites by layer, from all cells neighbouring the top boundary and
downwards. Second layer consists of sites neighbouring first layer, and so on.
"""
function _sort_by_layer_down(neighbours::Matrix{Int}, n_sites::Int)
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

"""
    reduce_layers(layers::Vector{Int})

Reduce vector containg layers. Returns a compressed vector where the element of
each index tells how many sites are in the layer equalling the index+1.
"""
function reduce_layers(layers::Vector{Int})
    reduced_layers = Vector{Int}(undef, maximum(layers)+1)

    reduced_layers[1] = 1

    layer = 2
    for i in eachindex(layers)
        if layers[i] == layer
            reduced_layers[layer] = i
            layer += 1
        end
    end

    reduced_layers[end] = length(layers)

    return reduced_layers
end

"""
    smallest_angle(position::Vector{<:Unitful.Length},
                        neighbours::Vector{Int},
                        k::Vector{Float64},
                        sites::VoronoiSites)

Calculates the dot product between the Delaunay lines connecting a site and
every neighbour, and the ray travelling in direction k. Keep in mind that k
is not toward the ray, but with the direction of the ray. Returns the two
Delaunay lines with the largest dot products, and the indices of their sites.
"""
function smallest_angle(position::Vector{<:Unitful.Length},
                        neighbours::Vector{Int},
                        k::Vector{Float64},
                        sites::VoronoiSites)

    dots = Vector{Float64}(undef, 2)
    fill!(dots, -1)

    indices = Vector{Int}(undef, 2)

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

            direction = p_n .- position
            norm_dir = direction/(norm(direction))
            # Two normalized direction vectors, denominator is 1
            # Save the dot product, don't bother calculating the angle
            dot_product = dot(k, norm_dir)

            if dot_product > dots[2]
                if dot_product > dots[1]
                    dots[1] = dot_product
                    indices[1] = neighbour
                else
                    dots[2] = dot_product
                    indices[2] = neighbour
                end
            end
        end
    end

    if dots[2] <= 0
        dots[2] = 0
        indices[2] = indices[1]
    end

    return dots, indices
end

"""
    Voronoi_to_Raster(sites::VoronoiSites,
                           atmos::Atmosphere,
                           S_λ::Array{<:UnitsIntensity_λ},
                           α_tot::Array{<:PerLength},
                           populations::Array{<:NumberDensity},
                           r_factor)
Intepolate quantities from an irregular grid back to a regular grid.
"""
function Voronoi_to_Raster(sites::VoronoiSites,
                           atmos::Atmosphere,
                           S_λ::Matrix{<:UnitsIntensity_λ},
                           α_tot::Matrix{<:PerLength},
                           populations::Matrix{<:NumberDensity},
                           r_factor::Any;
                           periodic=false)

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

    if periodic
        voronoi_atmos = Atmosphere(z,
                                   periodic_borders(x),
                                   periodic_borders(y),
                                   periodic_borders(temperature),
                                   periodic_borders(electron_density),
                                   periodic_borders(hydrogen_populations),
                                   periodic_borders(velocity_z),
                                   periodic_borders(velocity_x),
                                   periodic_borders(velocity_y))

        S_λ_grid = periodic_borders(S_λ_grid)
        α_grid = periodic_borders(α_grid)
        populations_grid = periodic_pops(populations_grid)
    else
        voronoi_atmos = Atmosphere(z, x, y, temperature,
                                   electron_density, hydrogen_populations,
                                   velocity_z, velocity_x, velocity_y)
    end

    return voronoi_atmos, S_λ_grid, α_grid, populations_grid
end

function Voronoi_to_Raster(sites::VoronoiSites,
                           atmos::Atmosphere,
                           S_λ::Array{<:UnitsIntensity_λ},
                           α_tot::Array{<:PerLength},
                           r_factor;
                           periodic=false)

    z = collect(LinRange(sites.z_min, sites.z_max,
                         floor(Int, r_factor*length(atmos.z))))
    x = collect(LinRange(sites.x_min, sites.x_max,
                         floor(Int, r_factor*length(atmos.x))))
    y = collect(LinRange(sites.y_min, sites.y_max,
                         floor(Int, r_factor*length(atmos.y))))

    nz = length(z)
    nx = length(x)
    ny = length(y)

    temperature = Array{Float64, 3}(undef, (nz, ny, nx))u"K"
    electron_density = Array{Float64, 3}(undef, (nz, ny, nx))u"m^-3"
    hydrogen_populations = Array{Float64, 3}(undef, (nz, ny, nx))u"m^-3"
    velocity_z = Array{Float64, 3}(undef, (nz, ny, nx))u"m*s^-1"
    velocity_x = copy(velocity_z)
    velocity_y = copy(velocity_z)
    S_λ_grid = Array{Float64, 3}(undef, (nz, nx, ny))u"kW*nm^-1*m^-2"
    α_grid = Array{Float64, 3}(undef, (nz, nx, ny))u"m^-1"

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
                S_λ_grid[k, i, j] = S_λ[idx]
                α_grid[k, i, j] = α_tot[idx]
            end
        end
    end

    if periodic
        voronoi_atmos = Atmosphere(z,
                                   periodic_borders(x),
                                   periodic_borders(y),
                                   periodic_borders(temperature),
                                   periodic_borders(electron_density),
                                   periodic_borders(hydrogen_populations),
                                   periodic_borders(velocity_z),
                                   periodic_borders(velocity_x),
                                   periodic_borders(velocity_y))

        S_λ_grid = periodic_borders(S_λ_grid)
        α_grid = periodic_borders(α_grid)
    else
        voronoi_atmos = Atmosphere(z, x, y, temperature,
                                   electron_density, hydrogen_populations,
                                   velocity_z, velocity_x, velocity_y)
    end

    return voronoi_atmos, S_λ_grid, α_grid
end

function Voronoi_to_Raster(sites::VoronoiSites,
                           atmos::Atmosphere,
                           r_factor;
                           periodic=false)

    # z = collect(LinRange(sites.z_min, sites.z_max,
                         # floor(Int, r_factor*length(atmos.z))))
    # x = collect(LinRange(sites.x_min, sites.x_max,
                         # floor(Int, r_factor*length(atmos.x))))
    # y = collect(LinRange(sites.y_min, sites.y_max,
                         # floor(Int, r_factor*length(atmos.y))))

    z = atmos.z
    x = atmos.x
    y = atmos.y

    nz = length(z)
    nx = length(x)
    ny = length(y)

    temperature = Array{Float64, 3}(undef, (nz, ny, nx))u"K"
    electron_density = Array{Float64, 3}(undef, (nz, ny, nx))u"m^-3"
    hydrogen_populations = Array{Float64, 3}(undef, (nz, ny, nx))u"m^-3"
    velocity_z = Array{Float64, 3}(undef, (nz, ny, nx))u"m*s^-1"
    velocity_x = copy(velocity_z)
    velocity_y = copy(velocity_z)

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
            end
        end
    end

    if periodic
        voronoi_atmos = Atmosphere(z,
                                   periodic_borders(x),
                                   periodic_borders(y),
                                   periodic_borders(temperature),
                                   periodic_borders(electron_density),
                                   periodic_borders(hydrogen_populations),
                                   periodic_borders(velocity_z),
                                   periodic_borders(velocity_x),
                                   periodic_borders(velocity_y))
    else
        voronoi_atmos = Atmosphere(z, x, y, temperature,
                                   electron_density, hydrogen_populations,
                                   velocity_z, velocity_x, velocity_y)
    end

    return voronoi_atmos
end


function _initialise(p_vec::Matrix{<:Unitful.Length}, atmos::Atmosphere)
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

function _initialiseII(p_vec::Matrix{<:Unitful.Length}, atmos::Atmosphere)
    println("---Interpolating quantities to new grid---")
    n_sites = length(p_vec[1,:])

    temperature_new = Vector{Float64}(undef, n_sites)u"K"
    N_e_new = Vector{Float64}(undef, n_sites)u"m^-3"
    N_H_new = Vector{Float64}(undef, n_sites)u"m^-3"
    velocity_z_new = Vector{Float64}(undef, n_sites)u"m*s^-1"
    velocity_x_new = Vector{Float64}(undef, n_sites)u"m*s^-1"
    velocity_y_new = Vector{Float64}(undef, n_sites)u"m*s^-1"

    for k in 1:n_sites
        p_k = p_vec[:, k]

        idz = searchsortedfirst(atmos.z, p_k[1]) - 1
        idx = searchsortedfirst(atmos.x, p_k[2]) - 1
        idy = searchsortedfirst(atmos.y, p_k[3]) - 1

        positions = [atmos.z[idz]   atmos.x[idx]    atmos.y[idy]
                     atmos.z[idz]   atmos.x[idx]    atmos.y[idy+1]
                     atmos.z[idz]   atmos.x[idx+1]  atmos.y[idy]
                     atmos.z[idz]   atmos.x[idx+1]  atmos.y[idy+1]
                     atmos.z[idz+1] atmos.x[idx]    atmos.y[idy]
                     atmos.z[idz+1] atmos.x[idx]    atmos.y[idy+1]
                     atmos.z[idz+1] atmos.x[idx+1]  atmos.y[idy]
                     atmos.z[idz+1] atmos.x[idx+1]  atmos.y[idy+1]]
        positions = transpose(positions)

        indices = [idz     idx     idy
                   idz     idx     idy+1
                   idz     idx+1   idy
                   idz     idx+1   idy+1
                   idz+1   idx     idy
                   idz+1   idx     idy+1
                   idz+1   idx+1   idy
                   idz+1   idx+1   idy+1]
        indices = transpose(indices)

        distances = Vector{Float64}(undef, 8)u"m"
        for i in 1:8
            distances[i] = euclidean(positions[:, i], p_k)
        end

        index = indices[:, argmin(distances)]

        temperature_new[k] = atmos.temperature[index...]
        N_e_new[k] = atmos.electron_density[index...]
        N_H_new[k] = atmos.electron_density[index...]
        velocity_z_new[k] = atmos.velocity_z[index...]
        velocity_x_new[k] = atmos.velocity_x[index...]
        velocity_y_new[k] = atmos.velocity_y[index...]
    end
    return temperature_new, N_e_new, N_H_new, velocity_z_new, velocity_x_new, velocity_y_new
end
