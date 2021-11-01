using Plots
using UnitfulRecipes
using NPZ

include("functions.jl")
include("voronoi_utils.jl")
include("irregular_ray_tracing.jl")

global my_seed = 1001
Random.seed!(my_seed)


function plot_sites(sites::VoronoiSites)
    pyplot()
    scatter(sites.x/1e6,
            sites.y/1e6,
            sites.z/1e6,
            dpi=300)
    savefig("sites.png")
end

function column_count_mass(atmos::Atmosphere, sites::VoronoiSites)
    # Not the actual mass, but column number density of hydrogen...
    column_mass = Array{Unitful.Mass, 2}(undef, (nx-1, ny-1))
    column_sites = Array{Int64, 2}(undef, (nx-1, ny-1))

    for i in 1:nx-1
        x_bounds = (atmos.x[i], atmos.x[i+1])
        for j in 1:ny-1
            print("column $((i-1)*(ny-1)+j)/$((nx-1)*(ny-1)) \r")
            cm = 0u"kg"
            y_bounds = (atmos.y[j], atmos.y[j+1])
            for k in 1:nz-1
                cm += mass_function(k, i, j, atmos)
            end
            hits = 0
            for l in 1:n_sites
                if x_bounds[1] < sites.x[l] < x_bounds[2] && y_bounds[1] < sites.y[l] < y_bounds[2]
                    hits += 1
                end
            end
            column_mass[i, j] = cm
            column_sites[i, j] = hits
        end
    end
    print("                                                                 \r")

    println("$(sum(column_sites)) and $n_sites should match")
    return column_mass, column_sites
end

function _initialise(p_vec, atmos::Atmosphere)
    n_sites = length(atmos.x)*length(atmos.z)*length(atmos.y)

    temperature_new = Vector{Unitful.Temperature}(undef, n_sites)
    N_e_new = Vector{NumberDensity}(undef, n_sites)
    N_H_new = Vector{NumberDensity}(undef, n_sites)

    for k in 1:n_sites
        temperature_new[k] = trilinear(p_vec[1, k], p_vec[2, k], p_vec[3, k], atmos, atmos.temperature)
        N_e_new[k] = trilinear(p_vec[1, k], p_vec[2, k], p_vec[3, k], atmos, atmos.electron_density)
        N_H_new[k] = trilinear(p_vec[1, k], p_vec[2, k], p_vec[3, k], atmos, atmos.hydrogen_populations)
    end
    return ustrip(p_vec[1,:]), ustrip(p_vec[2,:]), ustrip(p_vec[3,:]), ustrip(temperature_new), ustrip(N_e_new), ustrip(N_H_new)
end

function main()
    DATA = "../data/bifrost_qs006023_s525_quarter.hdf5"
    global atmos
    atmos = Atmosphere(get_atmos(DATA; periodic=false, skip=8)...)

    nx = length(atmos.x)
    nz = length(atmos.z)
    ny = length(atmos.y)

    n_sites = nz*nx*ny

    p_vec = rejection_sampling(n_sites, atmos, log10.(ustrip.(atmos.hydrogen_populations)))

    # sort z array
    p_vec = sort_array(p_vec; axis=1)


    sites_file = "../data/sites.txt"
    neighbours_file = "../data/neighbours.txt"

    # write sites to file
    write_arrays(ustrip(p_vec[1, :]),
                 ustrip(p_vec[2, :]),
                 ustrip(p_vec[3, :]),
                 sites_file)

    x_min = ustrip(atmos.x[1])
    x_max = ustrip(atmos.x[end])
    y_min = ustrip(atmos.y[1])
    y_max = ustrip(atmos.y[end])
    z_min = ustrip(atmos.z[1])
    z_max = ustrip(atmos.z[end])

    # export sites to voro++, and compute grid information
    println("---Preprocessing grid---")


    # compute neigbours
    run(`./voro.sh $sites_file $neighbours_file
                   $(x_min) $(x_max)
                   $(y_min) $(y_max)
                   $(z_min) $(z_max)`)

    # Create a tree
    tree = KDTree(ustrip(p_vec))

    # Voronoi grid
    sites = VoronoiSites(_initialise(p_vec, atmos)...,
                         z_min, z_max, x_min, x_max, y_min, y_max)

    # plot_sites(sites)
    # column_count_mass(atmos, sites)

    # Creates an array of VoronoiCell
    cells = read_cell(neighbours_file, n_sites, sites)

    #=
    site = [(z_max - z_min)/2,
            (x_max - x_min)/2,
            (y_max - y_min)/2]

    k = 30
    nearest_neighbours = knn(tree, site, k, true)
    =#

    #I_0 = 0
    #irregular_SC_up(sites, cells, tree, I_0)
    layer = zeros(Int, n_sites)

    lower_boundary = -5
    for cell in cells
        if any(i -> i==lower_boundary, cell.neighbours)
            layer[cell.ID] = 1
        end
    end

    lower_layer=1
    while true
        for cell in cells
            if layer[cell.ID] == 0
                for neighbor_ID in cell.neighbours
                    if neighbor_ID > 0 && layer[neighbor_ID] == lower_layer
                        layer[cell.ID] = lower_layer+1
                        break
                    end
                end
            end
        end

        if !any(i -> i==0, layer)
            break
        end

        lower_layer += 1
    end


    #= # This might be useful
    raster_size = 250
    x_range = LinRange(x_min,x_max,raster_size)
    z_range = LinRange(z_min,z_max,raster_size)
    y_range = LinRange(y_min,y_max,raster_size)

    layer_values = Array{Int, 3}(undef, (raster_size, raster_size, raster_size))

    for i in 1:raster_size
        for j in 1:raster_size
            for k in 1:raster_size
                raster_site = [z_range[i], x_range[j], y_range[k]]
                nearest_neighbour = nn(tree, raster_site)
                layer_values[i, j, k] = layer[nearest_neighbour[1]]
            end
        end
    end
    =#

    npzwrite("../python/p.npy", ustrip(p_vec))
    npzwrite("../python/layer.npy", layer)

end

main()
print("")

































#
