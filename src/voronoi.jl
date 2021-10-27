using NearestNeighbors

include("functions.jl")
include("voronoi_utils.jl")

global my_seed = 1001
Random.seed!(my_seed)

DATA = "../data/bifrost_qs006023_s525_quarter.hdf5"
atmos = Atmosphere(get_atmos(DATA; periodic=false, skip=3)...)

nz = length(atmos.z)
nx = length(atmos.x)
ny = length(atmos.y)

n_sites = nz*nx*ny

p_vec = rejection_sampling(n_sites, atmos)

# sort z array
p_vec = sort_array(p_vec; axis=1)

function plot_sites()
    pyplot()
    scatter(x_new/1e6,
            y_new/1e6,
            z_new/1e6,
            dpi=300)
    savefig("sites.png")
end

# plot_sites()

function column_count_mass()
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
                if x_bounds[1] < p_vec[2, l] < x_bounds[2] && y_bounds[1] < p_vec[3, l] < y_bounds[2]
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

 # column_count_mass()

#=
@time begin
    for k in 1:nz-1
        for i in 1:nx-1
            for j in 1:ny-1
                mass[(k-1)*nx*ny + (i-1)*ny + j] = mass_function(k, i, j, atmos)
                z_bounds = (atmos.z[k], atmos.z[k+1])
                x_bounds = (atmos.x[i], atmos.x[i+1])
                y_bounds = (atmos.y[j], atmos.y[j+1])
                sites_enclosed[(k-1)*nx*ny + (i-1)*ny + j] = find_sites_sorted(
                                            p_vec, z_bounds, x_bounds, y_bounds)
            end
        end
    end
end
=#

sites_file = "../data/sites.txt"
neighbours_file = "../data/neighbours.txt"

# write sites to file
write_arrays(ustrip(p_vec[1, :]),
             ustrip(p_vec[2, :]),
             ustrip(p_vec[3, :]),
             sites_file)

x_min = atmos.x[1]
x_max = atmos.x[end]
y_min = atmos.y[1]
y_max = atmos.y[end]
z_min = atmos.z[1]
z_max = atmos.z[end]
# export sites to voro++, and compute grid information
println("---Preprocessing grid---")
#=
run(`./voro.sh $sites_file $neighbours_file
               $(ustrip(x_min)) $(ustrip(x_max))
               $(ustrip(y_min)) $(ustrip(y_max))
               $(ustrip(z_min)) $(ustrip(z_max))`)
=#

temperature_new = Vector{Unitful.Temperature}(undef, n_sites)
N_e_new = Vector{NumberDensity}(undef, n_sites)
N_H_new = Vector{NumberDensity}(undef, n_sites)
for k in 1:n_sites
    temperature_new[k] = trilinear(p_vec[1, k], p_vec[2, k], p_vec[3, k], atmos, atmos.temperature)
    N_e_new[k] = trilinear(p_vec[1, k], p_vec[2, k], p_vec[3, k], atmos, atmos.electron_density)
    N_H_new[k] = trilinear(p_vec[1, k], p_vec[2, k], p_vec[3, k], atmos, atmos.hydrogen_populations)
end

# Create a tree
tree = KDTree(ustrip(p_vec))

# Voronoi grid
sites = VoronoiSites(p_vec[1,:], p_vec[2,:], p_vec[3,:], temperature_new, N_e_new, N_H_new)
cells = read_neighbours(neighbours_file, n_sites, sites)

raster = RasterDomain(collect(LinRange(z_min, z_max, 100)),
                      collect(LinRange(x_min, x_max, 100)),
                      collect(LinRange(y_min, y_max, 100)))

site = [raster.z[10], raster.x[25], raster.y[80]]

k = 30
nearest_neighbours = knn(tree, ustrip(raster_site), k)

i1 = inv_dist_itp(nearest_neighbours[1], nearest_neighbours[2]u"m", 1, sites)
i2 = inv_dist_itp(nearest_neighbours[1], nearest_neighbours[2]u"m", 2, sites)
i3 = inv_dist_itp(nearest_neighbours[1], nearest_neighbours[2]u"m", 3, sites)
i4 = inv_dist_itp(nearest_neighbours[1], nearest_neighbours[2]u"m", 4, sites)
i5 = inv_dist_itp(nearest_neighbours[1], nearest_neighbours[2]u"m", 5, sites)
i6 = inv_dist_itp(nearest_neighbours[1], nearest_neighbours[2]u"m", 6, sites)




































#
