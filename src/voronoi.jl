include("functions.jl")

DATA = "../data/bifrost_qs006023_s525_quarter.hdf5"
atmos = Atmosphere(get_atmos(DATA; periodic=false, skip=2)...)

nz = length(atmos.z)
nx = length(atmos.x)
ny = length(atmos.y)

n_sites = nz*nx*ny

p_vec, N_H_new = rejection_sampling(n_sites, atmos)

# sort z array
P = sortperm(p_vec[1, :])
p_vec = p_vec[:, P]

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

# Print boundaries for voro++ file
open("boundaries.txt", "w") do io
    println(io, "x_min = $(atmos.x[1])")
    println(io, "x_max = $(atmos.x[end])")
    println(io, "y_min = $(atmos.y[1])")
    println(io, "y_max = $(atmos.y[end])")
    println(io, "z_min = $(atmos.z[1])")
    println(io, "z_max = $(atmos.z[end])")
end

write_arrays(p_vec[1, :], p_vec[2, :], p_vec[3, :], "seeds.txt")

# initialise the system in LTE
temperature_new = Vector{Unitful.Temperature}(undef, n_sites)
N_e_new = Vector{NumberDensity}(undef, n_sites)
for k in 1:n_sites
    temperature_new[k] = trilinear(p_vec[1, k], p_vec[2, k], p_vec[3, k], atmos, atmos.temperature)
    N_e_new[k] = trilinear(p_vec[1, k], p_vec[2, k], p_vec[3, k], atmos, atmos.electron_density)
    # temperature_new[k] = trilinear(z_new[k], x_new[k], y_new[k], atmos, atmos.temperature)
    # N_e_new[k] = trilinear(z_new[k], x_new[k], y_new[k], atmos, atmos.electron_density)
end
