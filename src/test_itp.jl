using Plots
using Random
using NearestNeighbors

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

function sample_function(x, y)
    return exp(-x^2 - y^2)
end

n_samples = 100

samples = rand(2, n_samples).*2 .-1
x_samples = samples[1,:]
y_samples = samples[2,:]

# Create a tree
tree = KDTree(samples)

k = Int(n_samples/4)

function_values = sample_function.(x_samples, y_samples)

raster_size = 250
x_range = LinRange(-1,1,raster_size)
y_range = LinRange(-1,1,raster_size)

gr()

function plot_the_itp(p)

    raster_values = Matrix{Float64}(undef, (raster_size, raster_size))

    for i in 1:raster_size
        for j in 1:raster_size
            raster_site = [x_range[i], y_range[j]]
            nearest_neighbours = knn(tree, raster_site, k)
            raster_values[i, j] = inv_dist_itp_test(nearest_neighbours[1],
                                                    nearest_neighbours[2],
                                                    p,
                                                    function_values)
        end
    end
    heatmap(x_range,
            y_range,
            raster_values,
            title="p=$p, $k neighbours",
            dpi=400)
    savefig("../img/inv_dist_$(p)_$(k)")
end

p = [1, 2, 3, 4, 5, 7, 10, 15, 20, 30, 40, 50]
for i in p
    plot_the_itp(i)
end

true_values = Matrix{Float64}(undef, (raster_size, raster_size))
for i in 1:raster_size
    for j in 1:raster_size
        true_values[i, j] = sample_function(x_range[i], y_range[j])
    end
end

heatmap(x_range,
        y_range,
        true_values,
        title="comparison",
        dpi=400)
savefig("../img/inv_dist_comparison")
