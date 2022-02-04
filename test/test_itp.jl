using Plots
using Random
using NearestNeighbors

global my_seed = 1001
Random.seed!(my_seed)

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

function inv_dist_itp(idxs::AbstractVector,
                      dists::AbstractVector,
                      p::AbstractFloat,
                      values::AbstractVector)
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

function sample_function(x, y)
    return exp(-x^2 - y^2)
end

n_samples = 100

samples = rand(2, n_samples).*2 .-1
x_samples = samples[1,:]
y_samples = samples[2,:]

# Create a tree
tree = KDTree(samples)

k = 15

function_values = sample_function.(x_samples, y_samples)

raster_size = 500
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
            dpi=350,
            aspect_ratio=:equal,
            xlim=[-1,1])
    savefig("../img/InverseDistance/inv_dist_$(Int(p*100))_$(k)")
end

p = [3]
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
        title="True Function",
        dpi=350,
        aspect_ratio=:equal,
        xlim=[-1,1])
savefig("../img/InverseDistance/inv_dist_comparison")

scatter(x_samples,
        y_samples,
        title="samples",
        dpi=350,
        aspect_ratio=:equal,
        xlim=[-1,1],
        xlabel=raw"$x$",
        ylabel=raw"$y$")
savefig("../img/InverseDistance/sample_points")
