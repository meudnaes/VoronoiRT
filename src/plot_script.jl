include("plot_utils.jl")

plot_convergence("../data/voronoi_ul2n3_2.h5", "Irregular grid convergence")
# plot_convergence("../data/regular_ul2n3_C.h5", "Boost convergence")
# plot_convergence("../data/regular_ul2n3_zero_radiation_1.h5", "Regular grid convergence")

# plotter(read_quantities("../data/regular_ul2n3_zero_radiation_1.h5", periodic=true)..., 180.0, 0.0, "Regular-Line")
# plotter(read_quantities("../data/regular_ul2n3_C.h5", periodic=true)..., 180.0, 0.0, "Boost-Line")
plotter(read_irregular("../data/voronoi_ul2n3_2.h5")..., 180.0, 0.0, "Voronoi-Line")
