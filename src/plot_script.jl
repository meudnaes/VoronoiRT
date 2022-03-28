include("plot_utils.jl")

θ = 180.0
ϕ = 0.0
title = "boost_voronoi_12_extinction"

fname = "../data/voronoi_ul7n12_extinction.h5"

plot_convergence(fname, "Irregular Extinction 12 rays")


# atmos, line, S_λ, α_tot = plotter(read_irregular(fname)..., θ, ϕ, title)
# atmos, line, S_λ, α_tot = plotter(read_quantities(fname, periodic=true)..., θ, ϕ, title)

# plot_top_line(atmos, line, S_λ, α_tot, θ, ϕ, title)

# for idλ in line.λidx[1]+1 : line.λidx[2]
#     plot_top_intensity(atmos, line, S_λ, α_tot, θ, ϕ, idλ, title*string(idλ))
# end


# plot_convergence("../data/voronoi_ul2n3_2.h5", "Irregular grid convergence")
# plot_convergence("../data/regular_ul2n3_zero_radiation_1.h5", "Regular grid convergence")
