include("plot_utils.jl")

θ = 180.0
ϕ = 0.0

#title = "boost_regular_12"

fname = "../data/voronoi_ul7n12_2e6.h5"
title = split(fname, "/")[end]
title = split(title, ".")[1]

plot_convergence(fname, string(title))
write_convergence(fname, string(title))

atmos, line, S_λ, α_tot = plotter(read_irregular(fname)..., θ, ϕ, fname)
# atmos, line, S_λ, α_tot = plotter(read_quantities(fname, periodic=true)..., θ, ϕ, fname)
write_top_intensity(atmos, line, S_λ, α_tot, θ, ϕ, "$(title)_disk_centre_1")

# plot_top_line(atmos, line, S_λ, α_tot, θ, ϕ, title)

# for idλ in line.λidx[1]+1 : line.λidx[2]
#     plot_top_intensity(atmos, line, S_λ, α_tot, θ, ϕ, idλ, title*string(idλ))
# end


# plot_convergence("../data/voronoi_ul2n3_2.h5", "Irregular grid convergence")
# plot_convergence("../data/regular_ul2n3_zero_radiation_1.h5", "Regular grid convergence")
