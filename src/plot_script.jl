θ = 180.0
ϕ = 0.0

#title = "boost_regular_12"

fname = "../data/avg_ext_3e6.h5"
title = split(fname, "/")[end]
title = split(title, ".")[1]
title = string(title)

plot_convergence(fname, title)
write_convergence(fname, title)

atmos, line, S_λ, α_tot = plotter(read_irregular(fname)..., θ, ϕ, fname)
# atmos, line, S_λ, α_tot = plotter(read_quantities(fname, periodic=true)..., θ, ϕ, fname)
write_top_intensity(atmos, line, S_λ, α_tot, θ, ϕ, "$(title)_disk_centre_1")

# atmos, S_λ, populations = read_irregular(fname)
# atmos, S_λ, populations = read_quantities(fname, periodic=false)
# write_source_function(S_λ, atmos, title)

# write_to_original(fname, "../data/half_res_ul7n12.h5", "$(title)_2")
