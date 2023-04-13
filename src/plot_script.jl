using .VoronoiRT

θ = 180.0
ϕ = 0.0

# fname = "../data/half_res_ul7n12.h5"
fname = "../data/invNH_invT_3e6.h5"
title = split(fname, "/")[end]
title = split(title, ".")[1]
title = string(title)

# plot_convergence(fname, title)
# VoronoiRT.write_convergence(fname, title)

cont_waves = 0

atmos, line, S_λ, α_tot = VoronoiRT.plotter(VoronoiRT.read_irregular(fname; ncont=cont_waves)..., θ, ϕ, fname)
# atmos, line, S_λ, α_tot = VoronoiRT.plotter(VoronoiRT.read_quantities(fname; periodic=true, ncont=cont_waves)..., θ, ϕ, fname)

VoronoiRT.write_top_intensity(atmos, line, S_λ, α_tot, θ, ϕ, "$(title)_disk_centre_inv_dist_c")

# atmos, S_λ, populations = read_irregular(fname)
# atmos, S_λ, populations = read_quantities(fname, periodic=false)
# write_source_function(S_λ, atmos, title)

# VoronoiRT.write_to_original(fname, "../data/half_res_ul7n12.h5", "$(title)_2")
# VoronoiRT.write_tau_unity(fname)
