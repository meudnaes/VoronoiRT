using Plots
include("voronoi_utils.jl")
include("lambda_continuum.jl")
include("irregular_ray_tracing.jl")

global my_seed = 1998
Random.seed!(my_seed)

pyplot()

function plot_sites(sites::VoronoiSites)
    z = ustrip.(sites.positions[1,:])
    x = ustrip.(sites.positions[2,:])
    y = ustrip.(sites.positions[3,:])
    scatter(x/1e6,
            y/1e6,
            z/1e6,
            title="Sites",
            dpi=300)
    savefig("../img/sites.png")
end

function compare(DATA, quadrature)
    maxiter = 50
    ϵ = 1e-4

    θ = 170.0
    ϕ = 30.0

    n_skip = 1

    function regular()

        atmos = Atmosphere(get_atmos(DATA; periodic=true, skip=n_skip)...)

        @time J_mean, S_λ, α_tot = Λ_regular(ϵ, maxiter, atmos, quadrature)

        k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]
        I_top = short_characteristics_up(k, S_λ, α_tot, atmos, I_0=S_λ[1,:,:])

        I_top = ustrip(uconvert.(u"kW*nm^-1*m^-2", I_top[end, 2:end-1, 2:end-1]))

        # global min_lim, max_lim
        # min_lim = minimum(I_top)
        # max_lim = maximum(I_top)

        heatmap(ustrip.(atmos.x[2:end-1]),
                ustrip.(atmos.y[2:end-1]),
                transpose.(I_top),
                xaxis="x",
                yaxis="y",
                dpi=300,
                rightmargin=10Plots.mm,
                title="Regular Grid",
                aspect_ratio=:equal)
                # clim=(min_lim,max_lim))

        savefig("../img/compare_continuum/regular_tilted")

        return 0
    end


    function voronoi()

        atmos = Atmosphere(get_atmos(DATA; periodic=false, skip=n_skip)...)

        nx = length(atmos.x)
        nz = length(atmos.z)
        ny = length(atmos.y)

        x_min = ustrip(atmos.x[1])
        x_max = ustrip(atmos.x[end])
        y_min = ustrip(atmos.y[1])
        y_max = ustrip(atmos.y[end])
        z_min = ustrip(atmos.z[1])
        z_max = ustrip(atmos.z[end])

        n_sites = floor(Int, nz*nx*ny)
        positions = rand(3, n_sites)

        positions[1, :] = positions[1, :].*(z_max - z_min) .+ z_min
        positions[2, :] = positions[2, :].*(x_max - x_min) .+ x_min
        positions[3, :] = positions[3, :].*(y_max - y_min) .+ y_min

        positions = positions*1u"m"

        # rejection_sampling(n_sites, atmos, log10.(ustrip.(atmos.hydrogen_populations)))
        # sample_from_extinction(atmos, 500.0u"nm", n_sites)

        sites_file = "../data/sites_continuum.txt"
        neighbours_file = "../data/neighbours_continuum.txt"
        # write sites to file
        write_arrays(ustrip.(positions[2, :]),
                     ustrip.(positions[3, :]),
                     ustrip.(positions[1, :]),
                     sites_file)

        # export sites to voro++, and compute grid information
        println("---Preprocessing grid---")

        # compute neigbours
        run(`./voro.sh $sites_file $neighbours_file
                       $(x_min) $(x_max)
                       $(y_min) $(y_max)
                       $(z_min) $(z_max)`)

        # Voronoi grid
        sites = VoronoiSites(read_cell(neighbours_file, n_sites, positions)...,
                             _initialise(positions, atmos)...,
                             z_min*1u"m", z_max*1u"m",
                             x_min*1u"m", x_max*1u"m",
                             y_min*1u"m", y_max*1u"m",
                             n_sites)

        # plot_sites(sites)

        @time J_mean, S_λ, α_tot = Λ_voronoi(ϵ, maxiter, sites, quadrature)

        atmos_from_voronoi, S_λ_grid, α_grid = Voronoi_to_Raster(sites, atmos,
                                                                 S_λ, α_tot, 1;
                                                                 periodic=true)

        k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]
        I_top = short_characteristics_up(k, S_λ_grid, α_grid,
                                        atmos_from_voronoi, I_0=S_λ_grid[1,:,:])

        I_top = ustrip(uconvert.(u"kW*nm^-1*m^-2", I_top[50, 2:end-1, 2:end-1]))

        heatmap(ustrip.(atmos_from_voronoi.x[2:end-1]),
                ustrip.(atmos_from_voronoi.y[2:end-1]),
                transpose.(I_top),
                xaxis="x",
                yaxis="y",
                dpi=300,
                rightmargin=10Plots.mm,
                title="Irregular Grid",
                aspect_ratio=:equal)
                # clim=(min_lim,max_lim))

        savefig("../img/compare_continuum/irregular_50_n1")

        return 0
    end

    regular();
    # voronoi();

end

function LTE_ray(DATA)

    atmos = Atmosphere(get_atmos(DATA; periodic=false, skip=1)...)

    # choose a wavelength
    λ = 500u"nm"  # nm

    # Lte populations
    LTE_pops = LTE_ionisation(atmos)

    # Find continuum extinction (only with Thomson and Rayleigh)
    α_cont = α_continuum.(λ, atmos.temperature*1.0, atmos.electron_density*1.0,
                          LTE_pops[:,:,:,1]*1.0, LTE_pops[:,:,:,3]*1.0)

    # Planck function
    S_λ = blackbody_λ.(λ, atmos.temperature)

    θ = 170.0
    ϕ = 0.1

    k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]
    intensity = short_characteristics_up(k, S_λ, α_cont, atmos, I_0=S_λ[1,:,:])

    #=
    intensity = Array{Float64, 3}(undef, size(α_cont))u"kW*m^-2*nm^-1"
    for idx in eachindex(atmos.x)
        for idy in eachindex(atmos.y)
            intensity[:, idx, idy] = Transparency.piecewise_1D_linear(atmos.z,
                                              α_cont[:, idx, idy],
                                              B_λ[:, idx, idy];
                                              to_end=true,
                                              initial_condition=:source)
        end
    end
    =#

    I_top = ustrip(uconvert.(u"kW*nm^-1*m^-2", intensity[end, :, :]))
    μ = abs(k[1])

    heatmap(ustrip.(atmos.x),
            ustrip.(atmos.y),
            transpose.(I_top),
            xaxis="x",
            yaxis="y",
            dpi=300,
            rightmargin=10Plots.mm,
            title="Continuum at 500 nm",
            aspect_ratio=:equal)

    savefig("../img/compare_continuum/cont500_$(floor(Int, 100*μ))")

    return 0
end

function LTE_voronoi(DATA)

    atmos = Atmosphere(get_atmos(DATA; periodic=false, skip=1)...)

    nx = length(atmos.x)
    nz = length(atmos.z)
    ny = length(atmos.y)

    x_min = ustrip(atmos.x[1])
    x_max = ustrip(atmos.x[end])
    y_min = ustrip(atmos.y[1])
    y_max = ustrip(atmos.y[end])
    z_min = ustrip(atmos.z[1])
    z_max = ustrip(atmos.z[end])

    n_sites = floor(Int, nz*nx*ny)
    positions = rand(3, n_sites)

    positions[1, :] = positions[1, :].*(z_max - z_min) .+ z_min
    positions[2, :] = positions[2, :].*(x_max - x_min) .+ x_min
    positions[3, :] = positions[3, :].*(y_max - y_min) .+ y_min

    positions = positions*1u"m"

    # rejection_sampling(n_sites, atmos, log10.(ustrip.(atmos.hydrogen_populations)))
    # sample_from_extinction(atmos, 500.0u"nm", n_sites)

    sites_file = "../data/sites_continuum.txt"
    neighbours_file = "../data/neighbours_continuum.txt"
    # write sites to file
    write_arrays(ustrip.(positions[2, :]),
                 ustrip.(positions[3, :]),
                 ustrip.(positions[1, :]),
                 sites_file)

    # export sites to voro++, and compute grid information
    println("---Preprocessing grid---")

    # compute neigbours
    run(`./voro.sh $sites_file $neighbours_file
                   $(x_min-0.1) $(x_max+0.1)
                   $(y_min-0.1) $(y_max+0.1)
                   $(z_min-0.1) $(z_max+0.1)`)

    # Voronoi grid
    sites = VoronoiSites(read_cell(neighbours_file, n_sites, positions)...,
                         _initialise(positions, atmos)...,
                         z_min*1u"m", z_max*1u"m",
                         x_min*1u"m", x_max*1u"m",
                         y_min*1u"m", y_max*1u"m",
                         n_sites)

    # Lte populations
    LTE_pops = LTE_ionisation(sites)
    λ = 500u"nm"
    # Find continuum extinction (only with Thomson and Rayleigh)
    α_cont = α_continuum.(λ, sites.temperature*1.0, sites.electron_density*1.0,
                          LTE_pops[:,1]*1.0, LTE_pops[:,3]*1.0)

    S_λ = blackbody_λ.(λ, sites.temperature)

    k = [-1.0, 0.0, 0.0]
    bottom_layer = sites.layers_up[2] - 1
    bottom_layer_idx = sites.perm_up[1:bottom_layer]
    println("---Ray tracing---")
    I_0 = blackbody_λ.(500u"nm", sites.temperature[bottom_layer_idx])
    intensity = Delaunay_upII(k, S_λ, α_cont, sites, I_0, 3)

    x = collect(LinRange(sites.x_min, sites.x_max, 10*nx))
    y = collect(LinRange(sites.y_min, sites.y_max, 10*ny))
    top_z = atmos.z[50]

    tree = KDTree(ustrip(sites.positions))

    I_top = Matrix{Float64}(undef, (length(x), length(y)))u"kW*m^-2*nm^-1"
    for i in 1:length(x)
        for j in 1:length(y)
            position = [top_z, x[i], y[j]]
            idx, dist = nn(tree, ustrip(position))
            I_top[i, j] = intensity[idx]
        end
    end

    heatmap(ustrip.(x),
            ustrip.(y),
            ustrip.(transpose.(I_top)),
            xaxis="x",
            yaxis="y",
            dpi=300,
            rightmargin=10Plots.mm,
            title="Top Intensity, Irregular Grid",
            aspect_ratio=:equal)

    savefig("../img/compare_continuum/irregular_top_I")

    return 0
end

function test_interpolation(DATA, quadrature)
    atmos = Atmosphere(get_atmos(DATA; periodic=false, skip=1)...)

    nx = length(atmos.x)
    nz = length(atmos.z)
    ny = length(atmos.y)

    x_min = ustrip(atmos.x[1])
    x_max = ustrip(atmos.x[end])
    y_min = ustrip(atmos.y[1])
    y_max = ustrip(atmos.y[end])
    z_min = ustrip(atmos.z[1])
    z_max = ustrip(atmos.z[end])

    n_sites = floor(Int, nz*nx*ny)
    positions = rand(3, n_sites)

    positions[1, :] = positions[1, :].*(z_max - z_min) .+ z_min
    positions[2, :] = positions[2, :].*(x_max - x_min) .+ x_min
    positions[3, :] = positions[3, :].*(y_max - y_min) .+ y_min

    positions = positions*1u"m"

    # rejection_sampling(n_sites, atmos, log10.(ustrip.(atmos.hydrogen_populations)))
    # sample_from_extinction(atmos, 500.0u"nm", n_sites)

    sites_file = "../data/sites_continuum.txt"
    neighbours_file = "../data/neighbours_continuum.txt"
    # write sites to file
    write_arrays(ustrip.(positions[2, :]),
                 ustrip.(positions[3, :]),
                 ustrip.(positions[1, :]),
                 sites_file)

    # export sites to voro++, and compute grid information
    println("---Preprocessing grid---")

    # compute neigbours
    run(`./voro.sh $sites_file $neighbours_file
                   $(x_min) $(x_max)
                   $(y_min) $(y_max)
                   $(z_min) $(z_max)`)

    # Voronoi grid
    sites = VoronoiSites(read_cell(neighbours_file, n_sites, positions)...,
                         _initialise(positions, atmos)...,
                         z_min*1u"m", z_max*1u"m",
                         x_min*1u"m", x_max*1u"m",
                         y_min*1u"m", y_max*1u"m",
                         n_sites)

    atmos_from_voronoi = Voronoi_to_Raster(sites, atmos, 1; periodic=false)
    @time J_mean, S_λ, α_tot = Λ_regular(1e-4, 50, atmos_from_voronoi, quadrature)

    k = [-1.0, 0.0, 0.0]
    I_top = short_characteristics_up(k, S_λ, α_tot, atmos_from_voronoi, I_0=S_λ[1,:,:])

    I_top = ustrip(uconvert.(u"kW*nm^-1*m^-2", I_top[50, 2:end-1, 2:end-1]))

    # global min_lim, max_lim
    # min_lim = minimum(I_top)
    # max_lim = maximum(I_top)

    heatmap(ustrip.(atmos.x[2:end-1]),
            ustrip.(atmos.y[2:end-1]),
            transpose.(I_top),
            xaxis="x",
            yaxis="y",
            dpi=300,
            rightmargin=10Plots.mm,
            title="Regular Grid",
            aspect_ratio=:equal)
            # clim=(min_lim,max_lim))

    savefig("../img/compare_continuum/test_itp")

    return 0
end

function test_with_regular(DATA, quadrature)
    atmos = Atmosphere(get_atmos(DATA; periodic=false, skip=1)...)

    nx = length(atmos.x)
    nz = length(atmos.z)
    ny = length(atmos.y)

    x_min = ustrip(atmos.x[1])
    x_max = ustrip(atmos.x[end])
    y_min = ustrip(atmos.y[1])
    y_max = ustrip(atmos.y[end])
    z_min = ustrip(atmos.z[1])
    z_max = ustrip(atmos.z[end])

    n_sites = floor(Int, nz*nx*ny)
    positions = Matrix{Float64}(undef, (3, n_sites))u"m"
    temperature = Vector{Float64}(undef, n_sites)u"K"
    hydrogen_populations = Vector{Float64}(undef, n_sites)u"m^-3"
    electron_density = Vector{Float64}(undef, n_sites)u"m^-3"
    velocity_z = Vector{Float64}(undef, n_sites)u"m/s"
    velocity_x = Vector{Float64}(undef, n_sites)u"m/s"
    velocity_y = Vector{Float64}(undef, n_sites)u"m/s"


    for k in eachindex(atmos.z)
        for i in eachindex(atmos.x)
            for j in eachindex(atmos.y)
                p = [atmos.z[k], atmos.x[i], atmos.y[j]]
                m = k + (i-1)*nz + (j-1)*nz*nx
                positions[:,m] = p
                temperature[m] = atmos.temperature[k, i, j]
                hydrogen_populations[m] = atmos.hydrogen_populations[k, i, j]
                electron_density[m] = atmos.electron_density[k, i, j]
                velocity_z[m] = atmos.velocity_z[k, i, j]
                velocity_x[m] = atmos.velocity_x[k, i, j]
                velocity_y[m] = atmos.velocity_y[k, i, j]
            end
        end
    end

    sites_file = "../data/sites_continuum.txt"
    neighbours_file = "../data/neighbours_continuum.txt"
    # write sites to file
    write_arrays(ustrip.(positions[2, :]),
                 ustrip.(positions[3, :]),
                 ustrip.(positions[1, :]),
                 sites_file)

    # export sites to voro++, and compute grid information
    println("---Preprocessing grid---")

    # compute neigbours
    run(`./voro.sh $sites_file $neighbours_file
                   $(x_min) $(x_max)
                   $(y_min) $(y_max)
                   $(z_min) $(z_max)`)

    # Voronoi grid
    sites = VoronoiSites(read_cell(neighbours_file, n_sites, positions)...,
                         temperature, electron_density, hydrogen_populations,
                         velocity_z, velocity_x, velocity_y,
                         z_min*1u"m", z_max*1u"m",
                         x_min*1u"m", x_max*1u"m",
                         y_min*1u"m", y_max*1u"m",
                         n_sites)

    # Lte populations
    LTE_pops = LTE_ionisation(sites)
    λ = 500u"nm"
    # Find continuum extinction (only with Thomson and Rayleigh)
    α_cont = α_continuum.(λ, sites.temperature*1.0, sites.electron_density*1.0,
                       LTE_pops[:,1]*1.0, LTE_pops[:,3]*1.0)
    S_λ = blackbody_λ.(λ, sites.temperature)

    θ = 170.0
    ϕ = 30.0
    k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]
    bottom_layer = sites.layers_up[2] - 1
    bottom_layer_idx = sites.perm_up[1:bottom_layer]
    I_0 = blackbody_λ.(500u"nm", sites.temperature[bottom_layer_idx])
    intensity = Delaunay_upII(k, S_λ, α_cont, sites, I_0, 3)

    x = collect(LinRange(sites.x_min, sites.x_max, 5*nx))
    y = collect(LinRange(sites.y_min, sites.y_max, 5*ny))
    top_z = sites.z_max

    tree = KDTree(ustrip(sites.positions))

    I_top = Matrix{Float64}(undef, (length(x), length(y)))u"kW*m^-2*nm^-1"
    for i in 1:length(x)
        for j in 1:length(y)
            position = [top_z, x[i], y[j]]
            idx, dist = nn(tree, ustrip(position))
            I_top[i, j] = intensity[idx]
        end
    end

    heatmap(ustrip.(x),
            ustrip.(y),
            ustrip.(transpose.(I_top)),
            xaxis="x",
            yaxis="y",
            dpi=300,
            rightmargin=10Plots.mm,
            title="intensity",
            aspect_ratio=:equal,
            clim=(19.0, maximum(ustrip.(I_top))))

    savefig("../img/compare_continuum/irregular_top_I")

    return 0
end


# compare("../data/bifrost_qs006023_s525_quarter.hdf5", "../quadratures/n1.dat");
# LTE_ray("../data/bifrost_qs006023_s525_half.hdf5")
LTE_voronoi("../data/bifrost_qs006023_s525_quarter.hdf5")
# test_interpolation("../data/bifrost_qs006023_s525_quarter.hdf5", "../quadratures/n1.dat")
# test_with_regular("../data/bifrost_qs006023_s525_quarter.hdf5", "../quadratures/n1.dat")
print("")
