using NPZ
using Plots
include("plot_utils.jl")
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
    ϕ = 5.0

    n_skip = 1

    function regular()

        atmos = Atmosphere(get_atmos(DATA; periodic=true, skip=n_skip)...)

        @time J_mean, S_λ, α_tot = Λ_regular(ϵ, maxiter, atmos, quadrature)

        k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]
        I_top = short_characteristics_up(k, S_λ, α_tot, atmos, I_0=S_λ[1,:,:])

        I_top = ustrip(uconvert.(u"kW*nm^-1*m^-2", I_top[end, 2:end-1, 2:end-1]))

        global min_lim, max_lim
        min_lim = minimum(I_top)
        max_lim = maximum(I_top)

        heatmap(ustrip.(atmos.x[2:end-1]),
                ustrip.(atmos.y[2:end-1]),
                transpose.(I_top),
                xaxis="x",
                yaxis="y",
                dpi=300,
                rightmargin=10Plots.mm,
                title="Regular Grid",
                aspect_ratio=:equal,
                clim=(min_lim,max_lim))

        μ = abs(k[1])
        savefig("../img/compare_continuum/regular_cont_$(floor(Int, 100*μ))")

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
        positions = rejection_sampling(n_sites, atmos, log10.(ustrip.(atmos.hydrogen_populations)))

        # rand(3, n_sites)

        # positions[1, :] = positions[1, :].*(z_max - z_min) .+ z_min
        # positions[2, :] = positions[2, :].*(x_max - x_min) .+ x_min
        # positions[3, :] = positions[3, :].*(y_max - y_min) .+ y_min

        # positions = positions*1u"m"

        # rejection_sampling(n_sites, atmos, log10.(ustrip.(atmos.hydrogen_populations)))
        # sample_from_extinction(atmos, 500.0u"nm", n_sites)

        sites_file = "../data/sites_continuum_half.txt"
        neighbours_file = "../data/neighbours_continuum_half.txt"
        # write sites to file
        write_arrays(ustrip.(positions[2, :]),
                     ustrip.(positions[3, :]),
                     ustrip.(positions[1, :]),
                     sites_file)

        # export sites to voro++, and compute grid information
        println("---Preprocessing grid---")

        # compute neigbours
        # run(`./voro.sh $sites_file $neighbours_file
                       # $(x_min) $(x_max)
                       # $(y_min) $(y_max)
                       # $(z_min) $(z_max)`)

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
                                                                 S_λ, α_tot, 2;
                                                                 periodic=true)

        k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]
        I_top = short_characteristics_up(k, S_λ_grid, α_grid,
                                        atmos_from_voronoi, I_0=S_λ_grid[1,:,:])

        I_top = ustrip(uconvert.(u"kW*nm^-1*m^-2", I_top[end, 2:end-1, 2:end-1]))

        heatmap(ustrip.(atmos_from_voronoi.x[2:end-1]),
                ustrip.(atmos_from_voronoi.y[2:end-1]),
                transpose.(I_top),
                xaxis="x",
                yaxis="y",
                dpi=300,
                rightmargin=10Plots.mm,
                title="Irregular Grid",
                aspect_ratio=:equal,
                clim=(min_lim,max_lim))

        μ = abs(k[1])
        savefig("../img/compare_continuum/irregular_cont_$(floor(Int, 100*μ))")

        return 0
    end

    regular();
    # voronoi();

end

function LTE_compare(DATA::String, n_sites::Int)

    # atmos = Atmosphere(get_atmos(DATA; periodic=true, skip=1)...)

    # choose a wavelength
    λ = 500u"nm"  # nm

    """
    # Lte populations
    LTE_pops = LTE_populations(atmos)

    # Find continuum extinction (only with Thomson and Rayleigh)
    α_cont = α_absorption.(λ,
                           atmos.temperature,
                           atmos.electron_density*1.0,
                           LTE_pops[:,:,:,1].+LTE_pops[:,:,:,2],
                           LTE_pops[:,:,:,3]) .+
             α_scattering.(λ,
                           atmos.electron_density,
                           LTE_pops[:,:,:,1])

    # Planck function
    S_λ = blackbody_λ.(λ, atmos.temperature)

    θ = 180.0
    ϕ = 0.0

    k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]
    intensity = short_characteristics_up(k, S_λ, S_λ[1,:,:], α_cont, atmos)

    I_regular = transpose.(ustrip.(uconvert.(u"kW*nm^-1*m^-2", intensity[end, 2:end-1, 2:end-1])))
    x_regular = atmos.x[2:end-1]
    y_regular = atmos.y[2:end-1]
    """

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

    #positions = rand(3, n_sites)

    #positions[1, :] = positions[1, :].*(z_max - z_min) .+ z_min
    #positions[2, :] = positions[2, :].*(x_max - x_min) .+ x_min
    #positions[3, :] = positions[3, :].*(y_max - y_min) .+ y_min

    #positions = positions*1u"m"

    positions = sample_from_ionised_hydrogen(atmos, n_sites)
    # sample_from_extinction(atmos, λ, n_sites)
    # positions = sample_from_destruction(atmos, n_sites)
    # positions = sample_from_temp_gradient(atmos, n_sites)
    # positions = rejection_sampling(n_sites, atmos, ustrip.(atmos.temperature))


    sites_file = "../data/sites_continuum_$n_sites.txt"
    neighbours_file = "../data/neighbours_continuum_$n_sites.txt"
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

    run(`rm ../data/$sites_file ../data/$neighbours_file`)

    # Lte populations
    LTE_pops = LTE_populations(sites)

    S_λ = Matrix{Float64}(undef, (1, n_sites))u"kW*m^-2*nm^-1"
    S_λ[1,:] = blackbody_λ.(λ, sites.temperature)

    # println(typeof(S_λ))
    # println(typeof(LTE_pops))

    atmos_size = (nz, nx, ny).*2
    atmos_size = floor.(Int, atmos_size)

    atmos, S_λ_grid, populations_grid = Voronoi_to_Raster(sites, atmos_size,
                                                          S_λ, LTE_pops;
                                                          periodic=true)

    # Find continuum extinction (only with Thomson and Rayleigh)
    α_cont = α_absorption.(λ,
                           atmos.temperature,
                           atmos.electron_density*1.0,
                           populations_grid[:,:,:,1].+populations_grid[:,:,:,2],
                           populations_grid[:,:,:,3]) .+
             α_scattering.(λ,
                           atmos.electron_density,
                           populations_grid[:,:,:,1])

    θ = 180.0
    ϕ = 0.0

    k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]

    I_0 = S_λ_grid[1,1,:,:]
    intensity = short_characteristics_up(k, S_λ_grid[1,:,:,:], I_0, α_cont, atmos; n_sweeps=3)

    I_irregular = transpose.(ustrip.(uconvert.(u"kW*nm^-1*m^-2", intensity[end, 2:end-1, 2:end-1])))
    x_irregular = atmos.x[2:end-1]
    y_irregular = atmos.y[2:end-1]

    npzwrite("../data/LTE/I_irregular_$(n_sites)_ionised_hydrogen.npy", I_irregular)

end


function LTE_regular(DATA::String, n_skip::Integer)

    if n_skip == 4
        RES = "quarter"
    elseif n_skip == 3
        RES = "third"
    elseif n_skip == 2
        RES = "half"
    elseif n_skip == 1
        RES = "full"
    else
        RES = string(n_skip)
    end

    println("Continuum at $(RES) resolution, regular grid")

    atmos = Atmosphere(get_atmos(DATA; periodic=true, skip=n_skip)...)

    # choose a wavelength
    λ = 500u"nm"  # nm

    # Lte populations
    LTE_pops = LTE_populations(atmos)

    # Find continuum extinction (only with Thomson and Rayleigh)
    α_cont = α_absorption.(λ,
                           atmos.temperature,
                           atmos.electron_density*1.0,
                           LTE_pops[:,:,:,1].+LTE_pops[:,:,:,2],
                           LTE_pops[:,:,:,3]) .+
             α_scattering.(λ,
                           atmos.electron_density,
                           LTE_pops[:,:,:,1])

    # Planck function
    S_λ = blackbody_λ.(λ, atmos.temperature)

    θ = 180.0
    ϕ = 0.0

    k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]
    intensity = short_characteristics_up(k, S_λ, S_λ[1,:,:], α_cont, atmos)

    I_regular = transpose.(ustrip.(uconvert.(u"kW*nm^-1*m^-2", intensity[end, 2:end-1, 2:end-1])))
    x_regular = atmos.x[2:end-1]
    y_regular = atmos.y[2:end-1]

    npzwrite("../data/LTE/I_regular_$(RES).npy", I_regular)
    npzwrite("../data/LTE/y_regular_$(RES).npy", ustrip.(x_regular))
    npzwrite("../data/LTE/x_regular_$(RES).npy", ustrip.(y_regular))

end


function test_interpolation(DATA)
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
                   $(x_min - 0.1) $(x_max + 0.1)
                   $(y_min - 0.1) $(y_max + 0.1)
                   $(z_min - 0.1) $(z_max + 0.1)`)

    # Voronoi grid
    sites = VoronoiSites(read_cell(neighbours_file, n_sites, positions)...,
                         _initialise(positions, atmos)...,
                         z_min*1u"m", z_max*1u"m",
                         x_min*1u"m", x_max*1u"m",
                         y_min*1u"m", y_max*1u"m",
                         n_sites)

    # Lte populations
    LTE_pops = LTE_populations(sites)
    λ = 500u"nm"
    # Find continuum extinction
    α_cont = α_absorption.(λ,
                           sites.temperature,
                           sites.electron_density*1.0,
                           LTE_pops[:,1].+LTE_pops[:,2],
                           LTE_pops[:,3]) .+
             α_scattering.(λ,
                           sites.electron_density,
                           LTE_pops[:,1])

    S_λ = blackbody_λ.(λ, sites.temperature)

    θ = 175.0
    ϕ = 0.1

    k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]
    bottom_layer = sites.layers_up[2] - 1
    bottom_layer_idx = sites.perm_up[1:bottom_layer]
    println("---Ray tracing---")
    I_0 = blackbody_λ.(500u"nm", sites.temperature[bottom_layer_idx])
    intensity = Delaunay_upII(k, S_λ, α_cont, sites, I_0, 3)

    x = collect(LinRange(sites.x_min, sites.x_max, 10*nx))
    y = collect(LinRange(sites.y_min, sites.y_max, 10*ny))
    top_z = atmos.z[end]

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
            title="LTE Intensity (500 nm), Irregular Grid",
            aspect_ratio=:equal)

    savefig("../img/compare_continuum/irregular_top_I")

    atmos_from_voronoi = Voronoi_to_Raster(sites, atmos, 1; periodic=false)

    T_diff = abs.(1 .- atmos_from_voronoi.temperature ./ atmos.temperature)
    N_e_diff = abs.(1 .- atmos_from_voronoi.electron_density ./ atmos.electron_density)
    N_H_diff = abs.(1 .- atmos_from_voronoi.hydrogen_populations ./ atmos.hydrogen_populations)

    heatmap(ustrip.(atmos.x),
            ustrip.(atmos.y),
            ustrip.(transpose.(T_diff[end-20,:,:])),
            xaxis="x",
            yaxis="y",
            dpi=300,
            rightmargin=10Plots.mm,
            title="Temperature Difference",
            aspect_ratio=:equal)

    savefig("../img/compare_continuum/T_diff")

    heatmap(ustrip.(atmos.x),
            ustrip.(atmos.y),
            ustrip.(transpose.(N_e_diff[end-20,:,:])),
            xaxis="x",
            yaxis="y",
            dpi=300,
            rightmargin=10Plots.mm,
            title="Electron Density Difference",
            aspect_ratio=:equal)

    savefig("../img/compare_continuum/N_e_diff")

    heatmap(ustrip.(atmos.x),
            ustrip.(atmos.y),
            ustrip.(transpose.(N_H_diff[end-20,:,:])),
            xaxis="x",
            yaxis="y",
            dpi=300,
            rightmargin=10Plots.mm,
            title="Hydrogen Density Difference",
            aspect_ratio=:equal)

    savefig("../img/compare_continuum/N_H_diff")

    return 0
end

function test_with_regular_grid(DATA, quadrature)
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
                   $(x_min-0.1) $(x_max+0.1)
                   $(y_min-0.1) $(y_max+0.1)
                   $(z_min-0.1) $(z_max+0.1)`)

    # Voronoi grid
    sites = VoronoiSites(read_cell(neighbours_file, n_sites, positions)...,
                         temperature, electron_density, hydrogen_populations,
                         velocity_z, velocity_x, velocity_y,
                         z_min*1u"m", z_max*1u"m",
                         x_min*1u"m", x_max*1u"m",
                         y_min*1u"m", y_max*1u"m",
                         n_sites)

    # Lte populations
    LTE_pops = LTE_populations(sites)
    λ = 500u"nm"
    # Find continuum extinction (only with Thomson and Rayleigh)
    α_cont = α_absorption.(λ,
                           atmos.temperature,
                           atmos.electron_density*1.0,
                           LTE_pops[:,:,:,1].+LTE_pops[:,:,:,2],
                           LTE_pops[:,:,:,3]) .+
             α_scattering.(λ,
                           atmos.electron_density,
                           LTE_pops[:,:,:,1])
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

function compare_interpolations(DATA::String, n_sites::Int)

    atmos = Atmosphere(get_atmos(DATA; periodic=false, skip=1)...)

    # choose a wavelength
    λ = 500u"nm"  # nm

    nx = length(atmos.x)
    nz = length(atmos.z)
    ny = length(atmos.y)

    x_min = ustrip(atmos.x[1])
    x_max = ustrip(atmos.x[end])
    y_min = ustrip(atmos.y[1])
    y_max = ustrip(atmos.y[end])
    z_min = ustrip(atmos.z[1])
    z_max = ustrip(atmos.z[end])


    positions = sample_from_extinction(atmos, λ, n_sites)
    # positions = sample_from_temp_gradient(atmos, n_sites)
    # positions = rejection_sampling(n_sites, atmos, ustrip.(atmos.temperature))


    sites_file = "../data/sites_write.txt"
    neighbours_file = "../data/neighbours_write.txt"
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

    run(`rm ../data/$sites_file ../data/$neighbours_file`)

    # Lte populations
    populations = LTE_populations(sites)

    # LTE source function --> Planck function
    S_λ = Matrix{Float64}(undef, (1, n_sites))u"kW*m^-2*nm^-1"
    S_λ[1,:] = blackbody_λ.(λ, sites.temperature)

    CMAX = 60.52755753644262
    CMIN = 20.169713390042126

    atmos_size = (430, 256, 256)

    θ = 180.0
    ϕ = 0.0

    k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]

    # Inverse distance interpolation
    atmos, S_λ_grid, populations_grid = Voronoi_to_Raster_inv_dist(sites, atmos_size, S_λ, populations)

    # Find continuum extinction (only with Thomson and Rayleigh)
    α_cont = α_absorption.(λ,
                           atmos.temperature,
                           atmos.electron_density*1.0,
                           populations_grid[:,:,:,1].+populations_grid[:,:,:,2],
                           populations_grid[:,:,:,3]) .+
             α_scattering.(λ,
                           atmos.electron_density,
                           populations_grid[:,:,:,1])

    I_0 = S_λ_grid[1,1,:,:]
    intensity = short_characteristics_up(k, S_λ_grid[1,:,:,:], I_0, α_cont, atmos; n_sweeps=3)

    I_top = transpose.(ustrip.(uconvert.(u"kW*nm^-1*m^-2", intensity[end, 2:end-1, 2:end-1])))
    x = ustrip.(atmos.x[2:end-1])
    y = ustrip.(atmos.y[2:end-1])

    heatmap(x,
            y,
            I_top,
            xaxis="x",
            yaxis="y",
            dpi=300,
            rightmargin=10Plots.mm,
            title="Inverse Distance interpolation",
            aspect_ratio=:equal,
            clim=(CMIN,CMAX))

    savefig("../img/compare_continuum/inv_dist")

    # Nearest neighbour interpolation
    atmos, S_λ_grid, populations_grid = Voronoi_to_Raster(sites, atmos_size, S_λ, populations)

    # Find continuum extinction (only with Thomson and Rayleigh)
    α_cont = α_absorption.(λ,
                           atmos.temperature,
                           atmos.electron_density*1.0,
                           populations_grid[:,:,:,1].+populations_grid[:,:,:,2],
                           populations_grid[:,:,:,3]) .+
             α_scattering.(λ,
                           atmos.electron_density,
                           populations_grid[:,:,:,1])

    I_0 = S_λ_grid[1,1,:,:]
    intensity = short_characteristics_up(k, S_λ_grid[1,:,:,:], I_0, α_cont, atmos; n_sweeps=3)

    I_top = transpose.(ustrip.(uconvert.(u"kW*nm^-1*m^-2", intensity[end, 2:end-1, 2:end-1])))
    x = ustrip.(atmos.x[2:end-1])
    y = ustrip.(atmos.y[2:end-1])

    heatmap(x,
            y,
            I_top,
            xaxis="x",
            yaxis="y",
            dpi=300,
            rightmargin=10Plots.mm,
            title="Nearest Neighbour interpolation",
            aspect_ratio=:equal,
            clim=(CMIN,CMAX))

    savefig("../img/compare_continuum/nearest_neighbour")

end

# compare("../data/bifrost_qs006023_s525_half.hdf5", "../quadratures/ul2n3.dat");
n_list = [100_000, 250_000, 500_000, 1_000_000, 2_500_000, 5_000_000, 10_000_000, 15_000_000]

for num_sites in n_list
    LTE_compare("../data/bifrost_qs006023_s525.hdf5", num_sites)
end
# LTE_regular("../data/bifrost_qs006023_s525_quarter.hdf5")
# LTE_regular("../data/bifrost_qs006023_s525_half.hdf5")
# LTE_regular("../data/bifrost_qs006023_s525.hdf5", 1)
# LTE_regular("../data/bifrost_qs006023_s525.hdf5", 2)
# LTE_regular("../data/bifrost_qs006023_s525.hdf5", 3)
# LTE_regular("../data/bifrost_qs006023_s525.hdf5", 4)
# test_interpolation("../data/bifrost_qs006023_s525_half.hdf5", "../quadratures/n1.dat")
# test_with_regular("../data/bifrost_qs006023_s525_quarter.hdf5", "../quadratures/n1.dat")
# write_grid("../data/bifrost_qs006023_s525.hdf5", 10_000_000)
# compare_interpolations("../data/bifrost_qs006023_s525.hdf5", 1_000_000)
print("")
