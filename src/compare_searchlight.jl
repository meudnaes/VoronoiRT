using NPZ
using Plots
using NearestNeighbors

include("functions.jl")
include("plot_utils.jl")
include("voronoi_utils.jl")
include("characteristics.jl")
include("irregular_ray_tracing.jl")

pyplot()

function searchlight_irregular()
    nx = ny = nz = 51

    n_sites = nz*nx*ny

    bounds = [[0 1]
              [0 1]
              [0 1]]

    temperature = ones(n_sites)*1u"K"
    electron_density = zeros(n_sites)*1u"m^-3"
    hydrogen_density = zeros(n_sites)*1u"m^-3"

    velocity_z = zeros(n_sites)u"m/s"
    velocity_x = zeros(n_sites)u"m/s"
    velocity_y = zeros(n_sites)u"m/s"

    println("---Computing grid---")
    positions = rand(3, n_sites)u"m"

    v0 = [0, 0.5, 0.5]
    R0 = 0.1u"m"

    n_sweeps = 3
    # positions = sample_beam(n_sites, bounds, beam, v0, R0, k)

    sites_file = "../data/searchlight_sites.txt"
    neighbours_file = "../data/searchlight_neighbours.txt"

    # write sites to file
    write_arrays(ustrip.(positions[2, :]),
                 ustrip.(positions[3, :]),
                 ustrip.(positions[1, :]),
                 sites_file)

    # compute neigbours
    @time begin
    run(`./voro.sh $sites_file $neighbours_file
            $(bounds[2,1]) $(bounds[2,2])
            $(bounds[3,1]) $(bounds[3,2])
            $(bounds[1,1]) $(bounds[1,2])`)
    end

    # Voronoi grid
    sites = VoronoiSites(read_cell(neighbours_file, n_sites, positions)...,
                         temperature, electron_density, hydrogen_density,
                         velocity_z, velocity_x, velocity_y,
                         bounds[1,1]*1u"m", bounds[1,2]*1u"m",
                         bounds[2,1]*1u"m", bounds[2,2]*1u"m",
                         bounds[3,1]*1u"m", bounds[3,2]*1u"m",
                         n_sites)

    S = zeros(n_sites)u"kW*m^-2*nm^-1"
    α = zeros(n_sites)u"m^-1"

    tree = KDTree(ustrip(sites.positions))

    I_light = 1u"kW*m^-2*nm^-1"

    bottom_layer = sites.layers_up[2] - 1

    I_bottom = zeros(bottom_layer)u"kW*m^-2*nm^-1"
    for i in 1:bottom_layer
        idx = sites.perm_up[i]
        xi = sites.positions[2, idx]
        yi = sites.positions[3, idx]
        if sqrt((xi - 0.5u"m")^2 + (yi - 0.5u"m")^2) < R0
            I_bottom[i] = I_light
        end
    end

    top_layer = sites.layers_down[2] - 1

    I_top = zeros(top_layer)u"kW*m^-2*nm^-1"
    for i in 1:top_layer
        idx = sites.perm_down[i]
        xi = sites.positions[2, idx]
        yi = sites.positions[3, idx]
        if sqrt((xi - 0.5u"m")^2 + (yi - 0.5u"m")^2) < R0
            I_top[i] = I_light
        end
    end

    x = collect(LinRange(0,1,10*nx))
    y = collect(LinRange(0,1,10*ny))
    top_z = 1.0
    bottom_z = 0.0

    RES=500

    weights, θ_array, ϕ_array, n_angles = read_quadrature("../quadratures/ul7n12.dat")

    for i in eachindex(θ_array)
        θ = θ_array[i]
        ϕ = ϕ_array[i]
        println("$(floor(Int,θ)), $(floor(Int,ϕ))")
        # Unit vector pointing in the direction of the ray
        k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]
        if θ > 90
            @time I = Delaunay_upII(k, S, α, sites, I_bottom, n_sweeps)

            top_I = zeros(length(x), length(y))u"kW*nm^-1*m^-2"

            for i in 1:length(x)
                for j in 1:length(y)
                    position = [top_z, x[i], y[j]]
                    idx, dist = nn(tree, ustrip(position))
                    top_I[i, j] = I[idx]
                end
            end

            plot_searchlight(k, x, y, top_I, R0, "irregular_$(floor(Int,θ))_$(floor(Int,ϕ))")
        elseif θ < 90
            @time I = Delaunay_downII(k, S, α, sites, I_top, n_sweeps)

            bottom_z = 0
            bottom_I = zeros(length(x), length(y))u"kW*nm^-1*m^-2"

            for i in 1:length(x)
                for j in 1:length(y)
                    position = [bottom_z, x[i], y[j]]
                    idx, dist = nn(tree, ustrip(position))
                    bottom_I[i, j] = I[idx]
                end
            end

            # npzwrite("../data/searchlight_data/I_$(θ)_$(ϕ)_voronoi.npy", ustrip.(bottom_I))
            plot_searchlight(k, x, y, bottom_I, R0, "irregular_$(floor(Int,θ))_$(floor(Int,ϕ))")
        end
        # npzwrite("../data/searchlight_data/x_voronoi.npy", x)
        # npzwrite("../data/searchlight_data/y_voronoi.npy", y)
    end
end

function searchlight_regular()
    nx = ny = nz = 51

    z = collect(LinRange(0,1,nz))u"m"
    x = collect(LinRange(0,1,nx))u"m"
    y = collect(LinRange(0,1,ny))u"m"

    # npzwrite("../data/searchlight_data/x_regular.npy", ustrip.(x))
    # npzwrite("../data/searchlight_data/y_regular.npy", ustrip.(y))

    temperature = ones(nz, nx, ny)u"K"
    electron_density = zeros(nz, nx, ny)u"m^-3"
    hydrogen_density = zeros(nz, nx, ny)u"m^-3"

    velocity_z = zeros(nz, nx, ny)u"m/s"
    velocity_x = zeros(nz, nx, ny)u"m/s"
    velocity_y = zeros(nz, nx, ny)u"m/s"

    atmos = Atmosphere(z, x, y, temperature, electron_density, hydrogen_density,
                       velocity_z, velocity_x, velocity_y)

    S_0 = zeros(nz, nx, ny)u"kW*m^-2*nm^-1"
    α = zeros(nz, nx, ny)u"m^-1"

    I_light = 1u"kW*m^-2*nm^-1"

    I_0 = zero(S_0[1,:,:])
    R0 = 0.1
    for i in 1:nx
        for j in 1:ny
            xi = i/nx
            yi = j/ny
            if sqrt((xi - 0.5)^2 + (yi - 0.5)^2) < R0
                I_0[i, j] = I_light
            end
        end
    end

    weights, θ_array, ϕ_array, n_angles = read_quadrature("../quadratures/ul7n12.dat")

    for i in eachindex(θ_array)
        θ = θ_array[i]
        ϕ = ϕ_array[i]
        println("$(floor(Int,θ)), $(floor(Int,ϕ))")
        # Unit vector pointing in the direction of the ray
        k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]
        if θ > 90
            I = short_characteristics_up(k, S_0, α, atmos;
                                         I_0=I_0, pt=true, n_sweeps=3)[:, 2:end-1, 2:end-1]

            I = I[end, :, :]
            # npzwrite("../data/searchlight_data/I_$(θ)_$(ϕ)_regular.npy", ustrip.(I))
            plot_searchlight(k, x[2:end-1], y[2:end-1], I, R0, "regular_$(floor(Int,θ))_$(floor(Int,ϕ))")
            println("Bottom: $(I_light*80), Top: $(sum(I))")
        elseif θ < 90
            I = short_characteristics_down(k, S_0, α, atmos;
                                           I_0=I_0, pt=true, n_sweeps=3)[:, 2:end-1, 2:end-1]

            I = I[1, :, :]
            # npzwrite("../data/searchlight_data/I_$(θ)_$(ϕ)_regular.npy", ustrip.(I))
            plot_searchlight(k, x[2:end-1], y[2:end-1], I, R0, "regular_$(floor(Int,θ))_$(floor(Int,ϕ))")
            println("Top: $(I_light*80), Bottom: $(sum(I))")
        end
    end

    print("")
end

searchlight_irregular()
searchlight_regular()
