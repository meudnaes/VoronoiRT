using .VoronoiRT
using NearestNeighbors
using Plots
using Random
using Unitful
using NPZ

pyplot()

function searchlight_irregular()
    nx = ny = nz = 51

    n_sites = nz*nx*ny

    bounds = [[0.0 1.0]
              [0.0 1.0]
              [0.0 1.0]]

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
    voro(voro_exec, sites_file, neighbours_file,
         bounds[2,1], bounds[2,2],
         bounds[3,1], bounds[3,2],
         bounds[1,1], bounds[1,2])

    # println("$(bounds[2,1]u"m"), $(bounds[2,2]u"m"), $(bounds[3,1]u"m"), $(bounds[3,2]u"m")")

    # Voronoi grid
    sites = VoronoiSites(read_cell(neighbours_file, n_sites, positions,
                                   bounds[2,1]u"m", bounds[2,2]u"m",
                                   bounds[3,1]u"m", bounds[3,2]u"m")...,
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

    tot_time = 0.0

    weights, θ_array, ϕ_array, n_angles = VoronoiRT.read_quadrature("../quadratures/ul7n12.dat")

    for i in eachindex(θ_array)
        θ = θ_array[i]
        ϕ = ϕ_array[i]
        println("$(floor(Int,θ)), $(floor(Int,ϕ))")
        # Unit vector pointing in the direction of the ray
        k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]
        if θ > 90
            I, time = @timed VoronoiRT.Delaunay_upII(k, S, I_bottom, α, sites, n_sweeps)


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
            I, time = @timed VoronoiRT.Delaunay_downII(k, S, I_top, α, sites, n_sweeps)


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
        tot_time += time
        # npzwrite("../data/searchlight_data/x_voronoi.npy", x)
        # npzwrite("../data/searchlight_data/y_voronoi.npy", y)
    end
    println("Total time $tot_time s -- avg. time: $(tot_time/12) s")
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

    _, θ_array, ϕ_array, _ = VoronoiRT.read_quadrature("../quadratures/ul7n12.dat")

    tot_time = 0.0

    for i in eachindex(θ_array)
        θ = θ_array[i]
        ϕ = ϕ_array[i]
        println("$(floor(Int,θ)), $(floor(Int,ϕ))")
        # Unit vector pointing in the direction of the ray
        k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]
        if θ > 90
            I, time = @timed VoronoiRT.short_characteristics_up(k, S_0, I_0, α, atmos;
                                         pt=true, n_sweeps=3)[:, 2:end-1, 2:end-1]

            I = I[end, :, :]
            # npzwrite("../data/searchlight_data/I_$(θ)_$(ϕ)_regular.npy", ustrip.(I))
            plot_searchlight(k, x[2:end-1], y[2:end-1], I, R0, "regular_$(floor(Int,θ))_$(floor(Int,ϕ))")
            println("Bottom: $(I_light*80), Top: $(sum(I))")
        elseif θ < 90
            I, time = @timed VoronoiRT.short_characteristics_down(k, S_0, I_0, α, atmos;
                                           pt=true, n_sweeps=3)[:, 2:end-1, 2:end-1]

            I = I[1, :, :]
            # npzwrite("../data/searchlight_data/I_$(θ)_$(ϕ)_regular.npy", ustrip.(I))
            plot_searchlight(k, x[2:end-1], y[2:end-1], I, R0, "regular_$(floor(Int,θ))_$(floor(Int,ϕ))")
            println("Top: $(I_light*80), Bottom: $(sum(I))")
        end
        tot_time += time
    end

    println("Total time $tot_time -- avg. time: $(tot_time/12)")

    print("")
end

function compare_searchlight()
    # use same recipe as Hayek et al. (2010)
    nx = ny = nz = 100
    n_sites = nz*nx*ny
    θ = 180 - 28.1
    ϕ = 45

    corner = floor(Int, 0.3*nx)

    bounds = [[0.0 1.0]
              [0.0 1.0]
              [0.0 1.0]]

    temperature = ones(n_sites)*1u"K"
    electron_density = zeros(n_sites)*1u"m^-3"
    hydrogen_density = zeros(n_sites)*1u"m^-3"

    velocity_z = zeros(n_sites)u"m/s"
    velocity_x = zeros(n_sites)u"m/s"
    velocity_y = zeros(n_sites)u"m/s"

    positions = rand(3, n_sites)u"m"

    n_sweeps = 3

    sites_file = "../data/searchlight_sites.txt"
    neighbours_file = "../data/searchlight_neighbours.txt"

    # write sites to file
    write_arrays(ustrip.(positions[2, :]),
                 ustrip.(positions[3, :]),
                 ustrip.(positions[1, :]),
                 sites_file)

    # compute neigbours
    voro(voro_exec, sites_file, neighbours_file,
         bounds[2,1], bounds[2,2],
         bounds[3,1], bounds[3,2],
         bounds[1,1], bounds[1,2])

    # Voronoi grid
    sites = VoronoiSites(read_cell(neighbours_file, n_sites, positions,
                                   bounds[2,1]u"m", bounds[2,2]u"m",
                                   bounds[3,1]u"m", bounds[3,2]u"m")...,
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
        if (xi <= 0.3u"m") & (yi <= 0.3u"m")
            I_bottom[i] = I_light
        end
    end

    x = collect(LinRange(0,1,nx))
    y = collect(LinRange(0,1,ny))

    npzwrite("../data/searchlight_data/x_100_voronoi.npy", ustrip.(x))
    npzwrite("../data/searchlight_data/y_100_voronoi.npy", ustrip.(y))

    top_z = 1.0

    # Unit vector pointing in the direction of the ray
    k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]

    # calculate intensity
    I, time = @timed VoronoiRT.Delaunay_upII(k, S, I_bottom, α, sites, n_sweeps, 50.0)
    println("voronoi time: $time")

    top_I = zeros(length(x), length(y))u"kW*nm^-1*m^-2"

    for i in 1:length(x)
        for j in 1:length(y)
            position = [top_z, x[i], y[j]]
            idx, dist = nn(tree, ustrip(position))
            top_I[i, j] = I[idx]
        end
    end

    npzwrite("../data/searchlight_data/I_100_voronoi.npy", ustrip.(top_I))

    # regular grid
    z = collect(LinRange(0,1,nz))u"m"
    x = collect(LinRange(0,1,nx))u"m"
    y = collect(LinRange(0,1,ny))u"m"

    npzwrite("../data/searchlight_data/x_100_regular.npy", ustrip.(x))
    npzwrite("../data/searchlight_data/y_100_regular.npy", ustrip.(y))

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

    I_0 = zero(S_0[1,:,:])
    I_0[1:corner, 1:corner] .= I_light

    I, time = @timed VoronoiRT.short_characteristics_up(k, S_0, I_0, α, atmos;
                                                n_sweeps=3)[:, 2:end-1, 2:end-1]

    I = I[end, :, :]
    println("regular time: $time")
    npzwrite("../data/searchlight_data/I_100_regular.npy", ustrip.(I))
    print("")
end

function do_timing(data::String, quadrature::String)

    n_skip = 2

    function time_irregular(data::String, quadrature::String)

        atmos = Atmosphere(get_atmos(data; periodic=false, skip=n_skip)...)

        nx = length(atmos.x)
        nz = length(atmos.z)
        ny = length(atmos.y)

        x_min = ustrip(atmos.x[1])
        x_max = ustrip(atmos.x[end])
        y_min = ustrip(atmos.y[1])
        y_max = ustrip(atmos.y[end])
        z_min = ustrip(atmos.z[1])
        z_max = ustrip(atmos.z[end])

        n_sites = nx*ny*nz
        println("sites: $(n_sites)")

        positions = sample_from_invNH_invT(atmos, n_sites)

        sites_file = "../data/sites_compare_$(n_sites).txt"
        neighbours_file = "../data/neighbours_compare_$(n_sites).txt"

        # write sites to file
        write_arrays(positions[2, :],
                     positions[3, :],
                     positions[1, :],
                     sites_file)

        # export sites to voro++, and compute grid information
        println("---Preprocessing grid---")

        # compute neigbours
        voro(voro_exec, sites_file, neighbours_file, x_min, x_max, y_min, y_max, z_min, z_max)

        # Voronoi grid
        sites = VoronoiSites(read_cell(neighbours_file, n_sites, positions,
                                       x_min*1u"m", x_max*1u"m",
                                       y_min*1u"m", y_max*1u"m")...,
                             initialise(positions, atmos)...,
                             z_min*1u"m", z_max*1u"m",
                             x_min*1u"m", x_max*1u"m",
                             y_min*1u"m", y_max*1u"m",
                             n_sites)

        println("---grid calculations done---")
        run(`rm ../data/$sites_file ../data/$neighbours_file`)

        weights, θ_array, ϕ_array, n_angles = VoronoiRT.read_quadrature(quadrature)

        S_λ = rand(1, length(sites.temperature))u"kW*nm^-1*m^-2"
        α_tot = rand(length(sites.temperature))u"m^-1"

        J_λ = zero(S_λ)

        n_sweeps = 3

        time = @timed begin
            for i in 1:n_angles
                θ = θ_array[i]
                ϕ = ϕ_array[i]
                k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]

                if θ_array[i] > 90
                    bottom_layer = sites.layers_up[2] - 1
                    bottom_layer_idx = sites.perm_up[1:bottom_layer]
                    I_0 = B_λ.(100.0u"nm", sites.temperature[bottom_layer_idx])
                    J_λ[1,:] += weights[i]*VoronoiRT.Delaunay_upII(k, S_λ[1,:], I_0, α_tot, sites, n_sweeps, 7.0)

                elseif θ_array[i] < 90
                    top_layer = sites.layers_down[2] - 1
                    I_0 = zeros(top_layer)u"kW*nm^-1*m^-2"
                    J_λ[1,:] += weights[i]*VoronoiRT.Delaunay_downII(k, S_λ[1,:], I_0, α_tot, sites, n_sweeps)
                end

            end
        end

        return time

    end

    function time_regular(data::String, quadrature::String)
        atmos = Atmosphere(get_atmos(DATA; periodic=true, skip=n_skip)...)

        weights, θ_array, ϕ_array, n_angles = VoronoiRT.read_quadrature(quadrature)

        S_λ = rand(1, size(atmos.temperature)...)u"kW*nm^-1*m^-2"
        α_tot = rand(size(atmos.temperature)...)u"m^-1"

        J_λ = zero(S_λ)

        time = @timed begin
            for i in 1:n_angles
                θ = θ_array[i]
                ϕ = ϕ_array[i]
                k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]

                if θ > 90
                    I_0 =  B_λ.(100.0u"nm", atmos.temperature[1,:,:])
                    J_λ[1,:,:,:] .+= weights[i].*VoronoiRT.short_characteristics_up(k,
                                                                          S_λ[1,:,:,:],
                                                                          I_0,
                                                                          α_tot,
                                                                          atmos;
                                                                          n_sweeps=3)
                elseif θ < 90
                    I_0 = zero(S_λ[1, 1, :, :])
                    J_λ[1,:,:,:] .+= weights[i].*VoronoiRT.short_characteristics_down(k,
                                                                            S_λ[1,:,:,:],
                                                                            I_0,
                                                                            α_tot,
                                                                            atmos;
                                                                            n_sweeps=3)
                end
            end
        end

        return time

    end

    t_i = time_irregular(data, quadrature)
    t_r = time_regular(data, quadrature)

    println("Irregular: $t_i")
    println("Regular: $t_r")

end

# searchlight_irregular()
# searchlight_regular()
# compare_searchlight()

DATA = "../data/bifrost_qs006023_s525.hdf5"
QUADRATURE = "../quadratures/ul7n12.dat"

do_timing(DATA, QUADRATURE)
