using Plots
using UnitfulRecipes
using NPZ

include("functions.jl")
include("voronoi_utils.jl")
include("irregular_ray_tracing.jl")

global my_seed = 1001
Random.seed!(my_seed)


function plot_sites(sites::VoronoiSites)
    pyplot()
    z = sites.positions[1,:]
    x = sites.positions[2,:]
    y = sites.positions[3,:]
    scatter(x/1e6,
            y/1e6,
            z/1e6,
            dpi=300)
    savefig("sites.png")
end

function column_count_mass(atmos::Atmosphere, sites::VoronoiSites)
    # Not the actual mass, but column number density of hydrogen...
    column_mass = Array{Unitful.Mass, 2}(undef, (nx-1, ny-1))
    column_sites = Array{Int64, 2}(undef, (nx-1, ny-1))

    z = sites.positions[1,:]
    x = sites.positions[2,:]
    y = sites.positions[3,:]

    for i in 1:nx-1
        x_bounds = (atmos.x[i], atmos.x[i+1])
        for j in 1:ny-1
            print("column $((i-1)*(ny-1)+j)/$((nx-1)*(ny-1)) \r")
            cm = 0u"kg"
            y_bounds = (atmos.y[j], atmos.y[j+1])
            for k in 1:nz-1
                cm += mass_function(k, i, j, atmos)
            end
            hits = 0
            for l in 1:n_sites
                if x_bounds[1] < x[l] < x_bounds[2] && y_bounds[1] < y[l] < y_bounds[2]
                    hits += 1
                end
            end
            column_mass[i, j] = cm
            column_sites[i, j] = hits
        end
    end
    print("                                                                 \r")

    println("$(sum(column_sites)) and $n_sites should match")
    return column_mass, column_sites
end

function searchlight_irregular_2()
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
    global sites
    sites = VoronoiSites(read_cell(neighbours_file, n_sites, positions)...,
                         temperature, electron_density, hydrogen_density,
                         velocity_z, velocity_x, velocity_y,
                         bounds[1,1]*1u"m", bounds[1,2]*1u"m",
                         bounds[2,1]*1u"m", bounds[2,2]*1u"m",
                         bounds[3,1]*1u"m", bounds[3,2]*1u"m",
                         n_sites)

    S = zeros(n_sites)u"kW*m^-2*nm^-1"
    α = zeros(n_sites)u"m^-1"

    I_light = 1u"kW*m^-2*nm^-1"

    bottom_layer = sites.layers_up[2] - 1

    I_0 = zeros(bottom_layer)u"kW*m^-2*nm^-1"
    for i in 1:bottom_layer
        idx = sites.perm_up[i]
        xi = sites.positions[2, idx]
        yi = sites.positions[3, idx]
        if sqrt((xi - 0.5u"m")^2 + (yi - 0.5u"m")^2) < R0
            I_0[i] = I_light
        end
    end

    RES=500

    # Traces rays through an irregular grid
    θ = 30*π/180
    ϕ = 5*π/180

    # start at the bottom
    # shoot rays through every grid cell

    # Unit vector towards upwind direction of the ray
    k = -[cos(θ), cos(ϕ)*sin(θ), sin(ϕ)*sin(θ)]

    println("---Ray-tracing---")
    @time I = Delaunay_up(sites, I_0, S, α, k, n_sweeps)

    bottom_x = collect(0:0.001:1)
    bottom_y = collect(0:0.001:1)
    bottom_z = 0

    bottom_I = zeros(length(bottom_x), length(bottom_y))u"kW*nm^-1*m^-2"
    tree = KDTree(ustrip(sites.positions))

    for i in 1:length(bottom_x)
        for j in 1:length(bottom_y)
            position = [bottom_z, bottom_x[i], bottom_y[j]]
            idx, dist = nn(tree, ustrip(position))
            bottom_I[i, j] = I[idx]
        end
    end

    gr()
    heatmap(bottom_x, bottom_y, transpose(ustrip(bottom_I)),
            dpi=RES, title="Beam at the Bottom", xaxis="x", yaxis="y",
            aspect_ratio= :equal)
    plot!(circle_shape(0.5, 0.5, ustrip(R0)),
          aspect_ratio = :equal,
          linecolor=:red,
          lw=2,
          label="")
    savefig("../img/voronoi2/irregular_SL_bottom")

    top_x = collect(0:0.001:1)
    top_y = collect(0:0.001:1)
    top_z = 1

    top_I = zeros(length(top_x), length(top_y))u"kW*nm^-1*m^-2"

    for i in 1:length(top_x)
        for j in 1:length(top_y)
            position = [top_z, top_x[i], top_y[j]]
            idx, dist = nn(tree, ustrip(position))
            top_I[i, j] = I[idx]
        end
    end

    x_r = 0.5 + k[2]/k[1]
    if x_r < 0
        x_r = 1 - (ceil(x_r) - x_r)
    elseif x_r > 1
        x_r = x_r - floor(x_r)
    end

    y_r = 0.5 + k[3]/k[1]
    if y_r < 0
        y_r = 1 - (ceil(y_r) - y_r)
    elseif y_r > 1
        y_r = y_r - floor(y_r)
    end

    heatmap(top_x, top_y, transpose(ustrip(top_I)),
            dpi=RES, title="Beam at the Top", xaxis="x", yaxis="y",
            right_margin = 12Plots.mm, aspect_ratio = :equal,
            clim=(0., 1.))
    plot!(circle_shape(x_r, y_r, ustrip(R0)),
          aspect_ratio = :equal,
          linecolor=:red,
          lw=2,
          label="")
    savefig("../img/voronoi2/irregular_SL_top")

    # Top to bottom
    # Traces rays through an irregular grid
    θ = 150*π/180
    ϕ = 355*π/180

    top_layer = sites.layers_down[2] - 1

    I_0 = zeros(top_layer)u"kW*m^-2*nm^-1"
    for i in 1:top_layer
        idx = sites.perm_down[i]
        xi = sites.positions[2, idx]
        yi = sites.positions[3, idx]
        if sqrt((xi - 0.5u"m")^2 + (yi - 0.5u"m")^2) < R0
            I_0[i] = I_light
        end
    end

    # Unit vector towards upwind direction of the ray
    k = -[cos(θ), cos(ϕ)*sin(θ), sin(ϕ)*sin(θ)]
    @time I = Delaunay_down(sites, I_0, S, α, k, n_sweeps)

    bottom_I = zeros(length(bottom_x), length(bottom_y))u"kW*nm^-1*m^-2"

    for i in 1:length(bottom_x)
        for j in 1:length(bottom_y)
            position = [bottom_z, bottom_x[i], bottom_y[j]]
            idx, dist = nn(tree, ustrip(position))
            bottom_I[i, j] = I[idx]
        end
    end

    x_r = 0.5 - k[2]/k[1]
    if x_r < 0
        x_r = 1 - (ceil(x_r) - x_r)
    elseif x_r > 1
        x_r = x_r - floor(x_r)
    end

    y_r = 0.5 - k[3]/k[1]
    if y_r < 0
        y_r = 1 - (ceil(y_r) - y_r)
    elseif y_r > 1
        y_r = y_r - floor(y_r)
    end

    gr()
    heatmap(bottom_x, bottom_y, transpose(ustrip(bottom_I)),
            dpi=RES, title="Beam at the Bottom, going down", xaxis="x", yaxis="y",
            aspect_ratio= :equal, clim=(0.,1.))
    plot!(circle_shape(x_r, y_r, ustrip(R0)),
          aspect_ratio = :equal,
          linecolor=:red,
          lw=2,
          label="")
    savefig("../img/voronoi2/irregular_SL_bdonw")

    top_x = collect(0:0.001:1)
    top_y = collect(0:0.001:1)
    top_z = 1

    top_I = zeros(length(top_x), length(top_y))u"kW*nm^-1*m^-2"

    for i in 1:length(top_x)
        for j in 1:length(top_y)
            position = [top_z, top_x[i], top_y[j]]
            idx, dist = nn(tree, ustrip(position))
            top_I[i, j] = I[idx]
        end
    end

    heatmap(top_x, top_y, transpose(ustrip(top_I)),
            dpi=RES, title="Beam at the Top, going down", xaxis="x", yaxis="y",
            right_margin = 12Plots.mm, aspect_ratio = :equal)
    plot!(circle_shape(0.5, 0.5, ustrip(R0)),
          aspect_ratio = :equal,
          linecolor=:red,
          lw=2,
          label="")
    savefig("../img/voronoi2/irregular_SL_tdown")

end



# main()
searchlight_irregular_2()
print("")
