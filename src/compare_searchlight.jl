using Plots
using UnitfulRecipes
using NearestNeighbors

include("functions.jl")
include("voronoi_utils.jl")
include("characteristics.jl")
include("irregular_ray_tracing.jl")

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
    θ = 20*π/180
    ϕ = 15*π/180

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
    savefig("../img/compare_searchlight/irregular_SL_bottom")

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
    savefig("../img/compare_searchlight/irregular_SL_top")

    # Top to bottom
    # Traces rays through an irregular grid
    θ = 160*π/180
    ϕ = 330*π/180

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
    savefig("../img/compare_searchlight/irregular_SL_bdonw")

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
    savefig("../img/compare_searchlight/irregular_SL_tdown")

end

function searchlight_regular()
    nx = ny = nz = 51

    z = collect(LinRange(0,1,nx))u"m"
    x = collect(LinRange(0,1,ny))u"m"
    y = collect(LinRange(0,1,nz))u"m"

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

    θ_array = [175, 170]
    ϕ_array = [15,  345]

    for i in eachindex(θ_array)
        θ = θ_array[i]
        ϕ = ϕ_array[i]
        # Unit vector pointing in the direction of the ray
        k = -[cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]
        if θ > 90
            I = short_characteristics_up(θ, ϕ, S_0, α, atmos;
                            degrees=true, I_0=I_0, pt=true)[:, 2:end-1, 2:end-1]

            I = ustrip(I[end, :, :])/ustrip(I_light)
            heatmap(ustrip(x[2:end-1]),
                    ustrip(y[2:end-1]),
                    transpose(I),
                    xaxis="x",
                    yaxis="y",
                    dpi=300,
                    rightmargin=10Plots.mm,
                    title="Searchlight",
                    aspect_ratio=:equal,
                    clim=(0., 1.))

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
            plot!(circle_shape(x_r, y_r, 0.1),
                  aspect_ratio = :equal,
                  linecolor=:red,
                  lw=2)

            println("Bottom: $(I_light*80), Top: $(sum(I))")
        elseif θ < 90
            I = short_characteristics_down(θ, ϕ, S_0, α, atmos;
                            degrees=true, I_0=I_0, pt=true)[:, 2:end-1, 2:end-1]

            I = ustrip(I[1, :, :])/ustrip(I_light)
            heatmap(ustrip(x[2:end-1]),
                    ustrip(y[2:end-1]),
                    transpose(I),
                    xaxis="x",
                    yaxis="y",
                    dpi=300,
                    rightmargin=10Plots.mm,
                    title="Searchlight",
                    aspect_ratio=:equal,
                    clim=(0., 1.))

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
            plot!(circle_shape(x_r, y_r, 0.1),
                  aspect_ratio = :equal,
                  linecolor=:red,
                  lw=2)

            println("Top: $(I_light*80), Bottom: $(sum(I))")
        end
        savefig("../img/compare_searchlight/searchlight_$(θ)_$(ϕ)")
    end

    print("")
end

searchlight_regular()
searchlight_irregular()
