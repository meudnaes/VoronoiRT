using Plots
using NearestNeighbors

gr()

include("functions.jl")
include("lambda_iteration.jl")

global my_seed = 29
Random.seed!(my_seed)

function main()
    DATA = "../data/bifrost_qs006023_s525_quarter.hdf5"
    atmos = Atmosphere(get_atmos(DATA; periodic=true, skip=4)...)

    ϵ = 1e-3
    maxiter = 100
    J_mean, S_λ, α_tot = Λ_iteration(ϵ, maxiter, atmos; quadrature="../quadratures/ul2n3.dat")

    println(sum(isnan.(S_λ[:, 2:end-1, 2:end-1])))

    # calculate intensity at top for a small inclination
    ϕ = 170
    θ = 10

    I_top = short_characteristics_up(θ, ϕ, S_λ, α_tot, atmos; degrees=true)
    pyplot()
    heatmap(ustrip(atmos.x[2:end-1]),
            ustrip(atmos.y[2:end-1]),
            ustrip(I_top[end,2:end-1, 2:end-1]),
            c = :solar,
            dpi=300)
    savefig("intensity_500")

    #println(isnan.(J[:, 2:nx-1, 2:ny-1]))
    #=
    p_unitless = ustrip(p_vec)

    tree = KDTree(p_unitless; leafsize = 10)

    @time begin
        index, dist = nn(tree, [ustrip(z0), ustrip(x0), ustrip(y0)])
    end

    println("z0: $z0, x0: $x0, y0: $y0")
    println("z: $(p_vec[1, index]), x: $(p_vec[2, index]), y: $(p_vec[3, index])")
    println("Distance: $dist")
    =#

    #=
    pyplot()
    plot(mass,
         sites_enclosed,
         seriestype=:scatter,
         xscale=:identity,
         yscale=:identity,
         xlabel="cube mass",
         ylabel="sites inside cube",
         size=(800, 500),
         dpi=300,
         ylims=[0.1, :auto],
         legend=:bottomright)
    savefig("the_new1.png")
    =#
    print("")
end

function searchlight()
    nx = ny = nz = 50

    z = collect(LinRange(0,1,nx))u"m"
    x = collect(LinRange(0,1,ny))u"m"
    y = collect(LinRange(0,1,nz))u"m"

    temperature = ones(nz, nx, ny)u"K"
    electron_density = zeros(nz, nx, ny)u"m^-3"
    hydrogen_populations = zeros(nz, nx, ny)u"m^-3"

    atmos = Atmosphere(z, x, y, temperature, electron_density, hydrogen_populations)

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
                @inbounds I_0[i, j] = I_light
            end
        end
    end

    θ_array = [170, 120, 110, 70,  70,  10]
    ϕ_array = [10,  30,  50,  310, 330, 350]

    for (θ, ϕ) in zip(θ_array, ϕ_array)
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
                    aspect_ratio=:equal)

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

            println("Bottom: $(I_light*100), Top: $(sum(I))")
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
                    aspect_ratio=:equal)

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

            println("Top: $(I_light*100), Bottom: $(sum(I))")
        end
        savefig("../img/searchlight/searchlight_$(θ)_$(ϕ)")
    end

    print("")
end

# searchlight()
main()
