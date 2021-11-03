using Plots
using NearestNeighbors

include("functions.jl")
include("lambda_iteration.jl")

global my_seed = 29
Random.seed!(my_seed)

function main()
    DATA = "../data/bifrost_qs006023_s525_quarter.hdf5"
    atmos = Atmosphere(get_atmos(DATA; periodic=true, skip=4)...)

    global nz = length(atmos.z)
    global nx = length(atmos.x)
    global ny = length(atmos.y)

    global Δx = atmos.x[2] - atmos.x[1]
    global Δy = atmos.y[2] - atmos.y[1]

    ϵ = 1e-3
    maxiter = 100
    J_mean, S_λ, α_tot = Λ_iteration(ϵ, maxiter, atmos; quadrature="../quadratures/ul2n3.dat")

    # calculate intensity at top for a small inclination
    ϕ = 170
    θ = 10
    global I_top
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
    nx = ny = nz = 100

    z = Vector(LinRange(-10, 10, nz))u"m"
    x = Vector(LinRange(-10, 10, nx))u"m"
    y = Vector(LinRange(-10, 10, ny))u"m"

    temperature = ones(nz, nx, ny)u"K"
    electron_density = zeros(nz, nx, ny)u"m^-3"
    hydrogen_populations = zeros(nz, nx, ny)u"m^-3"

    atmos = Atmosphere(z, x, y, temperature, electron_density, hydrogen_populations)

    global nz = length(atmos.z)
    global nx = length(atmos.x)
    global ny = length(atmos.y)

    global Δx = atmos.x[2] - atmos.x[1]
    global Δy = atmos.y[2] - atmos.y[1]

    S_0 = zeros(nz, nx, ny)u"kW*m^-2*nm^-1"
    α = zeros(nz, nx, ny)u"m^-1"

    I_light = 1u"kW*m^-2*nm^-1"

    I_0 = zero(S_0[1,:,:])
    I_0[45:54,45:54] .= I_light

    θ_array = [170, 120, 110, 70,  70,  10]
    ϕ_array = [10,  30,  50,  310, 330, 350]

    for (θ, ϕ) in zip(θ_array, ϕ_array)

        if θ > 90
            I = short_characteristics_up(θ, ϕ, S_0, α, atmos;
                            degrees=true, I_0=I_0, pt=true)[:, 2:end-1, 2:end-1]
            heatmap(ustrip(x[2:end-1]),
                    ustrip(y[2:end-1]),
                    ustrip(I[end, :, :])/ustrip(I_light),
                    xaxis="x",
                    yaxis="y",
                    dpi=300,
                    rightmargin=10Plots.mm,
                    title="Searchlight",
                    aspect_ratio=:equal)

            println("Bottom: $(I_light*100), Top: $(sum(I[end,:,:]))")
        elseif θ < 90
            I = short_characteristics_down(θ, ϕ, S_0, α, atmos;
                            degrees=true, I_0=I_0, pt=true)[:, 2:end-1, 2:end-1]
            heatmap(ustrip(x[2:end-1]),
                    ustrip(y[2:end-1]),
                    ustrip(I[1, :, :])/ustrip(I_light),
                    xaxis="x",
                    yaxis="y",
                    dpi=300,
                    rightmargin=10Plots.mm,
                    title="Searchlight",
                    aspect_ratio=:equal)

            println("Top: $(I_light*100), Bottom: $(sum(I[1,:,:]))")
        end
        savefig("../img/searchlight/searchlight_$(θ)_$(ϕ)")
    end

    print("")
end

searchlight()
# main()
