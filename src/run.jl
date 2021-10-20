using Plots
using NearestNeighbors

include("functions.jl")
include("lambda_iteration.jl")

global my_seed = 29
Random.seed!(my_seed)

function run()
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

run()
