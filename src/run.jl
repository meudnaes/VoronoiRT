using HDF5
using Plots
using Printf
using Random
using NearestNeighbors

include("characteristics.jl")
include("functions.jl")

global my_seed = 29
Random.seed!(my_seed)

function run()
    DATA = "../data/bifrost_qs006023_s525_quarter.hdf5"
    atmos = Atmosphere(get_atmos(DATA)...)

    nz = length(atmos.z) - 1
    nx = length(atmos.x) - 1
    ny = length(atmos.y) - 1

    n_sites = 1_000 #nz*nx*ny

    #z_new, x_new, y_new
    p_vec, N_H_new = rejection_sampling(n_sites, atmos)

    #println("$(size(p_vec)[2]) new sites")

    #=
    pyplot()
    scatter(x_new/1e6,
            y_new/1e6,
            z_new/1e6,
            dpi=300)
    savefig("sites.png")
    =#

    mass = Array{Float32, 1}(undef, nx*ny*nz)
    sites_enclosed = Array{Int64, 1}(undef, nx*ny*nz)

    # sort z array
    P = sortperm(p_vec[1, :])
    p_vec = p_vec[:, P]
    # z_new = z_new[P]
    # x_new = x_new[P]
    # y_new = y_new[P]

    #=
    @time begin
        for k in 1:nz
            for i in 1:nx
                for j in 1:ny
                    mass[(k-1)*nx*ny + (i-1)*ny + j] = mass_function(k, i, j, atmos)
                    z_bounds = (atmos.z[k], atmos.z[k+1])
                    x_bounds = (atmos.x[i], atmos.x[i+1])
                    y_bounds = (atmos.y[j], atmos.y[j+1])
                    sites_enclosed[(k-1)*nx*ny + (i-1)*ny + j] = find_sites_sorted(
                                    z_new, x_new, y_new, z_bounds, x_bounds, y_bounds)
                end
            end
        end
    end
    =#

    # Print boundaries for voro++ file
    #=
    open("boundaries.txt", "w") do io
        println(io, "x_min = $(atmos.x[1])")
        println(io, "x_max = $(atmos.x[end])")
        println(io, "y_min = $(atmos.y[1])")
        println(io, "y_max = $(atmos.y[end])")
        println(io, "z_min = $(atmos.z[1])")
        println(io, "z_max = $(atmos.z[end])")
    end
    =#

    write_arrays(p_vec[1, :], p_vec[2, :], p_vec[3, :], "seeds.txt")

    # initialise the system in LTE
    temperature_new = Vector{Unitful.Temperature}(undef, n_sites)
    N_e_new = Vector{NumberDensity}(undef, n_sites)
    for k in 1:n_sites
        temperature_new[k] = trilinear(p_vec[1, k], p_vec[2, k], p_vec[3, k], atmos, atmos.temperature)
        N_e_new[k] = trilinear(p_vec[1, k], p_vec[2, k], p_vec[3, k], atmos, atmos.electron_density)
        # temperature_new[k] = trilinear(z_new[k], x_new[k], y_new[k], atmos, atmos.temperature)
        # N_e_new[k] = trilinear(z_new[k], x_new[k], y_new[k], atmos, atmos.electron_density)
    end

    # choose a wavelength
    λ = 500u"nm"  # nm

    # Only continuum
    η_ν = 0

    # constant destruction (simplification)
    ε_ν = 0.1

    # Find continuum extinction (only with Thomson and Rayleigh)
    α = αcont.(λ, atmos.temperature, atmos.electron_density,
                atmos.hydrogen_populations, atmos.hydrogen_populations)

    # Find τ continuum
    # τ_c = cumul_integrate(atmos.z, α)

    # Start with the source function as the Planck function
    S_0 = B_λ.(λ, atmos.temperature[:,:,:])

    # Do the Λ iteration to get the source function iteratively
    # S_new = Λ_iteration(τ_c, S_0, η_ν, ε_ν)[1]
    # S_new = uconvert.(u"kW*m^-2*nm^-1",S_new)

    # do a quadrature of rays in a grid cell ...
    # start at index 2, since we know the downwind intensity
    # idz = 2, idx = 10, idy = 10 ...

    # first ray kind of arbitrarily chosen
    # Ω = (θ, φ), space angle
    weights, θ_array, ϕ_array, n_points = read_quadrature("../quadratures/ul9n20.dat")

    # start shooting rays from the bottom of the domain. I have to start at the
    # bottom because the rays move upwards. Kind of obvious, but it took some
    # time to understand. If rays move downwards, I start at the top...
    # NB we start at 2nd index because 1st index is boundary condition
    #I = zeros(length(atmos.z), length(atmos.x), length(atmos.y))*u"kW*m^-2*nm^-1"
    # I_0 = S_0

    ############################ Test xy_up ####################################
    θ = 10
    ϕ = 10
    θ = θ*pi/180
    ϕ = ϕ*pi/180
    print("ray 1: ")
    @time I = short_characteristics_up(θ, ϕ, S_0, α, atmos, degrees=false)
    ########################### Test xy_down ###################################
    θ = 170
    ϕ = 10
    θ = θ*pi/180
    ϕ = ϕ*pi/180
    print("ray 2: ")
    @time I = short_characteristics_down(θ, ϕ, S_0, α, atmos, degrees=false)
    ############################ Test yz_up ####################################
    θ = 60
    ϕ = 10
    θ = θ*pi/180
    ϕ = ϕ*pi/180
    print("ray 3: ")
    @time I = short_characteristics_up(θ, ϕ, S_0, α, atmos, degrees=false)
    ########################### Test yz_down ###################################
    θ = 110
    ϕ = 10
    θ = θ*pi/180
    ϕ = ϕ*pi/180
    print("ray 4: ")
    @time I = short_characteristics_down(θ, ϕ, S_0, α, atmos, degrees=false)

    # for i in 1:n_points
        # print("Ray $i/$n_points \r")
        # I = short_characteristic_ray(θ_array[i], ϕ_array[i], S_0, α, atmos, degrees=true)
    # end

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
    plot(ustrip(atmos.z), ustrip(S_new), dpi=300)
    savefig("S_regular_500.png")

    I = zeros(length(atmos.x), length(atmos.y))u"kW*m^-2*nm^-1"

    for i=1:length(atmos.x), j=1:length(atmos.y)
        # Find continuum extinction (only with Thomson and Rayleigh)
        α = αcont.(λ, atmos.temperature[:,i,j], atmos.electron_density[:,i,j],
                    atmos.hydrogen_populations[:,i,j], atmos.hydrogen_populations[:,i,j])

        # Find τ continuum
        τ_c = cumul_integrate(atmos.z, α)

        # Start with the source function as the Planck function
        S_0 = B_λ.(λ, atmos.temperature[:,i,j])

        # Do the Λ iteration to get the source function iteratively
        S_new = Λ_iteration(τ_c, S_0, η_ν, ε_ν)[1]

        # integrate to find intensity at the top of the box
        I[i, j] = intensity(S_0[1], S_new, τ_c)
        progress = (length(atmos.x)*i+j)/(length(atmos.x)*(length(atmos.y)+1))
        @printf("progress: %.2f\r", progress)
    end

    pyplot()
    heatmap(ustrip(atmos.x), ustrip(atmos.y), ustrip(I), c = :solar, dpi=300)
    savefig("intensity_500.png")
    =#

    # plot(ustrip(atmos.z), τ_c, dpi=300)
    # savefig("tau_c.png")

    #println(S_new)

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
