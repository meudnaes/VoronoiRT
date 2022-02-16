using Plots
include("lambda_continuum.jl")

global my_seed = 1998
Random.seed!(my_seed)

function compare(DATA, quadrature)
    maxiter = 50
    ϵ = 1e-4

    θ = 0.0
    ϕ = 0.0

    n_skip = 1

    function regular()

        atmos = Atmosphere(get_atmos(DATA; periodic=true, skip=n_skip)...)

        J_mean, S_λ, α_tot = Λ_regular(ϵ, maxiter, atmos, quadrature)

        k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]
        I_top = short_characteristics_up(k, S_λ, α_tot, atmos, I_0=S_λ[1,:,:])

        I_top = ustrip(uconvert.(u"kW*nm^-1*m^-2", I_top[end, 2:end-1, 2:end-1]))

        # global min_lim, max_lim
        # min_lim = minimum(I_top)
        # max_lim = maximum(I_top)

        heatmap(ustrip(atmos.x[2:end-1]),
                ustrip(atmos.y[2:end-1]),
                transpose(I_top),
                xaxis="x",
                yaxis="y",
                dpi=300,
                rightmargin=10Plots.mm,
                title="Regular Grid",
                aspect_ratio=:equal)
                # clim=(min_lim,max_lim))

        savefig("../img/compare_continuum/regular_top_n2")

        return 0
    end


    function voronoi()

        atmos = Atmosphere(get_atmos(DATA; periodic=false, skip=n_skip)...)

        nx = length(atmos.x)
        nz = length(atmos.z)
        ny = length(atmos.y)

        n_sites = floor(Int, nz*nx*ny)
        positions = rejection_sampling(n_sites, atmos, log10.(ustrip.(atmos.hydrogen_populations)))

        sites_file = "../data/sites_continuum.txt"
        neighbours_file = "../data/neighbours_continuum.txt"
        # write sites to file
        write_arrays(ustrip.(positions[2, :]),
                     ustrip.(positions[3, :]),
                     ustrip.(positions[1, :]),
                     sites_file)

        x_min = ustrip(atmos.x[1])
        x_max = ustrip(atmos.x[end])
        y_min = ustrip(atmos.y[1])
        y_max = ustrip(atmos.y[end])
        z_min = ustrip(atmos.z[1])
        z_max = ustrip(atmos.z[end])

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

        J_mean, S_λ, α_tot = Λ_voronoi(ϵ, maxiter, sites, quadrature)

        atmos_from_voronoi, S_λ_grid, α_grid = Voronoi_to_Raster(sites, atmos, S_λ, α_tot, 1.5)

        k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]
        I_top = short_characteristics_up(k, S_λ_grid, α_grid, atmos_from_voronoi, I_0=S_λ_grid[1,:,:])

        I_top = ustrip(uconvert.(u"kW*nm^-1*m^-2", I_top[end, 2:end-1, 2:end-1]))

        heatmap(ustrip(atmos_from_voronoi.x[2:end-1]),
                ustrip(atmos_from_voronoi.y[2:end-1]),
                transpose(I_top),
                xaxis="x",
                yaxis="y",
                dpi=300,
                rightmargin=10Plots.mm,
                title="Irregular Grid",
                aspect_ratio=:equal)
                # clim=(min_lim,max_lim))

        savefig("../img/compare_continuum/irregular_top_n2")

        return 0
    end

    regular();
    voronoi();

end

function LTE_ray(DATA)

    atmos = Atmosphere(get_atmos(DATA; periodic=false, skip=1)...)
    line = HydrogenicLine(test_atom(1, 1)..., atmos.temperature)

    # choose a wavelength
    λ = 500u"nm"  # nm

    # Lte populations
    LTE_pops = LTE_ionisation(atmos)

    # Find continuum extinction (only with Thomson and Rayleigh)
    α_cont = α_continuum.(λ, atmos.temperature*1.0, atmos.electron_density*1.0,
                          LTE_pops[:,:,:,1]*1.0, LTE_pops[:,:,:,3]*1.0)

    # Planck function
    B_λ = blackbody_λ.(λ, atmos.temperature)

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

    I_top = ustrip(uconvert.(u"kW*nm^-1*m^-2", intensity[end, :, :]))

    heatmap(ustrip(atmos.x),
            ustrip(atmos.y),
            transpose(I_top),
            xaxis="x",
            yaxis="y",
            dpi=300,
            rightmargin=10Plots.mm,
            title="Continuum at 500 nm",
            aspect_ratio=:equal)

    savefig("../img/compare_continuum/cont500")

    return 0
end



compare("../data/bifrost_qs006023_s525_quarter.hdf5", "../quadratures/n2.dat");
# LTE_ray("../data/bifrost_qs006023_s525_quarter.hdf5")
print("")
