using Plots

include("line.jl")
include("functions.jl")
include("atmosphere.jl")
include("lambda_iteration.jl")

global my_seed = 2022
Random.seed!(my_seed)

function compare(DATA, quadrature)
    maxiter = 100
    ϵ = 1e-3

    θ = 10
    ϕ = 10

    function regular()
        atmos = Atmosphere(get_atmos(DATA; periodic=true, skip=6)...)
        line = HydrogenicLine(test_atom()..., atmos.temperature)

        global J_mean, S_λ, α_cont, populations
        J_mean, S_λ, α_cont, populations = Λ_regular(ϵ, maxiter, atmos, line, quadrature)

        ΔD = doppler_width.(line.λ0, line.atom_weight, atmos.temperature)
        γ = γ_constant(line,
                       atmos.temperature,
                       (populations[:, :, :, 1].+populations[:, :, :, 2]),
                       atmos.electron_density)

        # a = damping_constant.(γ, ΔD)
        damping_λ = damping

        k = -[cos(θ), cos(ϕ)*sin(θ), sin(ϕ)*sin(θ)]

        v_los = line_of_sight_velocity(atmos, k)

        v = Array{Float64, 4}(undef, (length(line.λ), size(v_los)...))
        for l in eachindex(line.λ)
            v[l, :, :, :] = (line.λ[l] .- line.λ0 .+ line.λ0.*v_los./c_0)./ΔD .|> Unitful.NoUnits
        end

        profile = Array{PerLength, 4}(undef, (length(line.λ), size(v_los)...))
        for l in eachindex(line.λ)
            damping_λ = damping.(γ, line.λ[l], line.ΔD)
            profile[l, :, :, :] = voigt_profile.(damping_λ, v[l, :, :, :], ΔD)
        end

        α_line = Array{Float64, 4}(undef, size(profile))u"m^-1"
        for l in eachindex(line.λ)
            α_line[l, :, :, :] = αline_λ(line,
                                         profile[l, :, :, :],
                                         populations[:, :, :, 1],
                                         populations[:, :, :, 2])
        end
        α_tot = α_line .+ α_cont

        I_top = short_characteristics_up(θ, ϕ, S_λ[6,:,:,:], α_tot[6,:,:,:],
                                            atmos, degrees=true, I_0=S_λ[6,1,:,:])

        I_top = ustrip(uconvert.(u"kW*nm^-1*m^-2", I_top[end, 2:end-1, 2:end-1]))

        heatmap(ustrip(atmos.x[2:end-1]),
                ustrip(atmos.y[2:end-1]),
                transpose(I_top),
                xaxis="x",
                yaxis="y",
                dpi=300,
                rightmargin=10Plots.mm,
                title="Regular Grid",
                aspect_ratio=:equal)#,
                #clim=(1.0,15.0))

        savefig("../img/compare_line/regular_top")

        return atmos, S_λ
    end


    function voronoi(atmos::Atmosphere)

        nx = length(atmos.x)
        nz = length(atmos.z)
        ny = length(atmos.y)

        n_sites = floor(Int, nz*nx*ny/8)
        positions = rejection_sampling(n_sites, atmos, log10.(ustrip.(atmos.hydrogen_populations)))

        sites_file = "../data/sites_compare.txt"
        neighbours_file = "../data/neighbours_compare.txt"
        # write sites to file
        write_arrays(ustrip(positions[2, :]),
                     ustrip(positions[3, :]),
                     ustrip(positions[1, :]),
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

        atmos_from_voronoi, S_λ_grid, α_grid = Voronoi_to_Raster(sites, atmos, S_λ, α_tot, 3)

        I_top = short_characteristics_up(θ, ϕ, S_λ_grid,
                                         α_grid, atmos_from_voronoi, I_0=S_λ_grid[1,:,:])

        I_top = ustrip(uconvert.(u"kW*nm^-1*m^-2", I_top[end, 2:end-1, 2:end-1]))

        heatmap(ustrip(atmos_from_voronoi.x[2:end-1]),
             ustrip(atmos_from_voronoi.y[2:end-1]),
             transpose(I_top),
             xaxis="x",
             yaxis="y",
             dpi=300,
             rightmargin=10Plots.mm,
             title="Irregular Grid",
             aspect_ratio=:equal,
             clim=(1.0,15.0))

        savefig("../img/compare_converged/irregular_top")

        return atmos_from_voronoi, S_λ_grid
    end

    atmos, S_regular = regular();
    # atmos_voronoi, S_voronoi = voronoi(atmos);

end

DATA = "../data/bifrost_qs006023_s525_quarter.hdf5"
QUADRATURE = "../quadratures/ul7n12.dat"

compare(DATA, QUADRATURE);
print("")
