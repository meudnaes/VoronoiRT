using NPZ
using Plots
using UnitfulRecipes

include("line.jl")
include("functions.jl")
include("radiation.jl")
include("atmosphere.jl")
include("broadening.jl")
include("populations.jl")
include("characteristics.jl")

pyplot()

"""
    function circle_shape(x, y, r)

Function for a circle with centre coordinates (x, y) and radius r.
"""
function circle_shape(x, y, r)
    θ = LinRange(0, 2π, 500)
    x .+ r*cos.(θ), y .+ r*sin.(θ)
end

function plot_searchlight(k::Vector{Float64},
                          x,
                          y,
                          I,
                          R0,
                          title::String)

    @assert norm(k) ≈ 1 "Wrong direction vector!"

    RES = 300
    x_r = 0.5 - sign(k[1])*k[2]/k[1]
    if x_r < 0
        x_r = 1 - (ceil(x_r) - x_r)
    elseif x_r > 1
        x_r = x_r - floor(x_r)
    end

    y_r = 0.5 - sign(k[1])*k[3]/k[1]
    if y_r < 0
        y_r = 1 - (ceil(y_r) - y_r)
    elseif y_r > 1
        y_r = y_r - floor(y_r)
    end

    heatmap(ustrip.(x), ustrip.(y), ustrip.(transpose(I)),
            dpi=RES, title=title, xaxis="x", yaxis="y",
            aspect_ratio= :equal, clim=(0.,1.))
    plot!(circle_shape(ustrip(x_r), ustrip(y_r), ustrip(R0)),
          aspect_ratio = :equal,
          linecolor=:red,
          lw=2,
          label="")
    savefig("../img/compare_searchlight/$title")

end

function plot_top_intensity(atmos::Atmosphere,
                            line::HydrogenicLine,
                            S_λ::Array{<:UnitsIntensity_λ, 4},
                            α_tot::Array{<:PerLength, 4},
                            θ::Float64,
                            ϕ::Float64,
                            idλ::Int,
                            title::String)

    k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]
    I_top = short_characteristics_up(k, S_λ[idλ, :, :, :], S_λ[idλ, 1, :, :],
                                     α_tot[idλ, :, :, :], atmos)

    I_top = ustrip(uconvert.(u"kW*nm^-1*m^-2", I_top[end, 2:end-1, 2:end-1]))

    l1 = line.λidx[1]+1
    l2 = line.λidx[2]
    if l1 <= idλ <= l2
        colors = :gist_yarg
    else
        colors = :magma
    end

    heatmap(ustrip(atmos.x[2:end-1]),
            ustrip(atmos.y[2:end-1]),
            transpose(I_top),
            xaxis="x",
            yaxis="y",
            dpi=300,
            rightmargin=10Plots.mm,
            title=title*" at $(round(ustrip(line.λ[idλ]); digits=4)) nm",
            aspect_ratio=:equal,
            c=colors,
            clim=(5,200))

    savefig("../img/compare_line/"*title)
end

function write_top_intensity(atmos::Atmosphere,
                             line::HydrogenicLine,
                             S_λ::Array{<:UnitsIntensity_λ, 4},
                             α_tot::Array{<:PerLength, 4},
                             θ::Float64,
                             ϕ::Float64,
                             fname::String)

    l1 = line.λidx[1]+1
    l2 = line.λidx[2]

    I_top = Array{Float64, 3}(undef, (l2-l1+1,
                                      length(atmos.x[2:end-1]),
                                      length(atmos.y[2:end-1])))

    k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]
    println("--Calculating intensity---")
    Threads.@threads for idλ in l1:l2
        intensity = short_characteristics_up(k, S_λ[idλ, :, :, :], S_λ[idλ, 1, :, :],
                                         α_tot[idλ, :, :, :], atmos)

        intensity = transpose(intensity[end, 2:end-1, 2:end-1])
        I_top[idλ,:,:] = ustrip(uconvert.(u"kW*nm^-1*m^-2", intensity))
    end

    println("---Writing to file---")
    npzwrite("../python/linedata/$(fname).npy", I_top)
end

function plot_top_line(atmos::Atmosphere,
                       line::HydrogenicLine,
                       S_λ::Array{<:UnitsIntensity_λ, 4},
                       α_tot::Array{<:PerLength, 4},
                       θ::Float64,
                       ϕ::Float64,
                       title::String)

    nλ, nz, nx, ny = size(S_λ)
    I_λ = zero(S_λ)

    k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]
    for iλ in 1:nλ
        I_λ[iλ, :, :, :] = short_characteristics_up(k, S_λ[iλ, :, :, :], S_λ[iλ, 1, :, :],
                                                    α_tot[iλ, :, :, :], atmos)
    end

    start = line.λidx[1]+1
    stop =  line.λidx[2]

    indices = sortperm(line.λ[start:stop])

    for idx in 1:10:size(S_λ)[end]
        for idy in 1:10:size(S_λ)[end-1]
            loc = "_"*string(idx)*"_"*string(idy)
            I_plot = ustrip(uconvert.(u"kW*nm^-1*m^-2", I_λ))[:, end, idx, idy]
            plot(ustrip.(line.λ[start:stop]),
                 ustrip.(I_plot[start:stop]),
                 xaxis="λ [nm]",
                 yaxis="Intensity [kW m^-2 nm^-1]",
                 dpi=300,
                 rightmargin=10Plots.mm,
                 title=title*loc)

            savefig("../img/compare_line/lines/"*title*loc)
        end
    end
end

function plot_top_cont(atmos::Atmosphere,
                       λ::Vector{<:Unitful.Length},
                       S_λ::Array{<:UnitsIntensity_λ, 4},
                       α_tot::Array{<:PerLength, 4},
                       θ::Float64,
                       ϕ::Float64,
                       title::String)

    nλ, nz, nx, ny = size(S_λ)
    I_λ = zero(S_λ)

    k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]
    for iλ in 1:nλ
        I_λ[iλ, :, :, :] = short_characteristics_up(k, S_λ[iλ, :, :, :], S_λ[iλ, 1, :, :],
                                                    α_tot[iλ, :, :, :], atmos)
    end

    for idx in 1:10:size(S_λ)[end]
        for idy in 1:10:size(S_λ)[end-1]
            loc = "_"*string(idx)*"_"*string(idy)
            I_plot = ustrip(uconvert.(u"kW*nm^-1*m^-2", I_λ))[:, end, idx, idy]
            plot(ustrip.(λ[:]),
                 ustrip.(I_plot[:]),
                 xaxis="λ [nm]",
                 yaxis="Intensity [kW m^-2 nm^-1]",
                 dpi=300,
                 rightmargin=10Plots.mm,
                 title=title*loc,
                 clim=(5,200))

            savefig("../img/compare_continuum/"*title*loc)
        end
    end
end

"""
    read_quantities(DATA::String; periodic=true)

read quantities from simulation from a hdf5 file
"""
function read_quantities(DATA::String; periodic=true)
    atmos = Atmosphere(get_atmos(DATA; periodic=periodic)...)
    global S_λ, populations
    h5open(DATA, "r") do file
        S_λ = read(file, "source_function")[:, :, :, :]*u"kW*m^-2*nm^-1"
        populations = read(file, "populations")[:, :, :, :]*u"m^-3"
    end

    if periodic
        S_λ = periodic_borders(S_λ)
        populations = periodic_pops(populations)
    end
    return atmos, S_λ, populations
end

"""
    read_irregular(DATA::String)

read quantities from irregular grid simulation from a hdf5 file
"""
function read_irregular(DATA::String)
    local positions, temperature, electron_density, hydrogen_populations, S_λ, populations
    local velocity_z, velocity_x, velocity_y, boundaries
    println("--Extracting simulation results---")
    h5open(DATA, "r") do file
        positions = read(file, "positions")[:, :]*u"m"
        boundaries = read(file, "boundaries")[:]u"m"

        temperature = read(file, "temperature")[:]*u"K"

        electron_density = read(file, "electron_density")[:]*u"m^-3"
        hydrogen_populations = read(file, "hydrogen_populations")[:]*u"m^-3"

        velocity_z = read(file, "velocity_z")[:]u"m*s^-1"
        velocity_x = read(file, "velocity_x")[:]u"m*s^-1"
        velocity_y = read(file, "velocity_y")[:]u"m*s^-1"

        S_λ = read(file, "source_function")[:, :]*u"kW*m^-2*nm^-1"
        populations = read(file, "populations")[:, :]*u"m^-3"
    end

    nn = Matrix{Int}(undef, (1, 1))
    ll = Vector{Int}(undef, 1)

    sites = VoronoiSites(positions, nn, ll, ll, ll, ll, temperature,
                         electron_density, hydrogen_populations, velocity_z,
                         velocity_x, velocity_y, boundaries...,
                         size(positions)[1])

    atmos_size = (430, 256, 256)
    atmos_size = floor.(Int, atmos_size.*1.0)

    println("--Converting grid---")
    atmos, S_λ_grid, populations_grid = Voronoi_to_Raster(sites, atmos_size,
                                                          S_λ, populations;
                                                          periodic=true)

    return atmos, S_λ_grid, populations_grid
end

"""
    plotter(atmos::Atmosphere,
            S_λ::Array{<:UnitsIntensity_λ, 4},
            populations::Array{<:NumberDensity, 4},
            θ::Float64,
            ϕ::Float64)

Compute simulations results for beam angle
"""
function plotter(atmos::Atmosphere,
                 S_λ::Array{<:UnitsIntensity_λ, 4},
                 populations::Array{<:NumberDensity, 4},
                 θ::Float64,
                 ϕ::Float64,
                 title::String)

    line = HydrogenicLine(test_atom(50, 20)..., atmos.temperature)

    LTE_pops = LTE_populations(line, atmos)

    # Find continuum extinction
    α_cont = α_absorption.(line.λ0,
                           atmos.temperature,
                           atmos.electron_density*1.0,
                           LTE_pops[:,:,:,1].+LTE_pops[:,:,:,2],
                           LTE_pops[:,:,:,3]) .+
             α_scattering.(line.λ0,
                           atmos.electron_density,
                           LTE_pops[:,:,:,1])

    γ = γ_constant(line,
                   atmos.temperature,
                   (populations[:, :, :, 1].+populations[:, :, :, 2]),
                   atmos.electron_density)

    damping_λ = Array{Float64, 4}(undef, size(S_λ))
    Threads.@threads for l in eachindex(line.λ)
        damping_λ[l, :, :, :] = damping.(γ, line.λ[l], line.ΔD)
    end

    k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]
    profile = compute_voigt_profile(line, atmos, damping_λ, k)

    α_tot = Array{Float64, 4}(undef, size(profile))u"m^-1"
    Threads.@threads for l in eachindex(line.λ)
        α_tot[l, :, :, :] = αline_λ(line,
                                    profile[l, :, :, :],
                                    populations[:, :, :, 2],
                                    populations[:, :, :, 1]) .+ α_cont
    end

    return atmos, line, S_λ, α_tot
end

function plot_convergence(DATA::String, title::String)
    local convergence
    h5open(DATA, "r") do file
        convergence = read(file, "convergence")[:]
    end


    converged = argmin(convergence)
    # println(convergence[1:converged-1])
    plot(convergence[1:converged-1],
         xlabel="iteration",
         ylabel="max rel. diff.",
         title=title,
         dpi=300,
         yscale=:log10)
    savefig("../img/$title")
    # println(title)
end

function write_convergence(DATA::String, title::String)
    local convergence
    h5open(DATA, "r") do file
        convergence = read(file, "convergence")[:]
    end

    fname = split(DATA, "/")[end]
    fname = split(fname, ".")[1]

    converged = argmin(convergence)
    convergence=convergence[1:converged-1]

    npzwrite("../python/convergence/$(title)_convergence.npy", convergence)
end
