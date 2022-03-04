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

    heatmap(ustrip(atmos.x[2:end-1]),
         ustrip(atmos.y[2:end-1]),
         transpose(I_top),
         xaxis="x",
         yaxis="y",
         dpi=300,
         rightmargin=10Plots.mm,
         title=title*" at $(round(ustrip(line.λ[idλ]); digits=2)) nm",
         aspect_ratio=:equal)

    savefig("../img/compare_line/"*title)
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
    stop =  line.λidx[2] #size(S_λ)[1]

    indices = sortperm(line.λ[start:stop])

    for idx in 1:5:size(S_λ)[end]
        for idy in 1:5:size(S_λ)[end-1]
    for idx in 1:10:size(S_λ)[end]
        for idy in 1:10:size(S_λ)[end-1]
            loc = "_"*string(idx)*"_"*string(idy)
            I_plot = ustrip(uconvert.(u"kW*nm^-1*m^-2", I_λ))[:, end, idx, idy]
            plot(ustrip.(line.λ[start+3:stop-3]),
                 ustrip.(I_plot[start+3:stop-3]),
                 xaxis="λ [nm]",
                 yaxis="Intensity [kW m^-2 nm^-1]",
                 dpi=300,
                 rightmargin=10Plots.mm,
                 title=title*loc)

            savefig("../img/compare_line/lines/"*title*loc)
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

    atmos_size = (72, 64, 64)
    atmos_size = floor.(Int, (1.5*72, 1.5*64, 1.5*64))

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

plot simulations results
"""
function plotter(atmos::Atmosphere,
                 S_λ::Array{<:UnitsIntensity_λ, 4},
                 populations::Array{<:NumberDensity, 4},
                 θ::Float64,
                 ϕ::Float64,
                 title::String)

    line = HydrogenicLine(test_atom(50, 20)..., atmos.temperature)

    γ = γ_constant(line,
                   atmos.temperature,
                   (populations[:, :, :, 1].+populations[:, :, :, 2]),
                   atmos.electron_density)

    damping_λ = Array{Float64, 4}(undef, size(S_λ))
    for l in eachindex(line.λ)
        damping_λ[l, :, :, :] = damping.(γ, line.λ[l], line.ΔD)
    end

    k = [cos(θ*π/180), cos(ϕ*π/180)*sin(θ*π/180), sin(ϕ*π/180)*sin(θ*π/180)]
    profile = compute_voigt_profile(line, atmos, damping_λ, k)

    α_line = Array{Float64, 4}(undef, size(profile))u"m^-1"
    for l in eachindex(line.λ)
        α_line[l, :, :, :] = αline_λ(line,
                                     profile[l, :, :, :],
                                     populations[:, :, :, 2],
                                     populations[:, :, :, 1])
    end

    LTE_pops = LTE_populations(line, atmos)*1.0

    # Find continuum extinction and absorption extinction (only with Thomson and Rayleigh)
    α_cont = Array{Float64, 4}(undef, (length(line.λ), size(atmos.temperature)...))u"m^-1"
    for l in eachindex(line.λ)
        α_cont[l, :, :, :] = α_continuum.(line.λ[l],
                                          atmos.temperature*1.0,
                                          atmos.electron_density*1.0,
                                          LTE_pops[:, :, :, 1]*1.0,
                                          LTE_pops[:, :, :, 3]*1.0)
    end

    α_tot = α_line .+ α_cont

    plot_top_line(atmos, line, S_λ, α_tot, θ, ϕ, title)

    # for i in 51:91
        # plot_top_intensity(atmos, line, S_λ, α_tot, θ, ϕ, i, "i_map/top_intensity_regular"*string(i))
    # end
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
         yscale=:log10,
         ylim=(0.5e-3, 1.0e1))
         yscale=:log10)
    savefig("../img/$(split(title)[1])")
end
