using Plots

include("line.jl")
include("functions.jl")
include("atmosphere.jl")
include("broadening.jl")
include("populations.jl")
include("characteristics.jl")

function plot_top_intensity(atmos::Atmosphere,
                            S_λ::Array{<:UnitsIntensity_λ, 4},
                            α_tot::Array{<:PerLength, 4},
                            populations::Array{<:NumberDensity, 4},
                            θ::Float64,
                            ϕ::Float64,
                            idλ::Int,
                            title::String)


    I_top = short_characteristics_up(θ, ϕ, S_λ[idλ, :, :, :], α_tot[idλ, :, :, :],
                                     atmos, degrees=true, I_0=S_λ[idλ, 1, :, :])

    I_top = ustrip(uconvert.(u"kW*nm^-1*m^-2", I_top[end, 2:end-1, 2:end-1]))

    heatmap(ustrip(atmos.x[2:end-1]),
         ustrip(atmos.y[2:end-1]),
         transpose(I_top),
         xaxis="x",
         yaxis="y",
         dpi=300,
         rightmargin=10Plots.mm,
         title=title,
         aspect_ratio=:equal)

    savefig("../img/compare_line/"*title)
end

function plot_top_line(atmos::Atmosphere,
                       line::HydrogenicLine,
                       S_λ::Array{<:UnitsIntensity_λ, 4},
                       α_tot::Array{<:PerLength, 4},
                       populations::Array{<:NumberDensity, 4},
                       θ::Float64,
                       ϕ::Float64,
                       idx::Int,
                       idy::Int,
                       title::String)

    nλ, nz, nx, ny = size(S_λ)

    I_λ = Vector{Float64}(undef, nλ)u"kW*m^-2*nm^-1"


    for iλ in 1:nλ
        I_λ[iλ] = short_characteristics_up(θ, ϕ, S_λ[iλ, :, :, :], α_tot[iλ, :, :, :],
                                           atmos, degrees=true, I_0=S_λ[iλ, 1, :, :])[end, idx, idy]
    end

    I_λ = ustrip(uconvert.(u"kW*nm^-1*m^-2", I_λ))

    scatter(ustrip.(line.λ),
         ustrip.(I_λ),
         xaxis="λ [nm]",
         yaxis="Iλ [kW m^-2 nm^-1]",
         dpi=300,
         rightmargin=10Plots.mm,
         title=title)

    savefig("../img/compare_line/"*title)
end

"""
    read_quantities(DATA::String)

read quantities from simulation from a hdf5 file
"""
function read_quantities(DATA::String)
    local S_λ, populations, atmos
    h5open(DATA, "r") do file
        z = read(file, "z")[:]*u"m"
        x = read(file, "x")[:]*u"m"
        y = read(file, "y")[:]*u"m"

        velocity_z = read(file, "velocity_z")[:, :, :]*u"m/s"
        velocity_x = read(file, "velocity_x")[:, :, :]*u"m/s"
        velocity_y = read(file, "velocity_y")[:, :, :]*u"m/s"

        temperature = read(file, "temperature")[:, :, :]*u"K"
        electron_density = read(file, "electron_density")[:, :, :]*u"m^-3"
        hydrogen_density = read(file, "hydrogen_density")[:, :, :]*u"m^-3"

        atmos = Atmosphere(z, x, y, temperature, electron_density, hydrogen_density, velocity_z, velocity_x, velocity_y)

        S_λ = read(file, "source_function")[:, :, :, :]*u"kW*m^-2*nm^-1"
        populations = read(file, "populations")[:, :, :, :]*u"m^-3"
    end
    return atmos, S_λ, populations
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
                 ϕ::Float64)

    title="Line at Top"

    idx = 15
    idy = 20

    line = HydrogenicLine(test_atom()..., atmos.temperature)
    println(line.λ)
    γ = γ_constant(line,
                   atmos.temperature,
                   (populations[:, :, :, 1].+populations[:, :, :, 2]),
                   atmos.electron_density)

    damping_λ = Array{Float64, 4}(undef, size(S_λ))
    for l in eachindex(line.λ)
        damping_λ[l, :, :, :] = damping.(γ, line.λ[l], line.ΔD)
    end

    profile = compute_voigt_profile(line, atmos, damping_λ, θ*π/180, ϕ*π/180)

    α_line = Array{Float64, 4}(undef, size(profile))u"m^-1"
    for l in eachindex(line.λ)
        α_line[l, :, :, :] = αline_λ(line,
                                     profile[l, :, :, :],
                                     populations[:, :, :, 1],
                                     populations[:, :, :, 2])
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

    plot_top_line(atmos, line, S_λ, α_tot, populations, θ, ϕ, idx, idy, title)
end

plotter(read_quantities("../data/voronoi_line.h5")..., 10., 10.)
