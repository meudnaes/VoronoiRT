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
                       title::String)

    nλ, nz, nx, ny = size(S_λ)
    I_λ = zero(S_λ)


    for iλ in 1:nλ
        I_λ[iλ, :, :, :] = short_characteristics_up(θ, ϕ, S_λ[iλ, :, :, :], α_tot[iλ, :, :, :],
                                                    atmos, degrees=true, I_0=S_λ[iλ, 1, :, :])
    end

    start = line.λidx[1]+5
    stop = line.λidx[2]-4

    for idx in 1:5:50
        for idy in 1:5:50
            loc = "_"*string(idx)*"_"*string(idy)
            I_plot = ustrip(uconvert.(u"kW*nm^-1*m^-2", I_λ))[:, end, idx, idy]
            scatter(ustrip.(line.λ)[start:stop],
                    ustrip.(I_plot)[start:stop],
                    xaxis="λ [nm]",
                    yaxis="Iλ [kW m^-2 nm^-1]",
                    dpi=300,
                    rightmargin=10Plots.mm,
                    title=title*loc)

            savefig("../img/compare_line/lines/"*title*loc)
        end
    end
end

"""
    read_quantities(DATA::String)

read quantities from simulation from a hdf5 file
"""
function read_quantities(DATA::String; periodic=true)
    atmos = Atmosphere(get_atmos(DATA; periodic=periodic, skip=1)...)
    local S_λ, populations
    h5open(DATA, "r") do file
        S_λ = read(file, "source_function")[:, :, :, :]*u"kW*m^-2*nm^-1"
        populations = read(file, "populations")[:, :, :, :]*u"m^-3"
    end

    if periodic
        S_λ_periodic = Array{Float64, 4}(undef, size(S_λ) .+ (0,0,2,2))u"kW*m^-2*nm^-1"
        populations_periodic = Array{Float64, 4}(undef, size(populations) .+ (0,2,2,0))u"m^-3"

        # fill inner box
        S_λ_periodic[:,:,2:end-1, 2:end-1] = S_λ
        # x-direction wall
        S_λ_periodic[:,:,1,2:end-1] = S_λ[:,:,end,:]
        S_λ_periodic[:,:,end,2:end-1] = S_λ[:,:,1,:]
        # y-direction wall
        S_λ_periodic[:,:,2:end-1,end] = S_λ[:,:,:,1]
        S_λ_periodic[:,:,2:end-1,1] = S_λ[:,:,:,end]
        # fix corners
        S_λ_periodic[:,:,1,1] .= S_λ[:,:,end,end]
        S_λ_periodic[:,:,1,end] .= S_λ[:,:,end,1]
        S_λ_periodic[:,:,end,1] .= S_λ[:,:,1,end]
        S_λ_periodic[:,:,end,end] .= S_λ[:,:,1,1]

        # fill inner box
        populations_periodic[:,2:end-1,2:end-1,:] = populations
        # x-direction wall
        populations_periodic[:,1,2:end-1,:] = populations[:,end,:,:]
        populations_periodic[:,end,2:end-1,:] = populations[:,1,:,:]
        # y-direction wall
        populations_periodic[:,2:end-1,end,:] = populations[:,:,1,:]
        populations_periodic[:,2:end-1,1,:] = populations[:,:,end,:]
        # fix corners
        populations_periodic[:,1,1,:] .= populations[:,end,end,:]
        populations_periodic[:,1,end,:] .= populations[:,end,1,:]
        populations_periodic[:,end,1,:] .= populations[:,1,end,:]
        populations_periodic[:,end,end,:] .= populations[:,1,1,:]

        return atmos, S_λ_periodic, populations_periodic
    else
        return atmos, S_λ, populations
    end
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

    idx = 12
    idy = 12

    line = HydrogenicLine(test_atom(50, 20)..., atmos.temperature)

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

    plot_top_line(atmos, line, S_λ, α_tot, populations, θ, ϕ, title)

    # for i in 1:51
        # plot_top_intensity(atmos, S_λ, α_tot, populations, θ, ϕ, i, "i_map/top_intensity_regular"*string(i))
    # end
end

# plotter(read_quantities("../data/linedata/voronoi_line.h5", periodic=true)..., 10., 10., "Voronoi-Line")
plotter(read_quantities("../data/linedata/regular_line_91.h5", periodic=true)..., 10., 10., "Regular-Line")
