using Transparency

"""
    function B_ν(ν, T)

Planck's law! Radiation in LTE. Takes frequency and temperature, returns
specific intensity
"""
function B_ν(ν, T)::AbstractFloat
    return 2*h*ν^3/c_0^2 * 1/(exp(h*ν/(k_B*T)) - 1)
end

"""
    function B_ν(λ, T)

Planck's law! Radiation in LTE. Takes wavelength and temperature, returns
specific intensity
"""
function B_λ(λ::Unitful.Length, T::Unitful.Temperature)
    return 2*h*c_0^2/λ^5 * 1/(exp(h*c_0/(λ*k_B*T)) - 1)
end

"""
    function α_cont(λ::Unitful.Length, temperature::Unitful.Temperature,
               electron_density::NumberDensity, h_ground_density::NumberDensity,
               proton_density::NumberDensity)

Total continuum extinction.
"""
function α_continuum(λ::Unitful.Length, temperature::Unitful.Temperature,
                     electron_density::NumberDensity, h_ground_density::NumberDensity,
                     proton_density::NumberDensity)

    α = max(0u"m^-1", Transparency.hminus_ff_stilley(λ, temperature, h_ground_density, electron_density))
    α += Transparency.hminus_bf_wbr(λ, temperature, h_ground_density, electron_density)
    α += hydrogenic_ff(c_0 / λ, temperature, electron_density, proton_density, 1)
    α += h2plus_ff(λ, temperature, h_ground_density, proton_density)
    α += h2plus_bf(λ, temperature, h_ground_density, proton_density)
    α += thomson(electron_density)
    α += rayleigh_h(λ, h_ground_density)
    return α
end

"""
    α_scattering(λ::Unitful.Length,
                 electron_density::NumberDensity,
                 h_ground_density::NumberDensity)

Extinction from scattering processes.
"""
function α_scattering(λ::Unitful.Length,
                      electron_density::NumberDensity,
                      h_ground_density::NumberDensity)

   α = thomson(electron_density)
   α += rayleigh_h(λ, h_ground_density)
   return α
end

"""
    function α_absorption(λ::Unitful.Length, temperature::Unitful.Temperature,
               electron_density::NumberDensity, h_ground_density::NumberDensity,
               proton_density::NumberDensity)

Extinction from photon destruction processes.
"""
function α_absorption(λ::Unitful.Length, temperature::Unitful.Temperature,
                      electron_density::NumberDensity, h_ground_density::NumberDensity,
                      proton_density::NumberDensity)

    α = max(0u"m^-1", Transparency.hminus_ff_stilley(λ, temperature, h_ground_density, electron_density))
    α += Transparency.hminus_bf_wbr(λ, temperature, h_ground_density, electron_density)
    α += hydrogenic_ff(c_0 / λ, temperature, electron_density, proton_density, 1)
    α += h2plus_ff(λ, temperature, h_ground_density, proton_density)
    α += h2plus_bf(λ, temperature, h_ground_density, proton_density)
    return α
end
