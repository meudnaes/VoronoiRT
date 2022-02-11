include("functions.jl")

"""
    c4_traving(line::AtomicLine)
Calculate the \$C_4\$ interaction constant (quadratic Stark effect) using the
recipe of Traving (1960), "Uber die Theorie der Druckverbreiterung von Spektrallinien",
p 93.
"""
function c4_traving(line::HydrogenicLine)
    n_eff_u = Transparency.n_eff(line.χ∞, line.χj, line.Z)
    n_eff_l = Transparency.n_eff(line.χ∞, line.χi, line.Z)
    C4 = (Transparency.e^2 * Transparency.inv_4πε0 * Transparency.a_0^3 * 2 * π / (h * 18 * line.Z^4) *
        ((n_eff_u * (5 * n_eff_u^2 + 1))^2 - (n_eff_l * (5 * n_eff_l^2 + 1))^2))
    return C4 |> u"m^4 / s"
end

"""
    function const_unsold(line::HydrogenicLine; H_scaling::Real=1, He_scaling::Real=1)
Compute atmosphere-independent constant for γ_unsold, to be used in function `γ_unsold`.
Based on expressions from RH broad.c, which uses formula in Mihalas (1978),
pp 282, 286-287, eq. (9-50) for v_rel, table 9-1 and eq. (9-76) for the interaction
coefficient C6.  Input is an `AtomicLine`, which contains the necessary energies and
element weight. The van der Waals broadening can be scaled for both H and He perturbers
using `H_scaling` and `He_scaling`.
"""
function const_unsold(line::HydrogenicLine; H_scaling::Real=1, He_scaling::Real=1)
    Δr = (Transparency.Ry^2 * (1 / (line.χ∞ - line.χj)^2 - 1 / (line.χ∞ - line.χi)^2)) |> u"J/J"
    C6 = ustrip((2.5 * Transparency.e^2 * Transparency.αp * Transparency.inv_4πε0^2 * 2 * π *
                 (line.Z * Transparency.a_0)^2 / h * Δr) |> u"C^2 * m^6 / (F * J * s)")
    v_rel_const = ustrip(8 * k_B / (π * line.atom_weight) |> u"J/(K * kg)")
    v_rel_H = v_rel_const * (1 + line.atom_weight / Transparency.mass_H)
    v_rel_He = v_rel_const * (1 + line.atom_weight / Transparency.mass_He)
    return 8.08 * (H_scaling * v_rel_H^0.3 + He_scaling * Transparency.abund_He * v_rel_He^0.3) * C6^0.4
end

"""
    const_quadratic_stark(line::AtomicLine;
                          mean_atomic_weight::Unitful.Mass=28 * m_u,
                          scaling::Real=1)
Calculate height-independent constant to use in `γ_quadratic_stark`, using the recipe
from RH, which is based on the following estimate:
\$\$
\\gamma = 11.37 \\cdot vrel^{1/3} * C_4^{2/3} * (n_e + n_{ion}),
\$\$
Using the estimate for \$C_4\$ from Traving (1960), "Uber die Theorie der
Druckverbreiterung von Spektrallinien", p 93., and \$n_{ion}\\approx n_e\$
(following Gray).
"""
function const_quadratic_stark(line::HydrogenicLine;
                               mean_atomic_weight::Unitful.Mass=28 * Transparency.m_u,
                               scaling::Real=1)
    C = ustrip(8 * k_B / (π * line.atom_weight) |> u"J/(K * kg)")
    Cm = ((1 + line.atom_weight / m_e)^(1/6) +
          (1 + line.atom_weight / mean_atomic_weight)^(1/6))
    C4 = ustrip(c4_traving(line) |> u"m^4 / s")
    cStark23 = 11.37u"m^3 / s" * (scaling * C4)^(2/3)
    return C^(1/6) * cStark23 * Cm
end

function γ_constant(line::HydrogenicLine,
                    temperature::Array{<:Unitful.Temperature},
                    neutral_hydrogen_density::Array{<:NumberDensity},
                    electron_density::Array{<:NumberDensity})

    u=2 #line.j
    l=1 #line.i

    unsold_const = const_unsold(line)
    quad_stark_const = const_quadratic_stark(line)

    γ = γ_unsold.(unsold_const, temperature, neutral_hydrogen_density)
    # γ .+= line.Aji
    γ .+= 4.702e8u"s^-1" # check
    γ .+= γ_linear_stark.(electron_density, u, l)
    γ .+= γ_quadratic_stark.(electron_density, temperature, stark_constant=quad_stark_const)

    return γ

end
