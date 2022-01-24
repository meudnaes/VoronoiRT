struct Atom
    density::Array{<:NumberDensity, 3}                  # (nz, nx, ny)
    n_levels::Integer
    n_lines::Integer
    χ::Vector{<:Unitful.Energy}
    g::Vector{Int64}
    Z::Integer
    f_value::Vector{Float64}
    λ::Vector{Unitful.Length}
    nλ::Integer
    iλbb
    iλbf
end

struct Line
    u::Int
    l::Int
    lineData::AtomicLine
    doppler_width::Array{<:Unitful.Length, 3}        # (n_lines, nz, nx, ny)
    damping_constant::Array{<:PerArea, 3}            # (n_lines, nx, ny)
end

"""
    sample_λ_line(nλ::Int64, χl::Unitful.Energy, χu::Unitful.Energy,
                            qwing::Float64, qcore::Float64)


Get sampling wavelengths. Bound free wavelengths are
linearly sampled, while the bound-bound follow the
log-sampling from github.com/ITA-Solar/rh.
"""
function sample_λ_line(nλ::Int64,
                       λ0::Unitful.Length,
                       qwing::Float64,
                       qcore::Float64)

    # Make sure odd # of bb wavelengths
    if nλ > 0 && nλ%2 == 0
        nλ += 1
    end

    # Either 1 or five or more wavelengths
    if 1 < nλ < 5
        nλ = 5
    end

    # Initialise wavelength array
    λ = Array{Float64,1}(undef, nλ)u"nm"

    # =================================================
    # Bound-bound transition
    # Follows github.com/ITA-Solar/rh/blob/master/getlambda.c
    # =================================================
    if nλ == 1
        λ[1] = λ0

    elseif nλ >= 5
        vmicro_char = 2.5u"km/s"

        n = nλ/2 # Questionable
        β = qwing/(2*qcore)
        y = β + sqrt(β*β + (β - 1.0)*n + 2.0 - 3.0*β)
        b = 2.0*log(y) / (n - 1)
        a = qwing / (n - 2.0 + y*y)

        center = (nλ÷2) + 1
        λ[center] = λ0
        q_to_λ = λ[center] * vmicro_char / c_0

        for w=1:(nλ÷2)
            Δλ = a*(w + (exp(b*w) - 1.0)) * q_to_λ
            λ[center-w] = λ[center] - Δλ
            λ[center+w] = λ[center] + Δλ
        end
    end

    return λ
end

"""
    sample_λ(nλ_bb::Int64, nλ_bf::Int64,
             χl::Unitful.Energy, χu::Unitful.Energy, χ∞::Unitful.Energy)

Get sampling wavelengths. Bound free wavelengths are
linearly sampled, while the bound-bound follow the
log-sampling from github.com/ITA-Solar/rh.
"""
function sample_λ_boundfree(nλ::Int64,
                            λ_min::Unitful.Length,
                            χl::Unitful.Energy,
                            χ∞::Unitful.Energy)


    λ_max  = transition_λ(χl, χ∞)

    # Initialise wavelength array
    λ = Array{Float64,1}(undef, nλ)u"nm"

    # =================================================
    # Bound-free transitions
    # Linear spacing
    # =================================================
    if nλ == 1

        λ[1] = λ_max

    elseif nλ > 1
        Δλ = (λ_max - λ_min)/(nλ-1)
        λ[1] = λ_min

        for w=2:nλ
            λ[w] = λ[w-1] + Δλ
        end
    end

    return λ
end


"""
    transition_λ(χ1::Unitful.Energy, χ2::Unitful.Energy)

Get the corresponding wavelength for
the energy difference between two levels.
"""
function transition_λ(χ1::Unitful.Energy, χ2::Unitful.Energy)
    ((h * c_0) / (χ2-χ1)) |> u"nm"
end


"""
    damping_constant(γ::Unitful.Frequency,
                     ΔλD::Unitful.Length)

Get daping constant to be multiplied with λ^2.
"""
function damping_constant(γ::Unitful.Frequency,
                          ΔλD::Unitful.Length)
    (γ / (4 * π * c_0 * ΔλD))
end



# ==================================================================
#  WRITE TO FILE
# ==================================================================


"""
    write_to_file(atom::Atom, output_path::String)

Writes wavelengths to the output file.
"""
function write_to_file(λ::Array{<:Unitful.Length,1}, output_path::String)
    h5open(output_path, "r+") do file
        write(file, "wavelength", ustrip(λ))
    end
end


"""
    write_to_file(atom::Atom, output_path::String)

Writes wavelength and number of bound-bound and
bound-free wavelengths to the output file.
"""
function write_to_file(λ, iλbf, iλbb, output_path::String)
    h5open(output_path, "r+") do file
        write(file, "wavelength", ustrip.(λ))
        #write(file, "iwbf", Arr
        #write(file, "iwbf", iλbb)
    end
end
