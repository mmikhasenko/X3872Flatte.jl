
abstract type Model end

"""
    FlatteModel(Ef, g, Γ₀)

parametrises the X(3872) lineshape according to Eq.(7) in arXiv: 0704.0605.
Parameters contained in the structure are
 - Ef_MeV: parameter of the mass peak
 - g : coupling to the Dˣ⁰D⁰
 - Γ₀_MeV: contribution to the width from other inelastic channels
"""
@with_kw struct FlatteModel <: Model
    Ef_MeV::Float64
    g::Float64
    Γ₀_MeV::Float64
    fρ::Float64
    fω::Float64
end


"""
Demominator of the X amplitude. Eq.(7) in arXiv: 0704.0605.
"""
function denominator(model::FlatteModel, E) # E is in MeV
    @unpack Ef_MeV, g, Γ₀_MeV = model
    @unpack fρ, fω = model
    #
    Bρ = BXρ(E)
    Bω = BXω(E)
    # 
    D = (E - Ef_MeV) * 1e-3 + 0.5im * (g * k1(E) + g * k2(E)) +
        0.5im * (Γ₀_MeV * 1e-3 + fρ * Bρ + fω * Bω)
    return D
end



"""
    AJψππ(model::FlatteModel, E)

Transition amplitude of X → J/ψ π π where the decay constant in the numerator is omitted.
The functional dependence is the same as for the Dˣ⁰ D̄⁰ → Dˣ⁰ D̄⁰.
"""
AJψππ(model::FlatteModel, E) = 1 / denominator(model::FlatteModel, E)

"""
    scattering_parameters(model::FlatteModel)

Scattering parameters of the Flatte model according to arXiv: 2108.11413
"""
function scattering_parameters(::Type{FlatteModel}, Ef_MeV, g)
    # expressions from arXiv: 2108.11413
    inva_GeV = (2 * Ef_MeV * 1e-3) / g + sqrt(2 * μ⁺ * δ⁺)  # Eq.18a 
    inva = inva_GeV * 1e3
    # 
    r_GeV⁻¹ = -2 / (μ * g) - sqrt(μ⁺ / (2 * μ^2 * δ⁺)) # Eq.18b  
    r = r_GeV⁻¹ * fm_times_GeV
    return (; inva, r)
end

scattering_parameters(model::FlatteModel) =
    scattering_parameters(typeof(model), model.Ef_MeV, model.g)

"""
    pole_position(model::FlatteModel)

Position of the pole singularity of the Flatte amplitude.
The gradient decent (BFGS) is started from [0, -model.Γ₀_MeV/2] point.
"""
function pole_position(model::FlatteModel, init = -1e3im * model.Γ₀_MeV / 10)
    fr = optimize(x -> abs2(denominator(model, x[1] + x[2] * 1im)), collect(reim(init)), BFGS())
    minimum_reached = (fr.minimum < 1e-8)
    !(minimum_reached) && error("Pole is not found: fr.minimum = $(fr.minimum)")
    Epole = complex(fr.minimizer...) # MeV
    return Epole
end
