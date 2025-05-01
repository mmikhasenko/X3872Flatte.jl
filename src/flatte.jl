"""
    FlatteModel(Ef, g, Γ₀)

parametrises the X(3872) lineshape according to Eq.(7) in arXiv: 0704.0605.
Parameters contained in the structure are
 - Ef_MeV: parameter of the mass peak
 - g : coupling to the Dˣ⁰D⁰
 - Γ₀_MeV: contribution to the width from other inelastic channels
"""
@with_kw struct FlatteModel
    Ef_MeV::Float64
    g::Float64
    Γ₀_MeV::Float64
    fρ::Float64
    fω::Float64
end

"""
    shift_Ef(g, Ef_corr)

Calculates the Ef value in MeV from the corrected energy parameter and coupling by adding dispersive contribution from the charged DxD channel to it. 

# Arguments
- `g`: Coupling parameter to the Dˣ⁰D⁰ channel
- `Ef_corr`: Energy parameter in MeV

# Returns
- `Ef_MeV::Float64`: Shifted effective energy parameter in MeV
"""
function shift_Ef(g, Ef_corr)
    _model = FlatteModel(; Ef_MeV = 0.0, g, Γ₀_MeV = 0.0, fρ = 0.0, fω = 0.0)
    Ef_GeV = denominator(_model, Ef_corr) |> real
    Ef_MeV = 1e3 * Ef_GeV
    return Ef_MeV
end

"""
    ReparametrizeFlatte(pars_corr)

Creates a FlatteModel instance using corrected Ef parameter instead of Ef_MeV along with the other parameters.

# Arguments
- `pars_corr::NamedTuple`: Contains `Ef_corr`, `g`, `Γ₀_MeV`, `fρ`, and `fω` parameters

# Returns
- `FlatteModel`: Model with physical parameters
"""
function ReparametrizeFlatte(pars_corr)
    @unpack Ef_corr, g, Γ₀_MeV, fρ, fω = pars_corr
    Ef_MeV = shift_Ef(g, Ef_corr)
    return FlatteModel(; Ef_MeV, g, Γ₀_MeV, fρ, fω)
end

"""
    compute_corrected_Ef(Ef, g, Ef_corr_guess = -0.04)

Computes the corrected energy parameter by numerically solving the inverse relationship
between physical Ef and corrected Ef_corr parameters.

# Arguments
- `Ef::Float64`: Target physical energy parameter in MeV
- `g::Float64`: Coupling parameter to the Dˣ⁰D⁰ channel
- `Ef_corr_guess::Float64`: Initial guess for the numerical solver (default: -0.04)

# Returns
- `NamedTuple`: Contains the solver result (`sol`) and the corrected energy parameter (`Ef_corr`)
"""
function compute_corrected_Ef(Ef, g, Ef_corr_guess = -0.04)
    sol = nlsolve(x -> (shift_Ef(g, x[1]) - Ef), [Ef_corr_guess])
    Ef_corr = sol.zero[1]
    (; sol, Ef_corr)
end

"""
    denominator(model::FlatteModel, E)

Calculates the denominator of the X(3872) amplitude according to Eq.(7) in arXiv: 0704.0605.

# Arguments
- `model::FlatteModel`: The Flatte model parameters
- `E::Float64`: Energy in MeV

# Returns
- Complex denominator value of the amplitude
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

Calculates the transition amplitude of X → J/ψ π π where the decay constant in the numerator is omitted.
The functional dependence is the same as for the Dˣ⁰ D̄⁰ → Dˣ⁰ D̄⁰.

# Arguments
- `model::FlatteModel`: The Flatte model parameters
- `E::Float64`: Energy in MeV

# Returns
- Complex amplitude value
"""
AJψππ(model::FlatteModel, E) = 1 / denominator(model::FlatteModel, E)

"""
    scattering_parameters(::Type{FlatteModel}, Ef_MeV, g)

Calculates the scattering parameters of the Flatte model according to arXiv: 2108.11413.
Returns the inverse scattering length and effective range.

# Arguments
- `::Type{FlatteModel}`: The FlatteModel type
- `Ef_MeV::Float64`: Effective energy parameter in MeV
- `g::Float64`: Coupling parameter

# Returns
- `NamedTuple`: Contains inverse scattering length (`inva`) and effective range (`r`)
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

"""
    scattering_parameters(model::FlatteModel)

Calculates the scattering parameters (inverse scattering length and effective range)
using the parameters from the provided FlatteModel instance.

# Arguments
- `model::FlatteModel`: The Flatte model parameters

# Returns
- `NamedTuple`: Contains inverse scattering length (`inva`) and effective range (`r`)
"""
scattering_parameters(model::FlatteModel) =
    scattering_parameters(typeof(model), model.Ef_MeV, model.g)

"""
    pole_position(model::FlatteModel, init = -1e3im * model.Γ₀_MeV / 10)

Finds the position of the pole singularity of the Flatte amplitude using numerical optimization.
The gradient descent method (BFGS) is used to locate the complex energy where the denominator
approaches zero.

# Arguments
- `model::FlatteModel`: The Flatte model parameters
- `init::Complex`: Initial guess for the optimization (default: scaled by model.Γ₀_MeV)

# Returns
- `Epole::Complex`: Complex energy position of the pole in MeV
"""
function pole_position(model::FlatteModel, init = -1e3im * model.Γ₀_MeV / 10)
    fr = optimize(x -> abs2(denominator(model, x[1] + x[2] * 1im)), collect(reim(init)), BFGS())
    minimum_reached = (fr.minimum < 1e-8)
    !(minimum_reached) && error("Pole is not found: fr.minimum = $(fr.minimum)")
    Epole = complex(fr.minimizer...) # MeV
    return Epole
end
