abstract type AbstractChannel end
abstract type ElasticChannel <: AbstractChannel end
abstract type InelasticChannel <: AbstractChannel end

struct JpsiRho <: InelasticChannel
    mJψ::Float64
    mπ::Float64
    mρ::Float64
    Γρ::Float64
end

struct JpsiOmega <: InelasticChannel
    mJψ::Float64
    mπ::Float64
    mω::Float64
    Γω::Float64
end

JpsiRho(particle_data::ParticleData) =
    JpsiRho(particle_data.mJψ, particle_data.mπ, particle_data.mρ, particle_data.Γρ)
JpsiOmega(particle_data::ParticleData) =
    JpsiOmega(particle_data.mJψ, particle_data.mπ, particle_data.mω, particle_data.Γω)

struct DxD <: ElasticChannel
    mDˣ::Float64
    mD::Float64
    DxD(mDˣ, mD) = new(Float64(mDˣ), Float64(mD))
end

struct Other <: InelasticChannel end

threshold(channel::JpsiRho) =
    channel.mJψ + 2 * channel.mπ
threshold(channel::JpsiOmega) =
    channel.mJψ + 3 * channel.mπ
threshold(channel::DxD) = channel.mDˣ + channel.mD
threshold(::Other) = -Inf

reduced_mass(channel::DxD) = channel.mD * channel.mDˣ / (channel.mD + channel.mDˣ)

k(E::Complex, channel::DxD, reference_mass) =
    1im * sqrt(-2 * reduced_mass(channel) * (E * 1e-3 - (threshold(channel) - reference_mass)))
k(E::Real, channel::DxD, reference_mass) = k(E + 1e-7im, channel, reference_mass)

"""
    FlatteModel((; Ef_MeV, g, Γ₀_MeV, fρ, fω); particle_data=ParticleData())

parametrises the X(3872) lineshape according to Eq.(7) in arXiv: 0704.0605.
The model stores the elastic and inelastic denominator contributions as channels.
When `fρ` and `fω` are supplied, the model contains five channels:
`DxD`, `DxD`, `Other`, `JpsiRho`, and `JpsiOmega`.
Without `fρ` and `fω`, the model contains the three channels `DxD`, `DxD`,
and `Other`.
"""
struct FlatteModel{C<:Tuple}
    Ef_MeV::Float64
    g::Float64
    Γ₀_MeV::Float64
    fρ::Float64
    fω::Float64
    particle_data::ParticleData
    channels::C
end

FlatteModel(; particle_data=ParticleData(), kwargs...) =
    FlatteModel((; kwargs...); particle_data)

function FlatteModel(
    pars::NamedTuple{(:Ef_MeV, :g, :Γ₀_MeV, :fρ, :fω)};
    particle_data=ParticleData(),
)
    @unpack Ef_MeV, g, Γ₀_MeV, fρ, fω = pars
    channels = (
        DxD(particle_data.mDˣ⁰, particle_data.mD⁰),
        DxD(particle_data.mDˣ⁺, particle_data.mD⁺),
        Other(),
        JpsiRho(particle_data),
        JpsiOmega(particle_data),
    )
    return FlatteModel(
        Float64(Ef_MeV), Float64(g), Float64(Γ₀_MeV), Float64(fρ), Float64(fω),
        particle_data, channels)
end

function FlatteModel(
    pars::NamedTuple{(:Ef_MeV, :g, :Γ₀_MeV)};
    particle_data=ParticleData(),
)
    @unpack Ef_MeV, g, Γ₀_MeV = pars
    channels = (
        DxD(particle_data.mDˣ⁰, particle_data.mD⁰),
        DxD(particle_data.mDˣ⁺, particle_data.mD⁺),
        Other(),
    )
    return FlatteModel(
        Float64(Ef_MeV), Float64(g), Float64(Γ₀_MeV), 0.0, 0.0,
        particle_data, channels)
end

neutral_threshold(model::FlatteModel) = threshold(model.channels[1])
threshold_offset_MeV(channel::DxD, neutral::DxD) =
    1e3 * (threshold(channel) - threshold(neutral))

contribution(model::FlatteModel, ::Type{JpsiRho}, E) =
    0.5im * model.fρ * BXρ(E, model.channels[4], neutral_threshold(model))
contribution(model::FlatteModel, ::Type{JpsiOmega}, E) =
    0.5im * model.fω * BXω(E, model.channels[5], neutral_threshold(model))
contribution_neutral(model::FlatteModel, E) =
    0.5im * model.g * k(E, model.channels[1], neutral_threshold(model))
contribution_charged(model::FlatteModel, E) =
    0.5im * model.g * k(E, model.channels[2], neutral_threshold(model))
contribution(model::FlatteModel, ::Type{Other}, E) =
    0.5im * model.Γ₀_MeV * 1e-3

"""
    shift_Ef(g, Ef_corr, particle_data)

Calculates the Ef value in MeV from the corrected energy parameter and coupling by adding dispersive contribution from the charged DxD channel to it. 

# Arguments
- `g`: Coupling parameter to the Dˣ⁰D⁰ channel
- `Ef_corr`: Energy parameter in MeV
- `particle_data::ParticleData`: Particle masses and widths

# Returns
- `Ef_MeV::Float64`: Shifted effective energy parameter in MeV
"""
function shift_Ef(g, Ef_corr, particle_data::ParticleData)
    _zero = zero(Ef_corr)
    model = FlatteModel((; Ef_MeV=_zero, g, Γ₀_MeV=_zero); particle_data)
    charged = DxD(particle_data.mDˣ⁺, particle_data.mD⁺)
    neutral = DxD(particle_data.mDˣ⁰, particle_data.mD⁰)
    Ef_GeV = threshold_offset_MeV(charged, neutral) > Ef_corr ?
             real(contribution_charged(model, Ef_corr)) :
             _zero
    Ef_MeV = 1e3 * Ef_GeV
    return Ef_MeV
end

"""
    ReparametrizeFlatte(pars; particle_data=ParticleData())

Creates a `FlatteModel` using corrected energy `Ef_corr` instead of `Ef_MeV`.
Dispatches on the named tuple keys the same way as [`FlatteModel`](@ref): five-channel
parameters include `fρ` and `fω`; the three-channel case omits them.
"""
ReparametrizeFlatte(; particle_data=ParticleData(), kwargs...) =
    ReparametrizeFlatte((; kwargs...); particle_data)

function ReparametrizeFlatte(
    pars::NamedTuple{(:Ef_corr, :g, :Γ₀_MeV, :fρ, :fω)};
    particle_data=ParticleData(),
)
    @unpack Ef_corr, g, Γ₀_MeV, fρ, fω = pars
    Ef_MeV = shift_Ef(g, Ef_corr, particle_data)
    return FlatteModel((; Ef_MeV, g, Γ₀_MeV, fρ, fω); particle_data)
end

function ReparametrizeFlatte(
    pars::NamedTuple{(:Ef_corr, :g, :Γ₀_MeV)};
    particle_data=ParticleData(),
)
    @unpack Ef_corr, g, Γ₀_MeV = pars
    Ef_MeV = shift_Ef(g, Ef_corr, particle_data)
    return FlatteModel((; Ef_MeV, g, Γ₀_MeV); particle_data)
end

"""
    compute_corrected_Ef(Ef, g, particle_data; Ef_corr_guess = -0.04)

Computes the corrected energy parameter by numerically solving the inverse relationship
between physical Ef and corrected Ef_corr parameters.

# Arguments
- `Ef::Float64`: Target physical energy parameter in MeV
- `g::Float64`: Coupling parameter to the Dˣ⁰D⁰ channel
- `particle_data::ParticleData`: Particle masses and widths used in the correction
- `Ef_corr_guess::Float64`: Initial guess for the numerical solver (default: -0.04)

# Returns
- `NamedTuple`: Contains the solver result (`sol`) and the corrected energy parameter (`Ef_corr`)
"""
function compute_corrected_Ef(Ef, g, particle_data::ParticleData; Ef_corr_guess=-0.04)
    sol = nlsolve(x -> (shift_Ef(g, x[1], particle_data) - Ef), [Ef_corr_guess])
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
    D = (E - model.Ef_MeV) * 1e-3 +
        contribution_neutral(model, E) +
        contribution_charged(model, E) +
        contribution(model, Other, E)
    if length(model.channels) == 5
        D += contribution(model, JpsiRho, E) + contribution(model, JpsiOmega, E)
    end
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
    scattering_parameters(::Type{FlatteModel}, Ef_MeV, g, particle_data)

Calculates the scattering parameters of the Flatte model according to arXiv: 2108.11413.
Returns the inverse scattering length and effective range.

# Arguments
- `::Type{FlatteModel}`: The FlatteModel type
- `Ef_MeV::Float64`: Effective energy parameter in MeV
- `g::Float64`: Coupling parameter
- `particle_data::ParticleData`: Particle masses and widths

# Returns
- `NamedTuple`: Contains inverse scattering length (`inva`) and effective range (`r`)
"""
function scattering_parameters(::Type{FlatteModel}, Ef_MeV, g, particle_data::ParticleData)
    # expressions from arXiv: 2108.11413
    μ = reduced_mass_neutral(particle_data)
    μ⁺ = reduced_mass_charged(particle_data)
    δ⁺ = charged_threshold_offset(particle_data)
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
    scattering_parameters(typeof(model), model.Ef_MeV, model.g, model.particle_data)

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
function pole_position(model::FlatteModel, init=-1e3im * model.Γ₀_MeV / 10)
    fr = optimize(x -> abs2(denominator(model, x[1] + x[2] * 1im)), collect(reim(init)), BFGS())
    minimum_reached = (fr.minimum < 1e-8)
    !(minimum_reached) && error("Pole is not found: fr.minimum = $(fr.minimum)")
    Epole = complex(fr.minimizer...) # MeV
    return Epole
end
