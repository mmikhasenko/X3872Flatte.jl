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

struct DxD{T<:Number} <: ElasticChannel
    mDˣ::T
    mD::T
end

function DxD(mDˣ::Number, mD::Number)
    _mDˣ, _mD = promote(mDˣ, mD)
    return DxD{typeof(_mDˣ)}(_mDˣ, _mD)
end

struct Other <: InelasticChannel end

threshold(channel::JpsiRho) =
    channel.mJψ + 2 * channel.mπ
threshold(channel::JpsiOmega) =
    channel.mJψ + 3 * channel.mπ
threshold(channel::DxD) = channel.mDˣ + channel.mD
threshold(::Other) = -Inf

reduced_mass(channel::DxD) = channel.mD * channel.mDˣ / (channel.mD + channel.mDˣ)

"""
    MomentumSheet(label)
    MomentumSheet(signs)

Riemann sheet for the two elastic Dˣ⁰D⁰ and Dˣ⁺D⁺ momenta. The first sign is
the neutral-channel momentum sign and the second sign is the charged-channel
momentum sign. Sheet labels follow `(I, II, III, IV) = ((+,+), (-,+), (-,-), (+,-))`.
"""
struct MomentumSheet{N}
    signs::NTuple{N,Int}
end

function MomentumSheet(signs::NTuple{N,<:Integer}) where {N}
    return MomentumSheet{N}(Tuple(sheet_sign.(signs)))
end

MomentumSheet(label::Symbol) = MomentumSheet(sheet_signs(Val(label)))

sheet_sign(sign::Integer) =
    sign in (-1, 1) ? Int(sign) : throw(ArgumentError("sheet signs must be -1 or +1"))

sheet_signs(::Val{:I}) = (+1, +1)
sheet_signs(::Val{:II}) = (-1, +1)
sheet_signs(::Val{:III}) = (-1, -1)
sheet_signs(::Val{:IV}) = (+1, -1)

momentum_argument(E, channel::DxD, reference_mass) =
    2 * reduced_mass(channel) * (E * 1e-3 - (threshold(channel) - reference_mass))

k(E::Complex, channel::DxD, reference_mass) =
    1im * sqrt(-momentum_argument(E, channel, reference_mass))
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
struct FlatteModel{C<:Tuple,P<:ParticleData}
    Ef_MeV::Float64
    g::Float64
    Γ₀_MeV::Float64
    fρ::Float64
    fω::Float64
    particle_data::P
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

function contribution_elastic(
    model::FlatteModel{<:Tuple{<:DxD,<:DxD,Other}},
    E,
    sheet::MomentumSheet{2},
)
    channels = (model.channels[1], model.channels[2])
    reference_mass = neutral_threshold(model)
    return 0.5im * model.g *
           sum(sign * k(E, channel, reference_mass) for (sign, channel) in zip(sheet.signs, channels))
end

"""
    shift_Ef(g, Ef_corr, particle_data)

Maps reparametrized energy `Ef_corr` to `Ef_MeV` via the calibration denominator
(`Ef = 0`, `Γ₀ = 0`, elastic DˣD channels only), evaluated at `E = Ef_corr`:

```math
E_{f,\\mathrm{MeV}} = E_{f,\\mathrm{corr}} + 10^3\\,\\mathrm{Re}\\big[\\Sigma(E_{f,\\mathrm{corr}})\\big],
```

with ``\\Sigma`` the sum of neutral and charged DˣD loop contributions.

# Arguments
- `g`: Coupling parameter to the Dˣ⁰D⁰ channel
- `Ef_corr`: Energy parameter in MeV
- `particle_data::ParticleData`: Particle masses and widths

# Returns
- `Ef_MeV::Float64`: Effective energy parameter in MeV
"""
function shift_Ef(g, Ef_corr, particle_data::ParticleData)
    model = FlatteModel((; Ef_MeV=0.0, g, Γ₀_MeV=0.0); particle_data)
    return 1e3 * real(denominator(model, Ef_corr))
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

function denominator(
    model::FlatteModel{<:Tuple{<:DxD,<:DxD,Other}},
    E,
    sheet::MomentumSheet{2},
)
    return (E - model.Ef_MeV) * 1e-3 +
           contribution_elastic(model, E, sheet) +
           contribution(model, Other, E)
end

AJψππ(model::FlatteModel{<:Tuple{<:DxD,<:DxD,Other}}, E, sheet::MomentumSheet{2}) =
    1 / denominator(model, E, sheet)

"""
    pole_parameters(pole_MeV, g; particle_data=ParticleData(), sheet=MomentumSheet(:II))

Return the standard three-channel Flatte parameters `(Ef_MeV, Γ₀_MeV)` whose
sheet-aware denominator vanishes at `pole_MeV`.

The pole energy is in MeV relative to the neutral Dˣ⁰D⁰ threshold. The returned
`Γ₀_MeV` follows the convention `Γ = -2 * imag(E_p)` when elastic loop
contributions are absent.
"""
function pole_parameters(
    pole_MeV::Number,
    g;
    particle_data=ParticleData(),
    sheet::MomentumSheet{2}=MomentumSheet(:II),
)
    calibration_model = FlatteModel((; Ef_MeV=0.0, g, Γ₀_MeV=0.0); particle_data)
    z_MeV = pole_MeV + 1e3 * contribution_elastic(calibration_model, pole_MeV, sheet)
    return (; Ef_MeV=real(z_MeV), Γ₀_MeV=-2 * imag(z_MeV))
end

"""
    PoleReparametrizeFlatte(pole_MeV, g; particle_data=ParticleData(), sheet=MomentumSheet(:II))
    PoleReparametrizeFlatte(; pole_MeV, g, kwargs...)
    PoleReparametrizeFlatte(; pole_re_MeV, pole_im_MeV, g, kwargs...)

Build the three-channel `FlatteModel` from a target pole position and elastic
coupling. The returned model has real `Ef_MeV` and `Γ₀_MeV` chosen so that
`denominator(model, pole_MeV, sheet) == 0`, up to floating-point precision.
"""
function PoleReparametrizeFlatte(; particle_data=ParticleData(), sheet=MomentumSheet(:II), kwargs...)
    pars = (; kwargs...)
    haskey(pars, :pole_MeV) &&
        return PoleReparametrizeFlatte((; pole_MeV=pars.pole_MeV, g=pars.g); particle_data, sheet)
    return PoleReparametrizeFlatte((;
        pole_re_MeV=pars.pole_re_MeV,
        pole_im_MeV=pars.pole_im_MeV,
        g=pars.g,
    ); particle_data, sheet)
end

PoleReparametrizeFlatte(
    pole_MeV::Number,
    g;
    particle_data=ParticleData(),
    sheet::MomentumSheet{2}=MomentumSheet(:II),
) =
    PoleReparametrizeFlatte((; pole_MeV, g); particle_data, sheet)

function PoleReparametrizeFlatte(
    pars::NamedTuple{(:pole_MeV, :g)};
    particle_data=ParticleData(),
    sheet::MomentumSheet{2}=MomentumSheet(:II),
)
    @unpack pole_MeV, g = pars
    (; Ef_MeV, Γ₀_MeV) = pole_parameters(pole_MeV, g; particle_data, sheet)
    return FlatteModel((; Ef_MeV, g, Γ₀_MeV); particle_data)
end

function PoleReparametrizeFlatte(
    pars::NamedTuple{(:pole_re_MeV, :pole_im_MeV, :g)};
    particle_data=ParticleData(),
    sheet::MomentumSheet{2}=MomentumSheet(:II),
)
    @unpack pole_re_MeV, pole_im_MeV, g = pars
    return PoleReparametrizeFlatte((; pole_MeV=pole_re_MeV + 1im * pole_im_MeV, g);
        particle_data, sheet)
end

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

function pole_position(
    model::FlatteModel{<:Tuple{<:DxD,<:DxD,Other}},
    sheet::MomentumSheet{2},
    init=-1e3im * model.Γ₀_MeV / 10,
)
    fr = optimize(
        x -> abs2(denominator(model, x[1] + x[2] * 1im, sheet)),
        collect(reim(init)),
        BFGS(),
    )
    minimum_reached = (fr.minimum < 1e-8)
    !(minimum_reached) && error("Pole is not found: fr.minimum = $(fr.minimum)")
    Epole = complex(fr.minimizer...) # MeV
    return Epole
end
