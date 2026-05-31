"""
    ParticleData(; ...)

Particle masses and widths used by the X(3872) Flatte model, in GeV.
Default values come from the PDG@2020.
"""
@with_kw struct ParticleData
    mJψ::Float64 = 3096.90e-3 # GeV
    mχc₁::Float64 = 3871.65e-3 # GeV mass from J/ψ mode
    mD⁰::Float64 = 1864.83e-3 # GeV
    mD⁺::Float64 = 1869.58e-3 # GeV
    mπ::Float64 = 139.57e-3 # GeV
    mρ::Float64 = 775.49e-3 # GeV Neutral only mass
    mω::Float64 = 782.65e-3 # GeV
    Γρ::Float64 = 149.1e-3 # GeV Neutral only width
    Γω::Float64 = 8.68e-3 # GeV
    mDˣ⁺::Float64 = 2010.26e-3 # GeV
    mDˣ⁰::Float64 = 2006.85e-3 # GeV
    ΓDˣ⁺::Float64 = 83.4e-6 # GeV
    ΓDˣ⁰::Float64 = 55.2e-6 # GeV
end

m2e(m, particle_data::ParticleData) =
    (m - particle_data.mDˣ⁰ - particle_data.mD⁰) * 1e3

e2m(E, particle_data::ParticleData) =
    E * 1e-3 + particle_data.mDˣ⁰ + particle_data.mD⁰

const fm_times_GeV = 197.3269804e-3

reduced_mass_neutral(particle_data::ParticleData) =
    particle_data.mD⁰ * particle_data.mDˣ⁰ / (particle_data.mD⁰ + particle_data.mDˣ⁰)

reduced_mass_charged(particle_data::ParticleData) =
    particle_data.mD⁺ * particle_data.mDˣ⁺ / (particle_data.mD⁺ + particle_data.mDˣ⁺)

charged_threshold_offset(particle_data::ParticleData) =
    (particle_data.mD⁺ + particle_data.mDˣ⁺) - (particle_data.mDˣ⁰ + particle_data.mD⁰)

k1(E::Complex, particle_data::ParticleData) =
    1im * sqrt(-2 * reduced_mass_neutral(particle_data) * (E * 1e-3))

k2(E::Complex, particle_data::ParticleData) =
    1im * sqrt(-2 * reduced_mass_charged(particle_data) * (E * 1e-3 - charged_threshold_offset(particle_data)))

k1(E::Real, particle_data::ParticleData) = k1(E + 1e-7im, particle_data)
k2(E::Real, particle_data::ParticleData) = k2(E + 1e-7im, particle_data)
