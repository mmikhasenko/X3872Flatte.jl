function breakup(m′, E, mJψ, reference_mass)
    m0 = E * 1e-3 + reference_mass
    return sqrt(m0 - (m′ + mJψ)) *
           sqrt(m0 + (m′ + mJψ)) *
           sqrt(m0 - (m′ - mJψ)) *
           sqrt(m0 + (m′ - mJψ)) / (2 * m0)
end

breakup(m′, E, particle_data::ParticleData) =
    breakup(m′, E, particle_data.mJψ, particle_data.mDˣ⁰ + particle_data.mD⁰)
# 
upperlimit(E::Number, mJψ, reference_mass) =
    E * 1e-3 + reference_mass - mJψ

upperlimit(E::Number, particle_data::ParticleData) =
    upperlimit(E, particle_data.mJψ, particle_data.mDˣ⁰ + particle_data.mD⁰)

function BXρ(E::Number, particle_data::ParticleData)
    return BXρ(E, JpsiRho(particle_data), particle_data.mDˣ⁰ + particle_data.mD⁰)
end

function BXρ(E::Number, channel::JpsiRho, reference_mass)
    integrand(m′) = breakup(m′, E, channel.mJψ, reference_mass) * channel.Γρ /
                    ((channel.mρ - m′)^2 + channel.Γρ^2 / 4)
    value = quadgk(integrand, 2 * channel.mπ, upperlimit(E, channel.mJψ, reference_mass))[1]
    return value / (2π)
end

function BXω(E::Number, particle_data::ParticleData)
    return BXω(E, JpsiOmega(particle_data), particle_data.mDˣ⁰ + particle_data.mD⁰)
end

function BXω(E::Number, channel::JpsiOmega, reference_mass)
    integrand(m′) =
        breakup(m′, E, channel.mJψ, reference_mass) * channel.Γω /
        ((channel.mω - m′)^2 + channel.Γω^2 / 4)
    #
    value = quadgk(integrand, 3 * channel.mπ, upperlimit(E, channel.mJψ, reference_mass))[1]
    return value / (2π)
end
# 

function dRρ(model::FlatteModel)
    _zero = zero(model.Ef_MeV)
    length(model.channels) < 4 && return _zero
    # 
    _one = one(model.Ef_MeV)
    E_cutoff = 20 * _one
    ρ = model.channels[4]
    value = quadgk(E -> BXρ(E, ρ, neutral_threshold(model)) / abs2(denominator(model, E)),
        -E_cutoff, E_cutoff)[1]
    MeV_dE_jacobian = 1e3
    return model.fρ * value / MeV_dE_jacobian
end

function dRω(model::FlatteModel)
    _zero = zero(model.Ef_MeV)
    length(model.channels) < 5 && return _zero
    # 
    ω = model.channels[5]
    E_cutoff = 20
    value = quadgk(E -> BXω(E, ω, neutral_threshold(model)) / abs2(denominator(model, E)),
        -E_cutoff, E_cutoff)[1]
    MeV_dE_jacobian = 1e3
    return model.fω * value / MeV_dE_jacobian
end

function dRDˣ⁰D⁰(model::FlatteModel)
    D⁰ = model.channels[1]
    E_cutoff = 20
    value = quadgk(E -> real(k(E, D⁰, neutral_threshold(model))) / abs2(denominator(model, E)),
        0.0, E_cutoff)[1]
    return model.g * value / 1e3  # 1e3 MeV is dE jacobian
end

function dR0(model::FlatteModel)
    E_cutoff = 20
    value = quadgk(E -> 1 / abs2(denominator(model, E)),
        -E_cutoff, E_cutoff)[1]
    return model.Γ₀_MeV * value / 1e3  # 1e3 MeV is dE jacobian
end
