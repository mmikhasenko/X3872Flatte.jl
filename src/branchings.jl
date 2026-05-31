function breakup(m′, E, particle_data::ParticleData)
    m0 = e2m(E, particle_data)
    @unpack mJψ = particle_data
    return sqrt(m0 - (m′ + mJψ)) *
           sqrt(m0 + (m′ + mJψ)) *
           sqrt(m0 - (m′ - mJψ)) *
           sqrt(m0 + (m′ - mJψ)) / (2 * m0)
end
# 
upperlimit(E::Number, particle_data::ParticleData) =
    e2m(E, particle_data) - particle_data.mJψ

function BXρ(E::Number, particle_data::ParticleData)
    integrand(m′) = breakup(m′, E, particle_data) * particle_data.Γρ /
                    ((particle_data.mρ - m′)^2 + particle_data.Γρ^2 / 4)
    value = quadgk(integrand, 2 * particle_data.mπ, upperlimit(E, particle_data))[1]
    return value / (2π)
end

function BXω(E::Number, particle_data::ParticleData)
    # 
    @unpack mω, Γω, mπ = particle_data
    integrand(m′) =
        breakup(m′, E, particle_data) * Γω /
        ((mω - m′)^2 + Γω^2 / 4)
    #
    value = quadgk(integrand, 3 * mπ, upperlimit(E, particle_data))[1]
    return value / (2π)
end
# 

function dRρ(model::FlatteModel)
    E_cutoff = 20
    @unpack particle_data = model
    value = quadgk(E -> BXρ(E, particle_data) / abs2(denominator(model, E)),
        -E_cutoff, E_cutoff)[1]
    return model.fρ * value / 1e3  # 1e3 MeV is dE jacobian
end

function dRω(model::FlatteModel)
    E_cutoff = 20
    @unpack particle_data = model
    value = quadgk(E -> BXω(E, particle_data) / abs2(denominator(model, E)),
        -E_cutoff, E_cutoff)[1]
    return model.fω * value / 1e3  # 1e3 MeV is dE jacobian
end

function dRDˣ⁰D⁰(model::FlatteModel)
    E_cutoff = 20
    @unpack particle_data = model
    value = quadgk(E -> real(k1(E, particle_data)) / abs2(denominator(model, E)),
        0.0, E_cutoff)[1]
    return model.g * value / 1e3  # 1e3 MeV is dE jacobian
end

function dR0(model::FlatteModel)
    E_cutoff = 20
    value = quadgk(E -> 1 / abs2(denominator(model, E)),
        -E_cutoff, E_cutoff)[1]
    return model.Γ₀ * value / 1e3  # 1e3 MeV is dE jacobian
end
