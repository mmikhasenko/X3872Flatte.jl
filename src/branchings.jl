

function breakup(m′, E)
    m0 = e2m(E)
    return sqrt(m0 - (m′ + mJψ)) *
           sqrt(m0 + (m′ + mJψ)) *
           sqrt(m0 - (m′ - mJψ)) *
           sqrt(m0 + (m′ - mJψ)) / (2 * m0)
end
# 
upperlimit(E::Number) = e2m(E) - mJψ

function BXρ(E::Number)
    integrand(m′) = breakup(m′, E) * Γρ / ((mρ - m′)^2 + Γρ^2 / 4)
    value = quadgk(integrand, 2mπ, upperlimit(E))[1]
    return value / (2π)
end

function BXω(E::Number)
    # 
    integrand(m′) =
        breakup(m′, E) * Γω / ((mω - m′)^2 + Γω^2 / 4)
    #
    value = quadgk(integrand, 3mπ, upperlimit(E))[1]
    return value / (2π)
end
# 

function dRρ(model::FlatteModel)
    E_cutoff = 20
    value = quadgk(E -> BXρ(E) / abs2(denominator(model, E)),
        -E_cutoff, E_cutoff)[1]
    return model.fρ * value / 1e3  # 1e3 MeV is dE jacobian
end

function dRω(model::FlatteModel)
    E_cutoff = 20
    value = quadgk(E -> BXω(E) / abs2(denominator(model, E)),
        -E_cutoff, E_cutoff)[1]
    return model.fω * value / 1e3  # 1e3 MeV is dE jacobian
end

function dRDˣ⁰D⁰(model::FlatteModel)
    E_cutoff = 20
    value = quadgk(E -> real(k1(E)) / abs2(denominator(model, E)),
        0.0, E_cutoff)[1]
    return model.g * value / 1e3  # 1e3 MeV is dE jacobian
end

function dR0(model::FlatteModel)
    E_cutoff = 20
    value = quadgk(E -> 1 / abs2(denominator(model, E)),
        -E_cutoff, E_cutoff)[1]
    return model.Γ₀ * value / 1e3  # 1e3 MeV is dE jacobian
end
