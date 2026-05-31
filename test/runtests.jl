using X3872Flatte
using Test

@testset "Testing amplitude" begin
    model = FlatteModel(
        Ef_MeV = -7.5, g = 0.1, Γ₀_MeV = 1.8,
        fρ = 0.0, fω = 0.0) # Ef is above the second threshold
    # 
    Ai = AJψππ(model, 1.1)
    @test Ai ≈ 154.17507563874798 - 179.56915419693846im
end

@testset "Pole position" begin
    model = FlatteModel(Ef_MeV = -8.8, g = 0.13, Γ₀_MeV = 1.8, fρ = 0.0, fω = 0.0)
    E_pole = pole_position(model)
    @test isapprox(E_pole, 0.0134358620 - 0.1152417753im, atol = 1e-8)
end


@testset "Reparametrize Flatte" begin
    g, Γ₀_MeV, Ef_MeV = 0.13, 1.8, -8.8
    _, Ef_corr = compute_corrected_Ef(Ef_MeV, g)
    model = ReparametrizeFlatte((; g, Γ₀_MeV, Ef_corr, fρ = 0, fω = 0))
    @test model.Ef_MeV ≈ Ef_MeV
end

@testset "Particle masses are particle data" begin
    particle_data = ParticleData(;
        mD⁰ = 1.80,
        mDˣ⁰ = 2.00,
        mD⁺ = 1.85,
        mDˣ⁺ = 2.05,
    )
    model = FlatteModel(;
        Ef_MeV = -7.5,
        g = 0.1,
        Γ₀_MeV = 1.8,
        fρ = 0.0,
        fω = 0.0,
        particle_data,
    )

    E = 1.1
    E_complex = E + 1e-7im
    μ = particle_data.mD⁰ * particle_data.mDˣ⁰ / (particle_data.mD⁰ + particle_data.mDˣ⁰)
    μ⁺ = particle_data.mD⁺ * particle_data.mDˣ⁺ / (particle_data.mD⁺ + particle_data.mDˣ⁺)
    δ⁺ = (particle_data.mD⁺ + particle_data.mDˣ⁺) - (particle_data.mDˣ⁰ + particle_data.mD⁰)
    k1 = 1im * sqrt(-2 * μ * (E_complex * 1e-3))
    k2 = 1im * sqrt(-2 * μ⁺ * (E_complex * 1e-3 - δ⁺))
    expected = (E - model.Ef_MeV) * 1e-3 +
               0.5im * model.g * (k1 + k2) +
               0.5im * model.Γ₀_MeV * 1e-3

    @test model.particle_data === particle_data
    @test denominator(model, E) ≈ expected
    @test denominator(model, E) != denominator(FlatteModel(;
        Ef_MeV = model.Ef_MeV,
        g = model.g,
        Γ₀_MeV = model.Γ₀_MeV,
        fρ = model.fρ,
        fω = model.fω,
    ), E)
end

@testset "Branching terms use particle data" begin
    default_model = FlatteModel(;
        Ef_MeV = -7.5,
        g = 0.1,
        Γ₀_MeV = 1.8,
        fρ = 0.02,
        fω = 0.03,
    )
    varied_particle_data = ParticleData(;
        mJψ = 3.000,
        mπ = 0.120,
        mρ = 0.720,
        mω = 0.760,
        Γρ = 0.120,
        Γω = 0.010,
    )
    varied_model = FlatteModel(;
        Ef_MeV = default_model.Ef_MeV,
        g = default_model.g,
        Γ₀_MeV = default_model.Γ₀_MeV,
        fρ = default_model.fρ,
        fω = default_model.fω,
        particle_data = varied_particle_data,
    )

    @test !isapprox(denominator(varied_model, 1.1), denominator(default_model, 1.1); rtol = 1e-8)
end

@testset "Reparametrize Flatte preserves particle data" begin
    particle_data = ParticleData(; mD⁺ = 1.90, mDˣ⁺ = 2.07)
    g, Γ₀_MeV, Ef_MeV = 0.13, 1.8, -8.8
    _, Ef_corr = compute_corrected_Ef(Ef_MeV, g, -0.04, particle_data)
    model = ReparametrizeFlatte((; g, Γ₀_MeV, Ef_corr, fρ = 0, fω = 0, particle_data))

    @test model.particle_data === particle_data
    @test model.Ef_MeV ≈ Ef_MeV
end
