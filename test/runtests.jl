using X3872Flatte
using Test

const default_particle_data = ParticleData()

@testset "Testing amplitude" begin
    model = FlatteModel(;
        Ef_MeV = -7.5, g = 0.1, Γ₀_MeV = 1.8,
        fρ = 0.0, fω = 0.0,
        particle_data = default_particle_data) # Ef is above the second threshold
    Ai = AJψππ(model, 1.1)
    @test Ai ≈ 154.17507563874798 - 179.56915419693846im
end

@testset "Channel construction" begin
    model = FlatteModel((;
        Ef_MeV = -7.5, g = 0.1, Γ₀_MeV = 1.8,
        fρ = 0.1, fω = 0.2);
        particle_data = default_particle_data)

    @test model.channels isa Tuple{DxD,DxD,Other,JpsiRho,JpsiOmega}
    @test threshold(model.channels[1]) ≈ default_particle_data.mDˣ⁰ + default_particle_data.mD⁰
    @test threshold(model.channels[2]) ≈ default_particle_data.mDˣ⁺ + default_particle_data.mD⁺
    @test X3872Flatte.neutral_threshold(model) ≈ threshold(model.channels[1])
    @test X3872Flatte.threshold_offset_MeV(model.channels[2], model.channels[1]) ≈
        X3872Flatte.charged_threshold_offset(default_particle_data) * 1e3

    elastic_model = FlatteModel((;
        Ef_MeV = -7.5, g = 0.1, Γ₀_MeV = 1.8);
        particle_data = default_particle_data)

    @test elastic_model.channels isa Tuple{DxD,DxD,Other}
    @test elastic_model.fρ == 0.0
    @test elastic_model.fω == 0.0
    @test_throws MethodError FlatteModel((;
        Ef_MeV = -7.5, g = 0.1, Γ₀_MeV = 1.8,
        fρ = 0.1);
        particle_data = default_particle_data)
end

@testset "Ef shift uses higher elastic thresholds" begin
    g = 0.13
    Ef_corr = -0.04
    model = FlatteModel((; Ef_MeV = 0.0, g, Γ₀_MeV = 0.0);
        particle_data = default_particle_data)
    expected = 1e3 * real(X3872Flatte.contribution_charged(model, Ef_corr))

    @test X3872Flatte.shift_Ef(g, Ef_corr, default_particle_data) ≈ expected
end

@testset "Pole position" begin
    model = FlatteModel(;
        Ef_MeV = -8.8, g = 0.13, Γ₀_MeV = 1.8, fρ = 0.0, fω = 0.0,
        particle_data = default_particle_data)
    E_pole = pole_position(model)
    @test isapprox(E_pole, 0.0134358620 - 0.1152417753im, atol = 1e-7)
end

@testset "Reparametrize Flatte" begin
    g, Γ₀_MeV, Ef_MeV = 0.13, 1.8, -8.8

    _, Ef_corr = compute_corrected_Ef(Ef_MeV, g, default_particle_data)
    model = ReparametrizeFlatte((;
        Ef_corr, g, Γ₀_MeV, fρ = 0.0, fω = 0.0);
        particle_data = default_particle_data)
    @test model.Ef_MeV ≈ Ef_MeV

    particle_data = ParticleData(; mD⁺ = 1.90, mDˣ⁺ = 2.07)
    _, Ef_corr = compute_corrected_Ef(Ef_MeV, g, particle_data; Ef_corr_guess = -0.04)
    model = ReparametrizeFlatte((; Ef_corr, g, Γ₀_MeV, fρ = 0.0, fω = 0.0); particle_data)
    @test model.Ef_MeV ≈ Ef_MeV
end

@testset "denominator with custom particle masses" begin
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

    @test denominator(model, E) ≈ expected
end
