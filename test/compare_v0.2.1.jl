#!/usr/bin/env julia
# Compare key observables against v0.2.1 (run from package root).
using Test
using Parameters

const PKG_ROOT = dirname(@__DIR__)
const OLD_ROOT = get(ENV, "X3872FLATTE_V021", joinpath(PKG_ROOT, "..", "X3872Flatte-v0.2.1"))
isdir(OLD_ROOT) || error("v0.2.1 worktree not found at $OLD_ROOT")

# Current package (this checkout)
using X3872Flatte

# v0.2.1 sources loaded under a separate module name
const OLD_SRC = joinpath(OLD_ROOT, "src")
module X3872Flatte_v021
    import Base: denominator
    using Parameters, NLsolve, QuadGK, Optim
    for f in ("masseswidths.jl", "flatte.jl", "branchings.jl", "plausible_parameters.jl")
        include(joinpath(Main.OLD_SRC, f))
    end
end
const Old = X3872Flatte_v021

@unpack Ef_MeV, g, Γ₀_MeV, fρ, fω = X3872Flatte.LHCb_point8

old_model = Old.FlatteModel(; Ef_MeV, g, Γ₀_MeV, fρ, fω)
new_model = FlatteModel(; Ef_MeV, g, Γ₀_MeV, fρ, fω, particle_data = ParticleData())

E_test = 1.1
E_grid = (-2.0, 0.0, 1.1, 2.0)

println("=== LHCb_point8 comparison (v0.2.1 vs current) ===")
println("denominator(E=$E_test): old=$(Old.denominator(old_model, E_test))")
println("                    new=$(denominator(new_model, E_test))")
println("AJψππ(E=$E_test):     old=$(Old.AJψππ(old_model, E_test))")
println("                    new=$(AJψππ(new_model, E_test))")

old_pole = Old.pole_position(old_model)
new_pole = pole_position(new_model)
println("pole:                old=$old_pole")
println("                    new=$new_pole")

old_br = (ρ = Old.dRρ(old_model), ω = Old.dRω(old_model), dd = Old.dRDˣ⁰D⁰(old_model))
new_br = (ρ = X3872Flatte.dRρ(new_model), ω = X3872Flatte.dRω(new_model), dd = X3872Flatte.dRDˣ⁰D⁰(new_model))
println("dRρ:  old=$(old_br.ρ)  new=$(new_br.ρ)")
println("dRω:  old=$(old_br.ω)  new=$(new_br.ω)")
println("dRD*: old=$(old_br.dd) new=$(new_br.dd)")

@testset "v0.2.1 regression" begin
    @test denominator(new_model, E_test) ≈ Old.denominator(old_model, E_test)
    @test AJψππ(new_model, E_test) ≈ Old.AJψππ(old_model, E_test)
    for E in E_grid
        @test denominator(new_model, E) ≈ Old.denominator(old_model, E)
    end
    @test pole_position(new_model) ≈ Old.pole_position(old_model) atol = 1e-7
    @test new_br.ρ ≈ old_br.ρ
    @test new_br.ω ≈ old_br.ω
    @test new_br.dd ≈ old_br.dd
end

println("All comparisons passed.")
