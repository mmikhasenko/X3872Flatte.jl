using X3872Flatte
using Parameters
using Plots

default_model = let
    @unpack Ef_MeV, g, Γ₀_MeV, fρ, fω = X3872Flatte.LHCb_point8
    FlatteModel(;
        Ef_MeV,     # MeV, binding energy
        g,      # Coupling to D*0D0
        Γ₀_MeV,     # MeV, decay width to other channels
        fρ,     # Coupling to J/ψρ
        fω,      # Coupling to J/ψω
    )
end

# Find pole position
pole = pole_position(default_model)
println("Pole position: $(real(pole)) - $(abs(imag(pole)))i MeV")

# Calculate amplitude over energy range
energies = range(-2, 2, length = 500) # MeV
amplitudes = [AJψππ(default_model, E) for E in energies]

theme(:boxed)
# Plot the amplitude
plot(energies, abs2.(amplitudes),
    xlabel = "m(J/ψπ⁺π⁻)-m(thr.) [MeV]",
    ylabel = "|A|²",
    title = "χc1(3872) Amplitude",
    fillalpha = 0.2, fill = 0)


# Naive calculation of the branching ratios
br_rho = dRρ(default_model)
br_omega = dRω(default_model)
br_dd = dRDˣ⁰D⁰(default_model)
println("Branching ratio to J/ψρ: $br_rho")
println("Branching ratio to J/ψω: $br_omega")
println("Branching ratio to D*0D0: $br_dd")
