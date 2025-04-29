# X3872Flatte.jl

A Julia package for modeling the $\chi_{c1}(3872)$ state (formerly known as X(3872)) using a non-relativistic Flatte parameterization.

## Description

This package implements a Flatte model for the $\chi_{c1}(3872)$ exotic hadron state, based on the parametrization described in [arXiv:0704.0605](https://arxiv.org/abs/0704.0605). The model accounts for the peculiar properties of this state, which lies extremely close to the $D^{*0}D^0$ threshold and exhibits both $J/\psi\pi\pi$ and $D^{*0}D^0$ decay modes.

Key features:
- Implementation of the Flatte model amplitude for $\chi_{c1}(3872)$
- Calculation of pole positions and scattering parameters
- Coupling to various decay channels including $J/\psi\rho$, $J/\psi\omega$, and $D^{*0}D^0$
- Branching ratio calculations for different decay modes

## Installation

Install the package using Julia's package manager:

```julia
using Pkg
Pkg.add("https://github.com/mmikhasenko/X3872Flatte.jl")
```

Or from the Julia REPL:

```julia
] add https://github.com/mikhailmikhasenko/X3872Flatte.jl
```

## Usage

Here's a basic example of how to use the package:

```julia
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
```

Naive calculation of the branching ratios is done as follows:

```julia
br_rho = dRρ(default_model)    # gives 0.49
br_omega = dRω(default_model)  # gives 0.45
br_dd = dRDˣ⁰D⁰(default_model) # gives 2.16
```

The pole position in the complex E plane is calculated using the `pole_position` function:

```julia
pole_E_MeV = pole_position(default_model)
pole_mass = real(pole_E_MeV)    # gives -0.02 (MeV)
pole_width = 2imag(-pole_E_MeV) # gives  0.05 (MeV)
```

## References

The model is 
- discussed in details in [[Hanhart:2007yq]](https://inspirehep.net/literature/747969),
- used in the LHCb analysis [[LHCb:2020xds]](https://inspirehep.net/literature/1798038), and
- in BESIII analysis [[BESIII:2023hml]](https://inspirehep.net/literature/2693487).

## License

MIT

