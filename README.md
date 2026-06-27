# X3872Flatte.jl

A Julia package for modeling the χc1(3872) state (formerly known as X(3872)) using a non-relativistic Flatte parameterization.

## Description

This package implements a Flatte model for the χc1(3872) exotic hadron state, based on the parametrization described in [arXiv:0704.0605](https://arxiv.org/abs/0704.0605).
It's a non-relativistic Breit-Wigner function that parametrizes the energy-dependent width by a sum of
- the D*0 D0 and D*+ D- channels with equal coupling strength,
- the J/ψ ρ and J/ψ ω channels with 2π and 3π resonances, and
- a constant inelastic contribution accumulating all other channels.

## Installation

Install the package using Julia's package manager:

```julia
using Pkg
Pkg.add("https://github.com/mmikhasenko/X3872Flatte.jl")
```

Or from the Julia REPL:

```julia
] add https://github.com/mmikhasenko/X3872Flatte.jl
```

## Usage

A minimal model can be built from the published LHCb point and evaluated directly
as a function of energy relative to the neutral D*0 D0 threshold:

```julia
using X3872Flatte
using Parameters

default_model = let
    @unpack Ef_MeV, g, Γ₀_MeV, fρ, fω = X3872Flatte.LHCb_point8
    FlatteModel(;
        Ef_MeV,     # MeV, binding energy
        g,          # Coupling to D*0D0
        Γ₀_MeV,     # MeV, decay width to other channels
        fρ,         # Coupling to J/ψρ
        fω,         # Coupling to J/ψω
        particle_data = ParticleData(),
    )
end

intensity_at_threshold = abs2(AJψππ(default_model, 0.0))
pole_E_MeV = pole_position(default_model)
rho_fraction = dRρ(default_model)
```

The rendered tutorial is available at
[mmikhasenko.github.io/X3872Flatte.jl](https://mmikhasenko.github.io/X3872Flatte.jl/).

## References

The model is 
- discussed in details in [[Hanhart:2007yq]](https://inspirehep.net/literature/747969),
- used in the LHCb analysis [[LHCb:2020xds]](https://inspirehep.net/literature/1798038), and
- in BESIII analysis [[BESIII:2023hml]](https://inspirehep.net/literature/2693487).

## License

MIT
