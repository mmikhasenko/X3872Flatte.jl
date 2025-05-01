module X3872Flatte

import Base: denominator
using Parameters
using NLsolve
using QuadGK
using Optim

export e2m, m2e
export fm_times_GeV
include("masseswidths.jl")

export FlatteModel
export ReparametrizeFlatte
export compute_corrected_Ef
export AJψππ, denominator
export scattering_parameters
export pole_position
include("flatte.jl")

export BXρ, BXω
export dRρ, dRω, dRDˣ⁰D⁰, dR0
include("branchings.jl")

include("plausible_parameters.jl")

end # module
