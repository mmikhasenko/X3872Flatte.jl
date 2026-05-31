module X3872Flatte

import Base: denominator
using Parameters
using NLsolve
using QuadGK
using Optim

include("masseswidths.jl")

export ParticleData
export FlatteModel
export AbstractChannel, ElasticChannel, InelasticChannel
export JpsiRho, JpsiOmega, DxD, Other, threshold
export ReparametrizeFlatte
export compute_corrected_Ef
export AJψππ, denominator
export scattering_parameters
export pole_position
include("flatte.jl")

include("branchings.jl")

include("plausible_parameters.jl")

end # module
