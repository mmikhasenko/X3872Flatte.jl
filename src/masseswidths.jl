# constants come from the PDG@2020
# checked
const mJψ = 3096.90e-3 # GeV 
const mχc₁ = 3871.65e-3 #GeV mass from J/ψ mode
const mD⁰ = 1864.83e-3  # GeV
const mD⁺ = 1869.58e-3  # GeV

#  light mesons
const mπ = 139.57e-3; # to be checked and moved to masses.jl
const mρ = 775.49e-3 # GeV Neutral only mass
const mω = 782.65e-3 #Gev
const Γρ = 149.1e-3 # GeV Neutral only width
const Γω = 8.68e-3 #GeV

const mDˣ⁺ = 2010.26e-3  # GeV 
const mDˣ⁰ = 2006.85e-3  # GeV

m2e(m) = (m - mDˣ⁰ - mD⁰) * 1e3
e2m(E) = E * 1e-3 + mDˣ⁰ + mD⁰

const ΓDˣ⁺ = 83.4e-6
const ΓDˣ⁰ = 55.2e-6

# 
const fm_times_GeV = 197.3269804e-3


const μ = mD⁰ * mDˣ⁰ / (mD⁰ + mDˣ⁰)  # GeV
const μ⁺ = mD⁺ * mDˣ⁺ / (mD⁺ + mDˣ⁺)  # GeV
const δ⁺ = (mD⁺ + mDˣ⁺) - (mDˣ⁰ + mD⁰) # GeV

k1(E::Complex) = 1im * sqrt(-2μ * (E * 1e-3))
k2(E::Complex) = 1im * sqrt(-2μ⁺ * (E * 1e-3 - δ⁺))

k1(E::Real) = k1(E + 1e-7im)
k2(E::Real) = k2(E + 1e-7im)

