@doc raw"""
    TaylorOrbits

Taylor methods applied to the orbital determination linkage problem. 
"""
module TaylorOrbits

# Required packages
using PlanetaryEphemeris, NEOs
using LinearAlgebra, JLD, JLD2, Dates, Printf
using TaylorSeries: set_variables, TaylorN
#using AstroTime: TDBEpoch, TDB, from_utc, j2000, value, julian
using AstroTime
using HTTP: get
using Roots: find_zeros

# Overwritten functions
import Base: hash, ==, show

# Exported objects / functions
export TO, k_gauss
export MPCRadec
export iterate_mpc_circulars, j2000_days, apophis_r_197, apophis_r_199, apophis_v_197,
       apophis_v_199, gauss, jet_transport

include("constants.jl")
include("parse_mpc.jl")
include("gauss.jl")
include("jet_transport.jl")
include("osculating.jl")

end # module TaylorOrbits
