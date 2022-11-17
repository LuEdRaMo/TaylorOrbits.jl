@doc raw"""
    TaylorOrbits

Taylor methods applied to the orbital determination linkage problem. 
"""
module TaylorOrbits

# Required packages
using PlanetaryEphemeris, NEOs, Healpix, SPICE
using LinearAlgebra, JLD, JLD2, Dates, Printf, Artifacts, DelimitedFiles
using TaylorSeries: set_variables, TaylorN
#using AstroTime: TDBEpoch, TDB, from_utc, j2000, value, julian
using AstroTime
using HTTP: get
using Roots: find_zeros

# Overwritten functions
import Base: hash, ==, show
import Dates: DateTime
import NEOs: w8sveres17

# Exported objects / functions
export TO, k_gauss
export MPCRadec
export ra, dec, read_mpc_file, parse_mpc_obs, search_mpc_circulars, write_mpc_file, 
       w8sveres17, debiasing
export j2000_days, apophis_r_197, apophis_r_199, apophis_v_197, apophis_v_199, gauss, jet_transport

include("constants.jl")
include("parse_mpc.jl")
include("gauss.jl")
include("jet_transport.jl")

end # module TaylorOrbits
