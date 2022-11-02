@doc raw"""
    TaylorOrbits

Taylor methods applied to the orbital determination linkage problem. 
"""
module TaylorOrbits

# Required packages
using PlanetaryEphemeris, NEOs
using LinearAlgebra
using AstroTime: TDBEpoch, TDB, from_utc, j2000, value
using HTTP: get
using Roots: find_zeros

#using Dates, SatelliteToolbox, DataFrames, RemoteFiles, DelimitedFiles
#using Artifacts, JLD
#using TaylorSeries: Taylor1, differentiate!, update!, constant_term

# Overwritten functions
import Base: show 

# Exported objects / functions
export TO, k_gauss
export MPCRadec
export gauss

include("constants.jl")
include("parse_mpc.jl")
include("gauss.jl")

end # module TaylorOrbits
