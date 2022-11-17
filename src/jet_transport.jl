@doc raw"""
    GaussJT{T <: AbstractFloat}

A jet transport orbit whose initial conditions were obtained by gauss method of IOD. 

# Fields 

- `obs::Vector{MPCRadec{T}}`: vector of optical observations. 
- `idx::Int`: index of Gauss method middle observation. 
- `t::T`: days since J2000 of middle observation. 
- `r::Vector{T}`: position at `t`. 
- `v::Vector{T}`: velocity at `t`.
- `bwd_dates::Vector{T}`: interpolation dates for backward ephemeris.
- `bwd_eph::PE.TaylorInterpolant{T, TaylorN{T}, 2}`: backward ephemeris.
- `fwd_dates::Vector{T}`: interpolation dates for forward ephemeris.
- `fwd_eph::PE.TaylorInterpolant{T, TaylorN{T}, 2}`: forward ephemeris. 
"""
struct GaussJT{T <: AbstractFloat}
    obs::Vector{MPCRadec{T}}
    idx::Int
    t::T
    r::Vector{T}
    v::Vector{T}
    bwd_dates::Vector{T}
    bwd_eph::PE.TaylorInterpolant{T, TaylorN{T}, 2}
    fwd_dates::Vector{T}
    fwd_eph::PE.TaylorInterpolant{T, TaylorN{T}, 2}
    # Inner constructor
    function GaussArc{T}(obs::Vector{MPCRadec{T}}, idx::Int, t::T, r::Vector{T}, v::Vector{T},
                         bwd_dates::Vector{T}, bwd_eph::PE.TaylorInterpolant{T, TaylorN{T}, 2},
                         fwd_dates::Vector{T}, fwd_eph::PE.TaylorInterpolant{T, TaylorN{T}, 2}) where {T <: Real}
        new{T}(obs, idx, t, r, v, bwd_dates, bwd_eph, fwd_dates, fwd_eph)
    end
end

# Outer constructor
function GaussJT(obs::Vector{MPCRadec{T}}, idx::Int, t::T, r::Vector{T}, v::Vector{T},
                  bwd_dates::Vector{T}, bwd_eph::PE.TaylorInterpolant{T, TaylorN{T}, 2},
                  fwd_dates::Vector{T}, fwd_eph::PE.TaylorInterpolant{T, TaylorN{T}, 2}) where {T <: Real}
    GaussArc{T}(obs, idx, t, r, v, bwd_dates, bwd_eph, fwd_dates, fwd_eph)
end

function show(io::IO, arc::GaussJT{T}) where {T <: AbstractFloat} 
    print(io, length(arc.obs), " observations from ", arc.obs[1].date, " to ", arc.obs[end].date)
end

function (jt::GaussJT{T})(t) where {T <: AbstractFloat}
    @assert jt.bwd_dates[end] ≤ t ≤ jt.fwd_dates[end]
    if t <= jt.t
        return jt.bwd_eph(t)
    else
        return jt.fwd_eph(t)
    end
end

function jet_transport(arc::GaussArc{T}, eph::PE.TaylorInterpolant; varorder::Int = 5, objname::String = "Apophis",
                   maxsteps::Int = 100, lyap::Bool = false, dense::Bool = false,
                   quadmath::Bool = false, dynamics::Function = NEOs.RNp1BP_pN_A_J23E_J2S_ng_eph_threads!,
                   order::Int = 25, abstol::T = 1.0E-20) where {T <: AbstractFloat}
    
    # Years for backward integration
    nyears_bwd = -ceil(abs(j2000_days(arc.obs[arc.idx]) - j2000_days(arc.obs[1])))/yr
    # Years for forward integration
    nyears_fwd = ceil(abs(j2000_days(arc.obs[end]) - j2000_days(arc.obs[arc.idx])))/yr

    # TaylorN variables setup
    # dq: perturbation to nominal initial condition (TaylorN jet transport)
    dq = set_variables("δx", order = varorder, numvars = 8)
    for i in 1:6
        dq[i][1][i] = 1e-8
    end
    dq[7][1][7] = 0.
    dq[8][1][8] = 0.
    
    # Initial conditions
    jd0 = julian_days(arc.obs[arc.idx])
    q00 = vcat(arc.r, arc.v)
    q0 = vcat(q00, 0.0, 0.0) .+ dq

    # Integrator warmup
    NEOs.propagate(objname, dynamics, 1, jd0, nyears_fwd, eph, output=false, dense=dense, q0=q0, quadmath=quadmath, lyap=lyap, order=order, abstol=abstol)
    println("*** Finished warmup")

    # Full jet transport integration
    NEOs.propagate(objname*"_bwd", dynamics, maxsteps, jd0, nyears_bwd, eph, dense=dense, q0=q0, quadmath=quadmath, lyap=lyap, order=order, abstol=abstol, tord=10, niter=5)
    NEOs.propagate(objname*"_fwd", dynamics, maxsteps, jd0, nyears_fwd, eph, dense=dense, q0=q0, quadmath=quadmath, lyap=lyap, order=order, abstol=abstol, tord=10, niter=5)
    println("*** Finished asteroid ephemeris integration")

    bwd = JLD.load(objname*"_bwd_jt.jld")
    fwd = JLD.load(objname*"_fwd_jt.jld")

    jt = GaussJT(
        arc.obs,
        arc.idx,
        arc.t,
        arc.r,
        arc.v,
        bwd["tv"],
        bwd["asteph"],
        fwd["tv"],
        fwd["asteph"]
    )

    JLD2.jldsave(objname*"_gauss_jt.jld2"; jt = jt)

end

# Sun ephemeris et_seconds -> km, km/s
function xvs(eph, t)
    pv = eph(t)[sundofs]
    pv[1:3] = pv[1:3] * au
    pv[4:6] = pv[4:6] * au / daysec
    return pv
end

# Earth ephemeris et_seconds -> km, km/s
function xve(t, eph)
    pv = eph(t)[earthdofs]
    pv[1:3] = pv[1:3] * au
    pv[4:6] = pv[4:6] * au / daysec
    return pv
end

# Asteroid ephemeris et_seconds -> km, km/s
function xva(jt::GaussJT{T}, t) where {T <: AbstractFloat}
    pv = jt(t/daysec)
    pv[1:3] = pv[1:3] * au
    pv[4:6] = pv[4:6] * au / daysec
    return pv
end

function observations_vector(jt::GaussJT{T}) where {T <: AbstractFloat}
    αs = ra.(jt.obs)
    δs = dec.(jt.obs)

    αs = αs .* cos.( NEOs.arcsec2rad.(δs) )

    return vcat(αs, δs)
end

function computed_radec(jt::GaussJT{T}, eph) where {T <: AbstractFloat}
    m = length(jt.obs)
    αs = Vector{TaylorN{T}}(undef, m)
    δs = Vector{TaylorN{T}}(undef, m)
    for i in 1:m
        α, δ = NEOs.radec(
            jt.obs[i].obs_code, 
            jt.obs[i].date;
            xve = t -> xve(eph, t), 
            xvs = t -> xvs(eph, t), 
            xva = t -> xva(jt, t)
        )
        αs[i] = α * cos(NEOs.arcsec2rad(δ))
        δs[i] = δ
    end

    return vcat(αs, δs)
end

function residuals(arc::GaussJT{T}, eph) where {T <: AbstractFloat}
    
    observed = observations_vector(arc)
    α_corr, δ_corr = debiasing(arc.obs)
    corrections = vcat(α_corr, δ_corr)
    computed = computed_radec(arc, eph)
    
    return observed - corrections - computed 
end