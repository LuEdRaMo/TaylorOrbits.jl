struct GaussJT{T <: Real}
    obs::Vector{MPCRadec}
    idx::Int
    t::T
    r::Vector{T}
    v::Vector{T}
    bwd_dates::Vector{T}
    bwd_eph::PE.TaylorInterpolant{T, TaylorN{T}, 2}
    fwd_dates::Vector{T}
    fwd_eph::PE.TaylorInterpolant{T, TaylorN{T}, 2}
    # Inner constructor
    function GaussArc{T}(obs::Vector{MPCRadec}, idx::Int, t::T, r::Vector{T}, v::Vector{T},
                         bwd_dates::Vector{T}, bwd_eph::PE.TaylorInterpolant{T, TaylorN{T}, 2},
                         fwd_dates::Vector{T}, fwd_eph::PE.TaylorInterpolant{T, TaylorN{T}, 2}) where {T <: Real}
        new{T}(obs, idx, t, r, v, bwd_dates, bwd_eph, fwd_dates, fwd_eph)
    end
end

function GaussJT(obs::Vector{MPCRadec}, idx::Int, t::T, r::Vector{T}, v::Vector{T},
                  bwd_dates::Vector{T}, bwd_eph::PE.TaylorInterpolant{T, TaylorN{T}, 2},
                  fwd_dates::Vector{T}, fwd_eph::PE.TaylorInterpolant{T, TaylorN{T}, 2}) where {T <: Real}
    GaussArc{T}(obs, idx, t, r, v, bwd_dates, bwd_eph, fwd_dates, fwd_eph)
end

function jet_transport(arc::GaussArc, eph::PE.TaylorInterpolant; varorder::Int = 5, objname::String = "Apophis",
                   maxsteps::Int = 100, lyap::Bool = false, dense::Bool = false,
                   quadmath::Bool = false, dynamics::Function = NEOs.RNp1BP_pN_A_J23E_J2S_ng_eph_threads!,
                   order::Int = 25, abstol::T = 1.0E-20) where {T <: Real}
    
    nyears_bwd = -ceil(abs(j2000_days(arc.obs[arc.idx]) - j2000_days(arc.obs[1])))/yr
    nyears_fwd = ceil(abs(j2000_days(arc.obs[end]) - j2000_days(arc.obs[arc.idx])))/yr

    # TaylorN variables setup
    # dq: perturbation to nominal initial condition (TaylorN jet transport)
    dq = set_variables("δx", order = varorder, numvars = 8)
    for i in 1:6
        dq[i][1][i] = 1e-8
    end
    dq[7][1][7] = 1e-14
    dq[8][1][8] = 1e-13
    
    # Initial conditions from Apophis JPL solution #197 at the Dec-17-2020.0 (TDB) epoch
    jd0 = value(julian(arc.obs[arc.idx].date))
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

function show(io::IO, arc::GaussJT{T}) where {T <: AbstractFloat} 
    print(io, length(arc.obs), " observations from ", arc.obs[1].date, " to ", arc.obs[end].date)
end

function (jt::GaussJT)(t)
    @assert jt.bwd_dates[end] ≤ t ≤ jt.fwd_dates[end]
    if t <= jt.t
        return jt.bwd_eph(t)
    else
        return jt.fwd_eph(t)
    end
end

function xvs(t, eph)
    pv = eph(t)[sundofs]
    pv[1:3] = pv[1:3] / au
    pv[4:6] = pv[4:6] / au * daysec
    return pv
end

function xve(t, eph)
    pv = eph(t)[earthdofs]
    pv[1:3] = pv[1:3] / au
    pv[4:6] = pv[4:6] / au * daysec
    return pv
end

function xva(t, jt::GaussJT)
    pv = jt(t/daysec)
    pv[1:3] = pv[1:3] / au
    pv[4:6] = pv[4:6] / au * daysec
    return pv
end

function observed_radec(jt::GaussJT)
    αs = getfield.(jt.obs, :α)
    δs = getfield.(jt.obs, :δ)

    return αs, δs
end

function computed_radec(jt::GaussJT, eph)
    m = length(jt.obs)
    αs = Vector{TaylorN{Float64}}(undef, m)
    δs = Vector{TaylorN{Float64}}(undef, m)
    for i in 1:m
        α, δ = NEOs.radec(
            jt.obs[i].obs_code, 
            DateTime(jt.obs[i].date);
            xve = t -> xve(t, eph), 
            xvs = t -> xvs(t, eph), 
            xva = t -> xva(t, jt))
        αs[i] = NEOs.arcsec2rad(α)
        δs[i] = NEOs.arcsec2rad(δ)
    end

    return αs, δs
end

function residuals(arc::GaussJT)
    # JPL Solution #197
    jpl_pos = apophis_r_197.(arc.obs)
    jpl_vel = apophis_v_197.(arc.obs)
    jpl_obs = vcat.(jpl_pos, jpl_vel)

    dates = j2000_days.(arc.obs)

    pre_ = arc.bwd_eph.(dates[1:arc.idx-1])
    pre = [i[1:6] for i in pre_]

    mid = arc.bwd_eph(dates[arc.idx])[1:6]

    pos_ = arc.fwd_eph.(dates[arc.idx+1:end])
    pos = [i[1:6] for i in pos_]

    com = vcat(pre, mid, pos)

    lin_obs = reduce(vcat, jpl_obs)
    lin_com = reduce(vcat, com)

    res = lin_obs - lin_com

    return res
end