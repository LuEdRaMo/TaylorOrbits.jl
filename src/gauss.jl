@doc raw"""
    GaussArc{T <: Real}

Gauss method of IOD applied to a set of optical observations. 

# Fields 

- `obs::Vector{MPCRadec}`: vector of optical observations. 
- `idx::Int`: index of Gauss method middle observation. 
- `t::T`: days since J2000 of middle observation. 
- `r::Vector{T}`: position at `t`. 
- `v::Vector{T}`: velocity at `t`. 
"""
struct GaussArc{T <: Real}
    obs::Vector{MPCRadec}
    idx::Int
    t::T
    r::Vector{T}
    v::Vector{T}
    # Inner constructor
    function GaussArc{T}(obs::Vector{MPCRadec}, idx::Int, t::T, r::Vector{T}, v::Vector{T}) where {T <: Real}
        new{T}(obs, idx, t, r, v)
    end
end

# Outer constructor
function GaussArc(obs::Vector{MPCRadec}, idx::Int, t::T, r::Vector{T}, v::Vector{T}) where {T <: Real}
    GaussArc{T}(obs, idx, t, r, v)
end

j2000_days(arc::GaussArc) = j2000_days.(arc.obs)

@doc raw"""
    topo_pos(α::T, δ::T) where {T <: Number}
    topo_pos(obs::MPCRadec)

Returns the topocentric unit vector corresponding to an optical observation.
"""
topo_pos(α::T, δ::T) where {T <: Number} = [cos(δ)*cos(α), cos(δ)*sin(α), sin(δ)]
topo_pos(obs::MPCRadec) = topo_pos(obs.α, obs.δ)

@doc raw"""
    sun_pv_bar(t::T, eph::PE.TaylorInterpolant) where {T <: Real}

Returns the barycentric position and velocity of the Sun at `t` days since J2000. 
"""
function sun_pv_bar(t::T, eph::PE.TaylorInterpolant) where {T <: Real}
    return eph(t*daysec)[sundofs]
end

@doc raw"""
    earth_pv_bar(t::T, eph::TaylorInterpolant) where {T <: Real}

Returns the barycentric position and velocity of the Earth at `t` days since J2000. 
"""
function earth_pv_bar(t::T, eph::TaylorInterpolant) where {T <: Real}
    return eph(t*daysec)[earthdofs]
end

@doc raw"""
    earth_pv_hel(t::T, eph::TaylorInterpolant) where {T <: Real}

Returns the heliocentric position and velocity of the Earth at `t` days since J2000. 
"""
function earth_pv_hel(t::T, eph::TaylorInterpolant) where {T <: Real}
    sun = sun_pv_bar(t, eph)
    earth = earth_pv_bar(t, eph)
    return earth - sun
end

@doc raw"""
    observer_position(t::T, obs::MPCRadec) where {T <: Real} 

Returns the geocentric position of the observer at `t` days since J2000.
"""
function observer_position(t::T, obs::MPCRadec) where {T <: Real} 
    return NEOs.observer_position(obs.obs_code, t*daysec)[1] / au
end

@doc raw"""
    apophis_r_197(t::T) where {T <: Real}
    apophis_r_197(obs::MPCRadec)    
    apophis_r_197(arc::GaussArc)

Returns the barycentric position of Apophis in JPL #197 solution at `t` days since J2000.
"""
apophis_r_197(t::T) where {T <: Real} = NEOs.apophis_pv_197(t*daysec)[1:3] / au
apophis_r_197(obs::MPCRadec) = apophis_r_197(j2000_days(obs))
apophis_r_197(arc::GaussArc) = apophis_r_197(arc.obs[arc.idx])

@doc raw"""
    apophis_r_199(t::T) where {T <: Real}
    apophis_r_199(obs::MPCRadec)    
    apophis_r_199(arc::GaussArc)

Returns the barycentric position of Apophis in JPL #199 solution at `t` days since J2000.
"""
apophis_r_199(t::T) where {T <: Real} = NEOs.apophis_pv_199(t*daysec)[1:3] / au
apophis_r_199(obs::MPCRadec) = apophis_r_199(j2000_days(obs))
apophis_r_199(arc::GaussArc) = apophis_r_199(arc.obs[arc.idx])

@doc raw"""
    apophis_v_197(t::T) where {T <: Real}
    apophis_v_197(obs::MPCRadec)    
    apophis_v_197(arc::GaussArc)

Returns the barycentric velocity of Apophis in JPL #197 solution at `t` days since J2000.
"""
apophis_v_197(t::T) where {T <: Real} = NEOs.apophis_pv_197(t*daysec)[4:6] / au * daysec
apophis_v_197(obs::MPCRadec) = apophis_v_197(j2000_days(obs))
apophis_v_197(arc::GaussArc) = apophis_v_197(arc.obs[arc.idx])

@doc raw"""
    apophis_v_199(t::T) where {T <: Real}
    apophis_v_199(obs::MPCRadec)    
    apophis_v_199(arc::GaussArc)

Returns the barycentric velocity of Apophis in JPL #199 solution at `t` days since J2000.
"""
apophis_v_199(t::T) where {T <: Real} = NEOs.apophis_pv_199(t*daysec)[4:6] / au * daysec
apophis_v_199(obs::MPCRadec) = apophis_v_199(j2000_days(obs))
apophis_v_199(arc::GaussArc) = apophis_v_199(arc.obs[arc.idx])

@doc raw"""
    lagrange(A::T, B::T, E::T, F::T, μ::T) where {T <: Real}

Returns the coefficients of Lagrange equation.
"""
function lagrange(A::T, B::T, E::T, F::T, μ::T) where {T <: Real}
    a = -(A^2 + A*E + F)
    b = -μ*(2*A*B + B*E)
    c = -(μ*B)^2
    return a, b, c
end

@doc raw"""
    f_Lagrange(u, τ)

Returns the Lagrange f coefficient.
"""
f_Lagrange(u, τ) = 1 - u * τ^2 / 2

@doc raw"""
    g_Lagrange(u, τ)

Returns the Lagrange g coefficient.
"""
g_Lagrange(u, τ) = τ - u * τ^3 / 6

@doc raw"""
    select_closer_197(t_2, x)

Returns the index of the element in `x` that is closest to JPL #197 solution for Apophis at 
`t_2` days since J2000.
"""
function select_closer_197(t_2, x)
    true_r = norm(apophis_r_197(t_2))
    errs = abs.(true_r .- x)
    return findmin(errs)[2]
end

@doc raw"""
    central(arc::Vector{MPCRadec})

Returns the index `j` of the element in `arc` that minimizes 

```
|(arc.obs[end] - arc.obs[j]) - (arc.obs[j] - arc.obs[1])|.
```
"""
function central(arc::Vector{MPCRadec})
    # Array of dates 
    dates = j2000_days.(arc)
    # Difference to first date 
    left = abs.(dates .- dates[1])[2:end-1]
    # Difference to last date 
    right = abs.(dates .- dates[end])[2:end-1]
    # |diff to last date - diff to first date|
    central = abs.(left .- right)   
    # Index of minimum index 
    j = findmin(central)[2] + 1

    return j
end

@doc raw"""
    gauss(obs::Vector{MPCRadec}, eph::PE.TaylorInterpolant, μ::T = 1.; 
          select_r::Function = select_closer_197, select_mid::Function = central,
          max_iter::Int = 1_000, abstol::T = 1e-8, r_range = (0, 10)) where {T <: AbstractFloat}

Gauss method of IOD.

# Arguments 

- `obs::Vector{MPCRadec}`: vector of optical observations.
- `eph::PE.TaylorInterpolant`: solar system ephemeris.
- `μ::T = 1.`: combined mass.
- `select_r::Function`: function to select a solution in case of multiple roots.
- `select_mid::Function`: function to select the middle observation.
- `max_iter::Int`: maximum number of refinement iterations.
- `abstol::T`: absolute tolerance for refinement.
- `r_range`: range where to look for solutions of the Lagrange equation. 
"""
function gauss(obs::Vector{MPCRadec}, eph::PE.TaylorInterpolant, μ::T = 1.; 
               select_r::Function = select_closer_197, select_mid::Function = central,
               max_iter::Int = 1_000, abstol::T = 1e-8, r_range = (0, 10)) where {T <: AbstractFloat}

    m = length(obs)
    @assert m >= 3 "Gauss method requires at least three observations"

    if m == 3
        mid_idx = 2
    else 
        mid_idx = select_mid(obs)
    end

    println("***Gauss method of orbit determination***")

    # Topocentric position of the body 
    L_1_vec = topo_pos(obs[1])
    L_2_vec = topo_pos(obs[mid_idx])
    L_3_vec = topo_pos(obs[end])
    # Days since j2000 epoch
    t_1 = j2000_days(obs[1])
    t_2 = j2000_days(obs[mid_idx])
    t_3 = j2000_days(obs[end])
    # Time intervals 
    τ_1 = k_gauss*(t_1 - t_2)
    τ_3 = k_gauss*(t_3 - t_2)
    τ = τ_3 - τ_1
    # Geocentric position of the observer 
    g_1_vec = observer_position(t_1, obs[1])
    g_2_vec = observer_position(t_2, obs[mid_idx])
    g_3_vec = observer_position(t_3, obs[end])
    # Heliocentric position of the earth 
    G_1_vec = earth_pv_hel(t_1, eph)[1:3]
    G_2_vec = earth_pv_hel(t_2, eph)[1:3]
    G_3_vec = earth_pv_hel(t_3, eph)[1:3]
    # Heliocentric position of the observer
    R_1_vec = -(g_1_vec + G_1_vec)
    R_2_vec = -(g_2_vec + G_2_vec)
    R_3_vec = -(g_3_vec + G_3_vec)
    # Gauss scalars
    D_0 = (L_1_vec × L_2_vec) ⋅ L_3_vec

    println("D_0 coefficient: ")
    println(D_0)
    D = [
        (R_1_vec × L_2_vec)⋅L_3_vec   (R_2_vec × L_2_vec)⋅L_3_vec   (R_3_vec × L_2_vec)⋅L_3_vec;
        (L_1_vec × R_1_vec)⋅L_3_vec   (L_1_vec × R_2_vec)⋅L_3_vec   (L_1_vec × R_3_vec)⋅L_3_vec;
        L_1_vec⋅(L_2_vec × R_1_vec)   L_1_vec⋅(L_2_vec × R_2_vec)   L_1_vec⋅(L_2_vec × R_3_vec)
    ]
    
    E = -2(L_2_vec ⋅ R_2_vec)
    F = R_2_vec ⋅ R_2_vec

    A_1 = τ_3 / τ
    B_1 = (1/6)*A_1*(τ^2 - τ_3^2)
    A_3 = -τ_1 / τ
    B_3 = (1/6)*A_3*(τ^2 - τ_1^2)

    A = -(A_1*D[2, 1] - D[2, 2] + A_3*D[2, 3])/D_0
    B = -(B_1*D[2, 1] + B_3*D[2, 3])/D_0

    # Lagrange equation coefficients 
    a, b, c = lagrange(A, B, E, F, μ)
    println("1.- Lagrange equation: ")
    println("r^8 + (", a, ")r^6 + (", b, ")r^3 + (", c, ") = 0")
    # Solutions to Lagrange equation 
    x = find_zeros(r -> r^8 + a*r^6 + b*r^3 + c, r_range)
    println("2- Solutions to the Lagrange equation: ")
    println(x)

    # Choose solution 
    n_x = length(x)
    if n_x == 0
        println("No solutions found, cannot proceed")
        return zero(L_1_vec), zero(L_1_vec)
    elseif n_x == 1
        println("One solution found, proceeding")
        r_2 = x[1]
    else
        j = select_r(t_2, x)
        println("More than one solution found, function picked: ", j)
        r_2 = x[j]
    end

    # f, g - Lagrange coefficients
    u_2 = μ / r_2^3

    f_1 = f_Lagrange(u_2, τ_1)
    f_3 = f_Lagrange(u_2, τ_3)

    g_1 = g_Lagrange(u_2, τ_1)
    g_3 = g_Lagrange(u_2, τ_3)

    # c - coefficients
    c_1 = g_3 / (f_1*g_3 - f_3*g_1)
    c_2 = -1
    c_3 = -g_1 / (f_1*g_3 - f_3*g_1)

    # Slant ranges (observer - body distance)
    p_1 = (c_1*D[1, 1] + c_2 * D[1, 2] + c_3 * D[1, 3])/c_1/D_0
    p_2 = (c_1*D[2, 1] + c_2 * D[2, 2] + c_3 * D[2, 3])/c_2/D_0
    p_3 = (c_1*D[3, 1] + c_2 * D[3, 2] + c_3 * D[3, 3])/c_3/D_0
    
    # Barycentric position of the body 
    r_1_vec = p_1 * L_1_vec - R_1_vec
    r_2_vec = p_2 * L_2_vec - R_2_vec
    r_3_vec = p_3 * L_3_vec - R_3_vec

    # Barycentric velocity of the body 
    v_2_vec = (-f_3*r_1_vec + f_1*r_3_vec)/(f_1*g_3 - f_3*g_1)

    # Barycentric position and velocity of the sun 
    pv_sun = sun_pv_bar(t_2, eph)

    println("3.- Initial solution: ")
    println("r_2_vec = ", r_2_vec + pv_sun[1:3])
    println("v_2_vec = ", k_gauss*v_2_vec + pv_sun[4:6])
    
    # Observation times corrected for light-time
    t_c1 = t_1 - p_1/c_au_per_day
    t_c2 = t_2 - p_2/c_au_per_day
    t_c3 = t_3 - p_3/c_au_per_day
    # Time intervals corrected for light-time
    τ_1 = k_gauss*(t_c1 - t_c2)
    τ_3 = k_gauss*(t_c3 - t_c2)
    τ = τ_3 - τ_1

    # Iterative procedure 
    println("4.- Solution refinement")
    # Number of iterations
    n = 0         
    # Slant ranges change 
    diff1 = 1
    diff2 = 1
    diff3 = 1
    # Old slant ranges 
    p_1_old = p_1
    p_2_old = p_2
    p_3_old = p_3

    while (diff1 > abstol) && (diff2 > abstol) && (diff3 > abstol) && (n < max_iter)
        n = n + 1
        # Slant range 
        r_2 = norm(r_2_vec)
        # f, g - Lagrange coefficients
        u_2 = μ / r_2^3

        ff_1 = f_Lagrange(u_2, τ_1)
        ff_3 = f_Lagrange(u_2, τ_3)

        gg_1 = g_Lagrange(u_2, τ_1)
        gg_3 = g_Lagrange(u_2, τ_3)
        
        # Average old and new values 
        f_1 = (f_1 + ff_1) / 2
        f_3 = (f_3 + ff_3) / 2

        g_1 = (g_1 + gg_1) / 2
        g_3 = (g_3 + gg_3) / 2

        # c - coefficients
        c_1 = g_3 / (f_1*g_3 - f_3*g_1)
        c_2 = -1
        c_3 = -g_1 / (f_1*g_3 - f_3*g_1)

        # Slant ranges (observer - body distance)
        p_1 = (c_1*D[1, 1] + c_2 * D[1, 2] + c_3 * D[1, 3])/c_1/D_0
        p_2 = (c_1*D[2, 1] + c_2 * D[2, 2] + c_3 * D[2, 3])/c_2/D_0
        p_3 = (c_1*D[3, 1] + c_2 * D[3, 2] + c_3 * D[3, 3])/c_3/D_0

        # Barycentric position of the body 
        r_1_vec = p_1 * L_1_vec - R_1_vec
        r_2_vec = p_2 * L_2_vec - R_2_vec
        r_3_vec = p_3 * L_3_vec - R_3_vec

        # Barycentric velocity of the body 
        v_2_vec = (-f_3*r_1_vec + f_1*r_3_vec)/(f_1*g_3 - f_3*g_1)

        # Slant ranges change 
        diff1 = abs(p_1 - p_1_old)
        diff2 = abs(p_2 - p_2_old)
        diff3 = abs(p_3 - p_3_old)

        p_1_old = p_1
        p_2_old = p_2
        p_3_old = p_3
        
        # Observation times corrected for light-time
        t_c1 = t_1 - p_1/c_au_per_day
        t_c2 = t_2 - p_2/c_au_per_day
        t_c3 = t_3 - p_3/c_au_per_day
        # Time intervals corrected for light-time
        τ_1 = k_gauss*(t_c1 - t_c2)
        τ_3 = k_gauss*(t_c3 - t_c2)
        τ = τ_3 - τ_1
    end

    # Barycentric position and velocity of the sun 
    pv_sun = sun_pv_bar(t_c2, eph)

    println("After ", n, " iterations, the slant ranges did not change by more than ", abstol)
    println("5.- Final solution: ")
    println("r_2_vec = ", r_2_vec + pv_sun[1:3])
    println("v_2_vec = ", k_gauss*v_2_vec + pv_sun[4:6])

    return GaussArc(
        obs,
        mid_idx,
        t_c2,
        r_2_vec + pv_sun[1:3],
        k_gauss*v_2_vec + pv_sun[4:6]
    )
end

function show(io::IO, arc::GaussArc{T}) where {T <: AbstractFloat} 
    print(io, length(arc.obs), " observations from ", arc.obs[1].date, " to ", arc.obs[end].date)
end