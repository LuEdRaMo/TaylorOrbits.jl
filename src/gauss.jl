
topo_pos(α::T, δ::T) where {T <: Number} = [cos(δ)*cos(α), cos(δ)*sin(α), sin(δ)]

@doc raw"""
    topo_pos(obs::MPCRadec)

Returns the topocentric unit vector of the body observed in `obs`.
"""
topo_pos(obs::MPCRadec) = topo_pos(obs.α, obs.δ)

function earth_position(eph::TaylorInterpolant, t::T) where {T <: AbstractFloat}
    return eph(t*daysec)[earthdofs[1:3]]
end

function lagrange(A::T, B::T, E::T, F::T, μ::T) where {T <: AbstractFloat}
    a = -(A^2 + A*E + F)
    b = -μ*(2*A*B + B*E)
    c = -(μ*B)^2
    return a, b, c
end

function observer_position(obs::MPCRadec, t::T) where {T <: AbstractFloat} 
    return NEOs.observer_position(obs.obs_code, t*daysec)[1] / au
end

f_Lagrange(u, τ) = 1 - u * τ^2 / 2
g_Lagrange(u, τ) = τ - u * τ^3 / 6

apophis_r_197(t::T) where {T <: Real} = NEOs.apophis_pv_197(t*daysec)[1:3] / au
apophis_r_199(t::T) where {T <: Real} = NEOs.apophis_pv_199(t*daysec)[1:3] / au

apophis_v_197(t::T) where {T <: Real} = NEOs.apophis_pv_197(t*daysec)[4:6] / au * daysec
apophis_v_199(t::T) where {T <: Real} = NEOs.apophis_pv_199(t*daysec)[4:6] / au * daysec

@doc raw"""
    gauss(obs::Vector{MPCRadec{T}}, eph::TaylorInterpolant, μ::T = 1.; 
          max_iter::Int = 1_000, tol::T = 1e-8) where {T <: AbstractFloat}

Gauss method. 
"""
function gauss(obs::Vector{MPCRadec}, eph::PE.TaylorInterpolant, μ::T = 1.; 
               max_iter::Int = 1_000, abstol::T = 1e-8, r_range = (0, 10)) where {T <: AbstractFloat}

    @assert length(obs) == 3 "Gauss method requires three observations"

    println("***Gauss method of orbit determination***")

    # Topocentric position of the body 
    L_1_vec = topo_pos(obs[1])
    L_2_vec = topo_pos(obs[2])
    L_3_vec = topo_pos(obs[3])
    # Days since j2000 epoch
    t_1 = j2000_days(obs[1])
    t_2 = j2000_days(obs[2])
    t_3 = j2000_days(obs[3])
    # Time intervals 
    τ_1 = k_gauss*(t_1 - t_2)
    τ_3 = k_gauss*(t_3 - t_2)
    τ = τ_3 - τ_1
    # Geocentric position of the observer 
    g_1_vec = observer_position(obs[1], t_1)
    g_2_vec = observer_position(obs[2], t_2)
    g_3_vec = observer_position(obs[3], t_3)
    # Barycentric position of the earth 
    G_1_vec = earth_position(eph, t_1)
    G_2_vec = earth_position(eph, t_2)
    G_3_vec = earth_position(eph, t_3)
    # Barycentric position of the observer
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
        println("More than one solution found, choose one:")
        j_ = readline()
        j = parse(Int, j_)
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

    println("3.- Initial solution: ")
    println("r_2_vec = ", r_2_vec)
    println("v_2_vec = ", k_gauss*v_2_vec)
    
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
        
    end

    println("After ", n, " iterations, the slant ranges did not change by more than ", abstol)
    println("5.- Final solution: ")
    println("r_2_vec = ", r_2_vec)
    println("v_2_vec = ", k_gauss*v_2_vec)

    return r_2_vec, k_gauss*v_2_vec
end