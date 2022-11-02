using TaylorOrbits, LinearAlgebra, JLD, TaylorSeries, SPICE, Roots

furnsh("naif0012.tls")

radecs = MPCRadec("data/99924.txt");

obs1 = radecs[1]
obs2 = radecs[2]

R_1_vec = TO.observer_position(obs1.obs_code, obs1.date*TO.daysec)[1] / TO.au
R_1 = norm(R_1_vec)
R_1_unit = R_1_vec / R_1
R_2_vec = TO.observer_position(obs2.obs_code, obs2.date*TO.daysec)[1] / TO.au
R_2 = norm(R_2_vec)
R_2_unit = R_2_vec / R_2

qs = load("data/sunearthmoon.jld")
t = qs["t"] ./ TO.daysec
n_1 = searchsortedlast(t, obs1.date)
n_2 = searchsortedlast(t, obs2.date)
δt_1 = obs1.date - t[n_1]
δt_2 = obs2.date - t[n_2]

earth1 = qs["earth"][n_1, :]
sun1 = qs["sun"][n_1, :]
earth2 = qs["earth"][n_2, :]
sun2 = qs["sun"][n_2, :]

order = qs["earth"][1].order
earth1 = map(x -> x(Taylor1(order)*TO.daysec), earth1)
sun1 = map(x -> x(Taylor1(order)*TO.daysec), sun1)
earth2 = map(x -> x(Taylor1(order)*TO.daysec), earth2)
sun2 = map(x -> x(Taylor1(order)*TO.daysec), sun2)

q_1_vec = (earth1[1:3] .- sun1[1:3])(δt_1)
q_1 = norm(q_1_vec)
q_1_unit = q_1_vec / q_1

q_2_vec = (earth2[1:3] .- sun2[1:3])(δt_2)
q_2 = norm(q_2_vec)
q_2_unit = q_2_vec / q_2

ρ_1_unit = TO.ρ(obs1)
ρ_2_unit = TO.ρ(obs2)

cos_0 = ρ_1_unit ⋅ R_1_unit

#= M_1 = [
    -0.225117 0.807296 0.545523;
    -0.433277 -0.584434 0.686081;
    0.872693 -0.081914 0.481348
] =#
M_1 = [
    1.0           0.000886599  0.000385238;
 -0.000886611   1.0          3.20075e-5;
 -0.000385209  -3.23491e-5   1.0
]
#= M_2 = [
    -0.27669 0.791934 0.54432;
    -0.447371 -0.607462 0.656391;
    0.850472 -0.061896 0.522367
] =#
M_2 = [
    1.0           0.000886601  0.000385239;
 -0.000886613   1.0          3.20076e-5;
 -0.00038521   -3.23491e-5   1.0
]

s_1_vec(ρ_1) = M_1 * (R_1_vec + ρ_1 * ρ_1_unit)
s_1(ρ_1) = norm(s_1_vec(ρ_1))
s_1_unit(ρ_1) = s_1_vec(ρ_1) / s_1(ρ_1)
cos_1(ρ_1) = q_1_unit ⋅ s_1_unit(ρ_1)
a2(ρ_1) = s_1(ρ_1)^2 - 2*q_1*s_1(ρ_1)*cos_1(ρ_1) + q_1^2
a(ρ_1) = sqrt(s_1(ρ_1)^2 - 2*q_1*s_1(ρ_1)*cos_1(ρ_1) + q_1^2)

s_2_vec(ρ_2) = M_2 * (R_2_vec + ρ_2 * ρ_2_unit)
s_2(ρ_2) = norm(s_2_vec(ρ_2))
s_2_unit(ρ_2) = s_2_vec(ρ_2) / s_2(ρ_2)
cos_2(ρ_2) = q_2_unit ⋅ s_2_unit(ρ_2)

ρ_2(ρ_1) = find_zeros(x -> s_2(x)^2 - 2*q_2*s_2(x)*cos_2(x) + q_2^2 - a2(ρ_1), (0, 10))[1]

r_1_vec(ρ_1) = q_1_vec + s_1(ρ_1)*s_1_unit(ρ_1)
r_2_vec(ρ_1) = q_2_vec + s_2(ρ_2(ρ_1))*s_2_unit(ρ_2(ρ_1))

T = TO.k_gauss * (obs2.date - obs1.date)

diff(ρ_1) = norm(r_1_vec(ρ_1)×r_2_vec(ρ_1)) - a2(ρ_1)*sin(T*a(ρ_1)^(-3/2))