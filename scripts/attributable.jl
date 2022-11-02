using TaylorOrbits, TaylorSeries, JLD

radecs = MPCRadec("data/99924.txt");

A1 = Attributable(radecs[1:5]; order = 2)
C1 = Curvature(A1)

qs = load("data/sunearthmoon.jld")
q_vec_earth, q_earth, v_vec_earth, v_earth = ephemeris(qs, A1)

#r(ρ, v_ρ) = sqrt(ρ^2 + 2*ρ*(v_vec_earth⋅C1.ρ_unit) + q_earth^2)
#v_r2(ρ, v_ρ) = v_ρ^2 + 2*v_ρ*(v_vec_earth⋅C1.ρ_unit) + ρ^2*C1.η^2 + 2*ρ*(A1.v_α*(v_vec_earth⋅C1.ρ_unit_α) + A1.v_δ*(v_vec_earth⋅C1.ρ_unit_δ)) + v_earth^2
#E(ρ, v_ρ) = v_r2(ρ, v_ρ)/2 - k_gauss^2/r(ρ, v_ρ)
#E_E(ρ, v_ρ) = v_ρ^2/2 + ρ^2*C1.η^2 - 2*k_gauss^2*μ_ES/ρ

#= 
x = 0:0.01:18
y = -1:0.01:1

z = [E(xi, yi)(ts[1]) for xi in x, yi in y]

c = Contour.contour(x,y,z,0)

l = first(lines(c))

xs, ys = coordinates(l)

plot(xs, ys)

=#

#=

moon = qs["moon"] .- qs["sun"]
q_vec_moon = moon[n_t, 1:3]
q_moon = sqrt(q_vec_moon[1]^2 + q_vec_moon[2]^2 + q_vec_moon[3]^2) 
v_vec_moon = moon[n_t, 4:6]
v_moon = sqrt(v_vec_moon[1]^2 + v_vec_moon[2]^2 + v_vec_moon[3]^2) 

ρ_vec_moon = q_vec_moon .- q_vec_earth
ρ_moon = sqrt(ρ_vec_moon[1]^2 + ρ_vec_moon[2]^2 + ρ_vec_moon[3]^2)
ρ_unit_moon = ρ_vec_moon ./ ρ_moon

r_moon(ρ, v_ρ) = sqrt(ρ^2 + 2*ρ*(v_vec_earth⋅ρ_unit_moon) + q_earth^2)
v_r2_moon(ρ, v_ρ) = v_ρ^2 + 2*v_ρ*(v_vec_earth⋅ρ_unit_moon) + ρ^2*C1.η^2 + 2*ρ*(A1.v_α*(v_vec_earth⋅ρ_unit_α) + A1.v_δ*(v_vec_earth⋅ρ_unit_δ)) + v_earth^2

=#

#=

const objname = "Apophis"
dynamics = TO.RNp1BP_pN_A_J23E_J2S_ng_eph!

dense = true
quadmath = false
const opticalobsfile = ""
const radarobsfile = ""
const debias_table = "2018"
const dynamics = RNp1BP_pN_A_J23E_J2S_ng_eph!
const maxsteps = 10000
const t0 = 0.0 # integration initial time
const nyears = 6.0
nyears = 1/365
ss_eph_file = load("/home/luiseduardo/Documentos/UNAM/Física/Servicio Social/sseph343ast016_p31y_et_MIZTLI.jld")["ss16ast_eph"]

r_2_vec = [0.09550965843229045, -0.18467188116225106, -0.1060382517993866]
v_2_vec = [0.09291637881166026, -0.18254939161332126, -0.1058214440445418]
M = PlanetaryEphemeris.t2c_jpl_de430(radecs[8193].date*daysec)
A1 = Attributable(radecs[8192:8194]; order = 2)
q_vec_earth(A1.t[2]) + M * r_2_vec

q0 = [-0.7809669915560347, -0.5476398185159349, -0.3122713504025452, -0.13738345751709308,
      0.16470411192740134, -0.003717355707863583, 0.0, 0.0]

qs = load("data/sunearthmoon.jld")
q_vec_earth, q_earth, v_vec_earth, v_earth = ephemeris(qs, A1)

q0 = vcat(q_vec_earth(A1.t[2]) + M * r_2_vec, v_vec_earth(A1.t[2]) + M * v_2_vec, [0, 0])
dq = set_variables("q", order=5, numvars=8)

q0 = q0 + dq

const μ_DE430 = PlanetaryEphemeris.μ
const μ_B16_DE430 = μ_DE430[12:27] # DE430 GM's of 16 most massive asteroids
const μ_ast343_DE430 = μ_DE430[12:end] # DE430 GM's of 343 main belt asteroids included in DE430 integration
const abstol = 1.0E-30

abstol = 1e-20 #the absolute tolerance of the integration

TO.propagate(objname, dynamics, maxsteps, t0, nyears, ss_eph_file, dense=dense, q0 = q0, quadmath=quadmath, radarobsfile=radarobsfile, opticalobsfile=opticalobsfile, debias_table=debias_table, μ_ast = μ_ast343_DE430, abstol = abstol, order = 25)
=#