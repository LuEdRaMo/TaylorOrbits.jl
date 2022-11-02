using TaylorOrbits, SPICE, JLD2, NEOs, Plots, LinearAlgebra

furnsh("naif0012.tls")
furnsh("a99942/a99942_s197.bsp")

radecs = MPCRadec("data/99924.txt");

eph = JLD2.load("sseph313ast016_p31y_et_MIZTLI.jld2", "ss16ast_eph")

f(t) = NEOs.apophis_pv_197(t*86_400)[1:3]/TO.au

Δs = []
t1s = []
t2s = []

N_obs = length(radecs)

for i_1 ∈ 8_000:N_obs, i_2 ∈ 8_000:N_obs, i_3 ∈ 8_000:N_obs
    if !(i_3 > i_2 > i_1)
        continue
    end

    if (i_1 == i_2) || (i_1 == i_3) || (i_2 == i_3)
        continue
    end

    arc = radecs[[i_1, i_2, i_3]]

    t1 = TO.j2000_days(arc[2]) - TO.j2000_days(arc[1])
    t2 = TO.j2000_days(arc[3]) - TO.j2000_days(arc[2])

    if (t1 > 5) || (t2 > 5)
        continue
    end

    push!(t1s, t1)
    push!(t2s, t2)

    gauss_pos = gauss(arc, eph)[1]

    t = TO.j2000_days(arc[2])
    jpl_pos = f(t)

    Δ = norm(gauss_pos - jpl_pos)

    push!(Δs, Δ)
end

# gauss(radecs[[19, 40, 44]], eph)