using TaylorOrbits, SPICE, JLD2

furnsh("naif0012.tls")

radecs = MPCRadec("data/99924.txt");

eph = JLD2.load("sseph313ast016_p31y_et_MIZTLI.jld2", "ss16ast_eph")

laplace(radecs[1:5], eph)