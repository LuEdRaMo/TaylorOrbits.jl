using TaylorOrbits, JLD2, SPICE, AstroTime

radecs = JLD2.load("data/apdc.jld2")["radecs"];
eph = JLD2.load("data/sseph313ast016_p31y_et_MIZTLI.jld2", "ss16ast_eph");

furnsh("data/naif0012.tls")
furnsh("data/a99942_s197.bsp")

last_N00hp15 = findlast(x -> x.temp_desig == "N00hp15", radecs)
arc_1 = gauss(radecs[1:last_N00hp15], eph; max_iter = 0)


