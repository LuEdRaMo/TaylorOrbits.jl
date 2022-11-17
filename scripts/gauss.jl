using ArgParse, TaylorOrbits, JLD2, SPICE

furnsh("data/naif0012.tls")
furnsh("data/a99942_s197.bsp")

function parse_commandline()
    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "gauss.jl"  
    # Desciption (for help screen)
    s.description = "Saves in a .jld2 file the gauss arcs obtained from a set of optical observations." 

    @add_arg_table! s begin
        "--obs_file"
            help = "File where to retrieve the optical observations (MPC formatted)"
            arg_type = String
            default = "data/N00hp15.txt"
        "--eph_file"
            help = "File where to retrieve the Solar System ephemeris (.jld2)"
            arg_type = String
            default = "data/sseph313ast016_p31y_et_MIZTLI.jld2"
        "--max_iter"
            help = "Maximum number of iterations in gauss method of IOD"
            arg_type = Int
            default = 0
        "--output_file"
            help = "Filename to save the arcs"
            arg_type = String
            default = "data/N00hp15.jld2"
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()

    println("Reading optical observations from file: ", parsed_args["obs_file"])
    radecs = read_mpc_file(parsed_args["obs_file"]);
    println("Reading SS ephemeris from file: ", parsed_args["eph_file"])
    eph = JLD2.load(parsed_args["eph_file"], "ss16ast_eph");

    N00hp15 = filter(x -> x.temp_desig == "N00hp15", radecs)
    arc = gauss(N00hp15, eph; max_iter = parsed_args["max_iter"])

    jldsave(parsed_args["output_file"]; arc = arc)
    println("Arcs saved to file: ", parsed_args["output_file"])

end

main()
