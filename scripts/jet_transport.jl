using ArgParse, TaylorOrbits, NEOs, Dates, TaylorIntegration, JLD2, 
      PlanetaryEphemeris, TaylorSeries, SPICE

function parse_commandline()
    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "jet_transport.jl"  
    # Desciption (for help screen)
    s.description = "Integrates via jet transport a gauss arc." 

    @add_arg_table! s begin
        "--arc_file"
            help = "File where to retrieve the gauss arc (.jld2)"
            arg_type = String
            default = "data/N00hp15.jld2"
        "--eph_file"
            help = "File where to retrieve the Solar System ephemeris (.jld2)"
            arg_type = String
            default = "data/sseph313ast016_p31y_et_MIZTLI.jld2"
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()

    println("Reading arc from file: ", parsed_args["arc_file"])
    arc = JLD2.load(parsed_args["arc_file"], "arc")

    println("Reading SS ephemeris from file: ", parsed_args["eph_file"])
    eph = JLD2.load(parsed_args["eph_file"], "ss16ast_eph");

    jet_transport(arc, eph)
    
    println("Jet transport saved to file: ", "Apophis_gauss_jt.jld2")

end

main()
