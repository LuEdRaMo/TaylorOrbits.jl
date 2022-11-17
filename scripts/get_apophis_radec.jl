using ArgParse, TaylorOrbits

function search_apophis(m::RegexMatch) 
    if (m["number"] == "99942") || (m["temp_desig"] == "N00hp15")
        if m["obs_code"] ∈ TO.space_obs
            return false
        else
            return true 
        end
    else 
        return false 
    end
end

function parse_commandline()
    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "get_apophis_radec.jl"  
    # Desciption (for help screen)
    s.description = "Search MPC circulars for Apophis (99942/N00hp15) optical observations." 

    @add_arg_table! s begin
        "--url_i"
            help = "URL of MPC circular where to start searching"
            arg_type = String
            default = "https://minorplanetcenter.net/mpec/K20/K20Y98.html"
        "--url_f"
            help = "URL of MPC circular where to finish searching"
            arg_type = String
            default = "https://minorplanetcenter.net/mpec/K21/K21JL0.html"
        "--max_iter"
            help = "Maximum number of circulars to check"
            arg_type = Int
            default = 10_000
        "--output"
            help = "Filename to save the observations found"
            arg_type = String
            default = "data/N00hp15.txt"
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()

    println("Searching Apophis (99942/N00hp15) observations from ")
    println(parsed_args["url_i"])
    println("to")
    println(parsed_args["url_f"])

    radecs = search_mpc_circulars(
        search_apophis, 
        parsed_args["url_i"],
        parsed_args["url_f"],
        max_iter = parsed_args["max_iter"]
    )

    println(length(radecs), " observations found")
    println("Saving output to: ", parsed_args["output"])

    write_mpc_file(radecs, parsed_args["output"])
    #jldsave("data/apdc.jld2"; radecs = radecs)
end

main()
