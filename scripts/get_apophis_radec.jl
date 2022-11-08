using TaylorOrbits, JLD2

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

radecs = iterate_mpc_circulars(
    search_apophis, 
    "https://minorplanetcenter.net/mpec/K20/K20Y98.html",
    "https://minorplanetcenter.net/mpec/K21/K21JL0.html"
)

jldsave("data/apdc.jld2"; radecs = radecs)