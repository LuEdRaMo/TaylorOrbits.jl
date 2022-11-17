using TaylorOrbits
using Test

is_all_spaces(x::String) = filter(!isspace, x) == ""

@testset "`MPCRadec`" begin
    # Test observations 
    obs_1 = "     N00hp15  C2020 12 21.63705611 32 36.184-11 27 00.73         18.42wU     F51"
    obs_2 = "99942         C2020 12 21.63705611 32 36.184-11 27 00.73         18.42wU     F51"

    # Parse observations 
    radec_1 = parse_mpc_obs(obs_1)
    radec_2 = parse_mpc_obs(obs_2)

    @test length(radec_1) == length(radec_2) == 1
    @test eltype(radec_1) == eltype(radec_2) == MPCRadec{Float64}
    @test radec_1 == radec_2
    # Number 
    @test is_all_spaces(radec_1[1].number)
    @test radec_2[1].number == "99942"
    # Temporary designation
    @test radec_1[1].temp_desig == "N00hp15"
    @test is_all_spaces(radec_2[1].temp_desig)

    # MPCRadec -> String -> MPCRadec
    radec_1_ = parse_mpc_obs.(TO.mpc_obs_str.(radec_1))[1]
    radec_2_ = parse_mpc_obs.(TO.mpc_obs_str.(radec_2))[1]
    @test radec_1 == radec_1_
    @test radec_2 == radec_2_
end

@testset "Read / write mpc files" begin
    # Test files 
    local source_file = "../data/99924.txt"
    local target_file = "../data/99924_.txt"

    # Real file 
    radecs = read_mpc_file(source_file)
    # Wrote file 
    write_mpc_file(radecs, target_file)
    # Read file 
    radecs_ = read_mpc_file(target_file)
    # Compare 
    @test radecs == radecs_

    # Test url
    local source_url = "https://minorplanetcenter.net/mpec/K20/K20Y98.html"

    # Get raw text 
    txt = TO.get_raw_html(source_url)
    # Parse text 
    radecs = parse_mpc_obs(txt)
    # Wrote file 
    write_mpc_file(radecs, target_file)
    # Read file 
    radecs_ = read_mpc_file(target_file)
    # Compare 
    @test radecs == radecs_
end