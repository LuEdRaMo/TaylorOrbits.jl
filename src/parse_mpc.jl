@doc raw"""
    MPCRadec{T <: AbstractFloat} 

An optical measurement in MPC format. See https://minorplanetcenter.net/iau/info/OpticalObs.html
for a detailed description of the format. 

# Fields 

- `number::String`: object's number. 
- `temp_desig::String`: temporary designation.
- `disc_astk::String`: discovery asterisk.
- `note_1::String`: note 1.
- `note_2::String`: note 2.
- `date::DateTime`: date of observation.
- `α::T`: right ascension [arcsec].
- `δ::T`: declination [arcsec].
- `ref::String`: reference. 
- `mag::String`: observed magnitude. 
- `band::String`: magnitude band. 
- `catalog::String`: catalog. 
- `info_2::String`: additional information. 
- `obs_code::String`: observatory code. 
"""
struct MPCRadec{T <: AbstractFloat} 
    number::String 
    temp_desig::String
    disc_astk::String
    note_1::String
    note_2::String
    date::DateTime
    α::T
    δ::T
    ref::String
    mag::String
    band::String
    catalog::String
    info_2::String
    obs_code::String
    # Inner constructor 
    function MPCRadec{T}(number::String, temp_desg::String, disc_astk::String, note_1::String, 
                         note_2::String, date::DateTime, α::T, δ::T, ref::String, mag::String, 
                         band::String, catalog::String, info_2::String, obs_code::String) where {T <: AbstractFloat}
        new{T}(number, temp_desg, disc_astk, note_1, note_2, date, α, δ, ref, mag, band, 
               catalog, info_2, obs_code)
    end
end

# Outer constructor
function MPCRadec(number::String, temp_desg::String, disc_astk::String, note_1::String, 
                  note_2::String, date::DateTime, α::T, δ::T, ref::String, mag::String, 
                  band::String, catalog::String, info_2::String, obs_code::String) where {T <: AbstractFloat} 
    MPCRadec{T}(number, temp_desg, disc_astk, note_1, note_2, date, α, δ, ref, mag, band,
                catalog, info_2, obs_code)
end

# Two MPCRadec are equal if ther date, α, δ and observatory code are equal
function hash(a::MPCRadec{T}, h::UInt) where {T <: AbstractFloat}
    return hash((a.date, a.α, a.δ, a.obs_code), h)
end

function ==(a::MPCRadec{T}, b::MPCRadec{T}) where {T <: AbstractFloat}
    return hash(a) == hash(b)
end

# Print method for MPCRadec
# Example: 
# id: 99942 α: 2.0482180537085553 δ: 0.4716727578016939 t: 2021-05-12T06:29:45.089 TDB obs: F51
function show(io::IO, m::MPCRadec{T}) where {T <: AbstractFloat} 
    # If there is no number, use temporary designation
    if filter(!isspace, m.number) == ""
        print(io, "id: ", m.temp_desig, " α: ", m.α, " δ: ", m.δ, " t: ", m.date,
              " obs: ", m.obs_code)
    else
        print(io, "id: ", m.number, " α: ", m.α, " δ: ", m.δ, " t: ", m.date,
              " obs: ", m.obs_code)
    end
end

# Regular expression to parse an optical measurement in MPC format
const mpc_line_regex = Regex(join(
    [
        # Number regex (columns 1-5)
        raw"(?P<number>.{5})",
        # Temporary designation regex (columns 6-12)
        raw"(?P<temp_desig>.{7})",
        # Discovery asterisk regex (column 13)
        raw"(?P<disc_astk>.{1})",
        # Note 1 regex (column 14)
        raw"(?P<note_1>.{1})",
        # Note 2 regex (column 15)
        raw"(?P<note_2>.{1})",
        # Year regex + space (columns 16-20)
        raw"(?P<year>\d{4}) ",
        # Month regex + space (columns 21-23)
        raw"(?P<month>\d{2}) ",
        # Day regex (columns 24-25)
        raw"(?P<day>\d{2})",
        # fraction of days regex (columns 26-32) 
        raw"(?P<utc>\.[\d\s]{6})",
        # α hours regex + space (columns 33-35)
        raw"(?P<α_hrs>\d{2}) ",
        # α minutes regex + space (columns 36-38)
        raw"(?P<α_min>\d{2}) ",
        # α seconds regex (columns 39-44)
        raw"(?P<α_sec>\d{2}\.[\d\s]{3})",
        # δ sign regex (column 45)
        raw"(?P<δ_sgn>\+|\-)",
        # δ degrees regex + space (columns 46-48)
        raw"(?P<δ_deg>\d{2}) ",
        # δ minutes regex + space (columns 49-51)
        raw"(?P<δ_min>\d{2}) ",
        # δ seconds regex (columns 52-56)
        raw"(?P<δ_sec>\d{2}\.[\d\s]{2})",
        # Reference regex (columns 57-65)
        raw"(?P<ref>.{9})",
        # Magnitude regex (columns 66-70)
        raw"(?P<mag>.{5})",
        # Band regex (column 71)
        raw"(?P<band>.{1})",
        # Catalog regex (column 72)
        raw"(?P<catalog>.{1})",
        # Info 2 regex (columns 73-77)
        raw"(?P<info_2>.{5})",
        # Observatory code regex (columns 78-80)
        raw"(?P<obs_code>.{3})"
    ]
))

@doc raw"""
    intndec(x::T) where {T <: AbstractFloat}

Returns the integer and decimal part of `x`.
"""
function intndec(x::T) where {T <: AbstractFloat}
    # Integer part 
    x_ = floor(Int, x)
    # Return integer + decimal part 
    return x_, x - x_
end

@doc raw"""
    DateTime(year::Int, month::Int, day::Int, utc::T) where {T <: Real}

Construct a `DateTime` type by parts. `utc` is the fraction of day. 
"""
function DateTime(year::Int, month::Int, day::Int, utc::T) where {T <: Real}
    return DateTime(year, month, day) + Microsecond( round(1e6*86400*utc) )
end

@doc raw"""
    ra(hrs::Int, min::Int, sec::T) where {T <: Real}

Returns the right ascension in arcsec. 
"""
function ra(hrs::Int, min::Int, sec::T) where {T <: Real}
    # Convert hours minutes seconds to arcsec
    α_arcsec = (54_000 * hrs) + (900 * min) + (15 * sec)
    # Convert deg to arcsec
    #α_arcsec = 3_600 * α_deg
    return α_arcsec
end

@doc raw"""
    ra(obs::MPCRadec{T}) where {T <: AbstractFloat}

Returns the right ascension of `obs`.
"""
ra(obs::MPCRadec{T}) where {T <: AbstractFloat} = getfield(obs, :α)

@doc raw"""
    dec(sgn::String, deg::Int, min::Int, sec::T) where {T <: Real}

Returns the declination in arcsec. 
"""
function dec(sgn::SubString{String}, deg::Int, min::Int, sec::T) where {T <: Real}
    # Convert degrees minutes seconds to arcsec
    if sgn == "+"
        δ_arcsec = (3_600 * deg) + (60 * min) + sec
    elseif sgn == "-"
        δ_arcsec = -( (3_600 * deg) + (60 * min) + sec ) 
    end
    return δ_arcsec
end

@doc raw"""
    dec(obs::MPCRadec{T}) where {T <: AbstractFloat}

Returns the declination of `obs`.
"""
dec(obs::MPCRadec{T}) where {T <: AbstractFloat} = getfield(obs, :δ)

@doc raw"""
    MPCRadec(m::RegexMatch)

Returns a `MPCRadec` object. `m` must be a match of `TO.mpc_line_regex`.
"""
function MPCRadec(m::RegexMatch)
    date = DateTime(
        Meta.parse(m["year"]), 
        Meta.parse(m["month"]), 
        Meta.parse(m["day"]),
        Meta.parse(m["utc"])
    )
    α = ra(
        Meta.parse(m["α_hrs"]), 
        Meta.parse(m["α_min"]), 
        Meta.parse(m["α_sec"])
    )
    δ = dec(
        m["δ_sgn"], 
        Meta.parse(m["δ_deg"]), 
        Meta.parse(m["δ_min"]), 
        Meta.parse(m["δ_sec"])
    )
    
    return MPCRadec(
        string(m["number"]),
        string(m["temp_desig"]),
        string(m["disc_astk"]),
        string(m["note_1"]),
        string(m["note_2"]),
        date, 
        α, 
        δ, 
        string(m["ref"]),
        string(m["mag"]),
        string(m["band"]),
        string(m["catalog"]),
        string(m["info_2"]),
        string(m["obs_code"])
    )
end

@doc raw"""
    read_mpc_file(filename::String)

Returns a `Vector{MPCRadec}` corresponding to the matches of `TO.mpc_line_regex` in `filename`.
"""
function read_mpc_file(filename::String)
    # Read lines of mpc formatted file 
    lines = readlines(filename)
    # Apply regular expressions
    matches = match.(mpc_line_regex, lines)
    # Eliminate nothings
    filter!(!isnothing, matches)
    # Convert matches to MPCRadec
    radecs = MPCRadec.(matches)
    # Sort observations by date
    sort!(radecs, by = x -> x.date)
    # Eliminate repeated observations
    unique!(radecs)
    
    return radecs
end

@doc raw"""
    parse_mpc_obs(text::String)
    parse_mpc_obs(f::Function, text::String)

Returns the matches of `TO.mpc_line_regex` in `text`. A function `f` can be passed to filter
the matches. 
"""
function parse_mpc_obs(f::Function, text::String)

    # Vector of observations 
    radecs = Vector{MPCRadec{Float64}}(undef, 0)
    # Iterate over the matches 
    for m in eachmatch(mpc_line_regex, text)
        # Filter by f
        if f(m)
            push!(radecs, MPCRadec(m))
        end
    end
    # If there is at least one observation
    if length(radecs) > 0
        # Sort observations by date
        sort!(radecs, by = x -> x.date)
        # Eliminate repeated observations
        unique!(radecs)
    end
    
    return radecs
end
parse_mpc_obs(text::String) = parse_mpc_obs(t -> true, text)

@doc raw"""
    j2000_days(obs::MPCRadec{T}) where {T <: AbstractFloat}    

Returns the days since J2000 epoch. 
"""
function j2000_days(obs::MPCRadec{T}) where {T <: AbstractFloat}    
    return datetime2julian(obs.date) - PE.J2000
end

function julian_days(obs::MPCRadec{T}) where {T <: AbstractFloat}
    return datetime2julian(obs.date)
end

@doc raw"""
    get_raw_html(url::String = "https://minorplanetcenter.net/mpec/K20/K20YA9.html")

Returns the raw html text of webpage `url`.
"""
function get_raw_html(url::String = "https://minorplanetcenter.net/mpec/K20/K20YA9.html")
    # Get raw html 
    resp = get(url)
    # Convert to string 
    text = String(resp.body)
    return text
end

# MPC main page url 
const mpc_url = "https://minorplanetcenter.net"

# Regex for next circular url 
const next_circular_regex = r"<a href=\"(?P<next>.*)\"><img src=\"/iau/figs/RArrow.gif\""

@doc raw"""
    search_mpc_circulars(
        f::Function,
        url1::String = "https://minorplanetcenter.net/mpec/K20/K20YA9.html",
        url2::String = "https://minorplanetcenter.net/mpec/K21/K21JL0.html";
        max_iter::Int = 10_000
    )


Iterates MPC circulars from `url1` to `url2` and returns the matches of `TO.mpc_line_regex` 
filtered by `f`. If the function do not reach `url2` before `max_iter` iterations, it will 
print a warning and return the matches found so far. 
"""
function search_mpc_circulars(
    f::Function,
    url1::String = "https://minorplanetcenter.net/mpec/K20/K20YA9.html",
    url2::String = "https://minorplanetcenter.net/mpec/K21/K21JL0.html";
    max_iter::Int = 10_000
)

    # Vector of observations
    obs = Vector{MPCRadec{Float64}}(undef, 0)
    
    # Number of urls checked 
    n = 0
    # First url 
    u = url1

    while true 
        n += 1
        if n > max_iter
            @warn("$n pages checked before getting to $url2")
            break 
        end
        # Raw html text of webpage u 
        text = get_raw_html(u)
        # Observations found in text 
        obs_ = parse_mpc_obs(f, text)
        # Add new observations 
        obs = vcat(obs, obs_) 
        # Final url 
        if u == url2
            break
        end
        # Next circular url 
        next = match(next_circular_regex, text)["next"]
        u = mpc_url * next
    end
    # If there is at least one observation
    if length(obs) > 0
        # Sort observations by date 
        sort!(obs, by = x -> x.date)
        # Eliminate repeated observations
        unique!(obs)
    end

    return obs
end

@doc raw"""
    mpc_date_str(date::DateTime)

Returns the date in MPC format. 
"""
function mpc_date_str(date::DateTime)
    
    # Year string 
    year_s = lpad(Dates.year(date), 4)
    # Month string 
    month_s = lpad(Dates.month(date), 2, "0")
    # Hours [days] 
    hrs = Dates.hour(date) / 24
    # Minutes [days] 
    min = Dates.minute(date) / 24 / 60
    # Seconds [days] 
    sec = Dates.second(date) / 24 / 60 / 60
    # Milliseconds [days] 
    mls = Dates.millisecond(date) / 24 / 60 / 60 / 1_000
    # Days 
    day_val = Dates.day(date) + hrs + min + sec + mls
    # Days string 
    day_s = @sprintf("%09.6f", day_val)
    # Join everything 
    date_s = join([
        year_s,
        " ",
        month_s,
        " ",
        day_s,
    ])

    return date_s
end

@doc raw"""
    mpc_α_str(α::T) where {T <: Number}

Returns the right ascension [arcsec] in MPC format. 
"""
function mpc_α_str(α::T) where {T <: Number}
    # Convert arcsec to deg 
    α_deg = α / 3_600
    # Hours 
    hrs, hrs_ = intndec(α_deg / 15)
    # Hours string 
    hrs_s = lpad(hrs, 2, "0")
    # Minutes 
    min, min_ = intndec(60 * hrs_)
    # Minutes string 
    min_s = lpad(min, 2, "0")
    # Seconds 
    sec = 60 * min_
    # Seconds string 
    sec_s = @sprintf("%06.3f", sec)
    # Join everything
    α_s = join([
        hrs_s,
        " ",
        min_s,
        " ",
        sec_s,
    ])

    return α_s 
end

@doc raw"""
    mpc_δ_str(δ::T) where {T <: Number}

Returns the declination [arcsec] in MPC format. 
"""
function mpc_δ_str(δ::T) where {T <: Number}
    # Sign string 
    sgn_s = δ >= 0 ? "+" : "-"
    # Convert arcsec to deg 
    δ_deg = abs(δ / 3_600)
    # Degrees 
    deg, deg_ = intndec(δ_deg)
    # Degrees string 
    deg_s = lpad(deg, 2, "0")
    # Minutes 
    min, min_ = intndec(60 * deg_)
    # Minutes string 
    min_s = lpad(min, 2, "0")
    # Seconds
    sec = 60 * min_
    # Seconds string 
    sec_s = @sprintf("%05.2f", sec)
    # Join everything
    δ_s = join([
        sgn_s,
        deg_s,
        " ",
        min_s,
        " ",
        sec_s,
    ])

    return δ_s 
end

@doc raw"""
    mpc_line_str(obs::MPCRadec{T}) where {T <: AbstractFloat}

Returns an observation in MPC format. 
"""
function mpc_obs_str(obs::MPCRadec{T}) where {T <: AbstractFloat}
    # Date string 
    date_s = mpc_date_str(obs.date)
    # Right ascension string 
    α_s = mpc_α_str(obs.α)
    # Declination string 
    δ_s = mpc_δ_str(obs.δ)
    # Join everything
    obs_s = join([
        obs.number,
        obs.temp_desig,
        obs.disc_astk,
        obs.note_1,
        obs.note_2,
        date_s,
        α_s,
        δ_s,
        obs.ref,
        obs.mag,
        obs.band,
        obs.catalog,
        obs.info_2,
        obs.obs_code,
        "\n"
    ])

    return obs_s
end

@doc raw"""
    write_mpc_file(obs::Vector{MPCRadec{T}}, filename::String) where {T <: AbstractFloat}

Writes `obs` to `filename` in MPC format. 
"""
function write_mpc_file(obs::Vector{MPCRadec{T}}, filename::String) where {T <: AbstractFloat}
    open(filename, "w") do file
        for i in eachindex(obs)
            line = mpc_obs_str(obs[i])
            write(file, line)
        end 
    end
end

function w8sveres17(obs::MPCRadec{T}) where {T <: AbstractFloat}
   return w8sveres17(obs.obs_code, obs.date, obs.catalog)
end

function  debiasing(obs::Vector{MPCRadec{T}}, debias_table::String="2018") where {T <: AbstractFloat}

    m = length(obs)
    α_corr = Vector{Float64}(undef, m)
    δ_corr = Vector{Float64}(undef, m)
    
    # Select debiasing table: 
    # - 2014 corresponds to https://doi.org/10.1016/j.icarus.2014.07.033
    # - 2018 corresponds to https://doi.org/10.1016/j.icarus.2019.113596
    # Debiasing tables are loaded "lazily" via Julia artifacts, according to rules in Artifacts.toml
    if debias_table == "2018"
        debias_path = artifact"debias_2018"
        mpc_catalog_codes_201X = NEOs.mpc_catalog_codes_2018
        # The healpix tesselation resolution of the bias map from https://doi.org/10.1016/j.icarus.2019.113596
        NSIDE= 64 
        # In 2018 debias table Gaia DR2 catalog is regarded as the truth
        truth = "V" 
    elseif debias_table == "hires2018"
        debias_path = artifact"debias_hires2018"
        mpc_catalog_codes_201X = NEOs.mpc_catalog_codes_2018
        # The healpix tesselation resolution of the high-resolution bias map from https://doi.org/10.1016/j.icarus.2019.113596
        NSIDE= 256 
        # In 2018 debias table Gaia DR2 catalog is regarded as the truth
        truth = "V" 
    elseif debias_table == "2014"
        debias_path = artifact"debias_2014"
        mpc_catalog_codes_201X = NEOs.mpc_catalog_codes_2014
        # The healpix tesselation resolution of the bias map from https://doi.org/10.1016/j.icarus.2014.07.033
        NSIDE= 64 
        # In 2014 debias table PPMXL catalog is regarded as the truth
        truth = "t" 
    else
        @error "Unknown bias map: $(debias_table). Possible values are `2014`, `2018` and `hires2018`."
    end
    # Debias table file 
    bias_file = joinpath(debias_path, "bias.dat")
    # Read bias matrix 
    bias_matrix = readdlm(bias_file, comment_char='!', comments=true)
    # Initialize healpix Resolution variable
    resol = Resolution(NSIDE) 
    # Compatibility between bias matrix and resolution 
    @assert size(bias_matrix) == (resol.numOfPixels, 4length(mpc_catalog_codes_201X)) "Bias table file $bias_file dimensions do not match expected parameter NSIDE=$NSIDE and/or number of catalogs in table."

    for i in eachindex(obs)
        α_i_rad = NEOs.arcsec2rad(obs[i].α)
        δ_i_rad = NEOs.arcsec2rad(obs[i].δ)

        if (obs[i].catalog ∉ mpc_catalog_codes_201X) && obs[i].catalog != "Y"
            # Handle case: if star catalog not present in debiasing table, then set corrections equal to zero
            if haskey(NEOs.mpc_catalog_codes, obs[i].catalog)
                if obs[i].catalog != truth
                    catalog_not_found = NEOs.mpc_catalog_codes[obs[i].catalog]
                    @warn "Catalog not found in $(debias_table) table: $(catalog_not_found). Setting debiasing corrections equal to zero."
                end
            elseif obs[i].catalog == " "
                @warn "Catalog information not available in observation record. Setting debiasing corrections equal to zero."
            else
                @warn "Catalog code $(obs[i].catalog) does not correspond to MPC catalog code. Setting debiasing corrections equal to zero."
            end
            α_corr[i] = 0.0
            δ_corr[i] = 0.0
            continue
        else
            # Otherwise, if star catalog is present in debias table, compute corrections
            # get pixel tile index, assuming iso-latitude rings indexing, which is the formatting in tiles.dat
            # substracting 1 from the returned value of `ang2pixRing` corresponds to 0-based indexing, as in tiles.dat
            # not substracting 1 from the returned value of `ang2pixRing` corresponds to 1-based indexing, as in Julia
            # since we use pix_ind to get the corresponding row number in bias.dat, it's not necessary to substract 1
            # @show α_i_rad δ_i_rad π/2-δ_i_rad
            pix_ind = ang2pixRing(resol, π/2-δ_i_rad, α_i_rad)
            # Healpix.pix2angRing(resol, pix_ind)
            # @show Healpix.pix2angRing(resol, pix_ind)
            # Handle edge case: in new MPC catalog nomenclature, "UCAC-5"->"Y"; but in debias tables "UCAC-5"->"W"
            if obs[i].catalog == "Y"
                cat_ind = findfirst(x->x=="W", mpc_catalog_codes_201X)
            else
                cat_ind = findfirst(x->x==obs[i].catalog, mpc_catalog_codes_201X)
            end
            # @show pix_ind, cat_ind
            # read dRA, pmRA, dDEC, pmDEC data from bias.dat
            # dRA, position correction in RA*cos(DEC) at epoch J2000.0 [arcsec];
            # dDEC, position correction in DEC at epoch J2000.0 [arcsec];
            # pmRA, proper motion correction in RA*cos(DEC) [mas/yr];
            # pmDEC, proper motion correction in DEC [mas/yr].
            dRA, dDEC, pmRA, pmDEC = bias_matrix[pix_ind, 4*cat_ind-3:4*cat_ind]
            # @show dRA, dDEC, pmRA, pmDEC
            # utc_i = DateTime(obs_t.yr[i], obs_t.month[i], obs_t.day[i]) + Microsecond( round(1e6*86400*obs_t.utc[i]) )
            et_secs_i = str2et(string(obs[i].date))
            tt_secs_i = et_secs_i - NEOs.ttmtdb(et_secs_i)
            yrs_J2000_tt = tt_secs_i/(daysec*yr)
            # Total debiasing correction in right ascension (arcsec)
            α_corr[i] = dRA + yrs_J2000_tt*pmRA/1000 
            # Total debiasing correction in declination (arcsec)
            δ_corr[i] = dDEC + yrs_J2000_tt*pmDEC/1000 
        end
    end

    return α_corr, δ_corr
end