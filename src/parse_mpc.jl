@doc raw"""
    MPCRadec{T <: AbstractFloat} 

An optical measurement obtained from a MPC formatted file. See https://minorplanetcenter.net/iau/info/OpticalObs.html
for a detailed description of the format. 

# Fields 

- `number::String`: object's number. 
- `temp_desig::String`: temporary designation.
- `disc_astk::String`: discovery asterisk.
- `note_1::String`: note 1.
- `note_2::String`: note 2.
- `date::TDBEpoch{T}`: days since j2000 epoch.
- `α::T`: right ascension [rad].
- `δ::T`: declination [rad].
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
    date::TDBEpoch{T}
    α::T
    δ::T
    ref::String
    mag::String
    band::String
    catalog::String
    info_2::String
    obs_code::String
    # Inner constructor 
    function MPCRadec{T}(
        number::String, temp_desg::String, disc_astk::String, note_1::String, 
        note_2::String, date::TDBEpoch{T}, α::T, δ::T, ref::String, mag::String, 
        band::String, catalog::String, info_2::String, obs_code::String
    ) where {T <: AbstractFloat}
        new{T}(number, temp_desg, disc_astk, note_1, note_2, date, α, δ, ref, mag, band, 
               catalog, info_2, obs_code)
    end
end

# Outer constructor
function MPCRadec(
    number::String, temp_desg::String, disc_astk::String, note_1::String, 
    note_2::String, date::TDBEpoch{T}, α::T, δ::T, ref::String, mag::String, 
    band::String, catalog::String, info_2::String, obs_code::String
) where {T <: AbstractFloat} 
    MPCRadec{T}(number, temp_desg, disc_astk, note_1, note_2, date, α, δ, ref, mag, band,
                catalog, info_2, obs_code)
end

# Regular expressions to parse each of MPCRadec's fields
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
    parse_date(year::Int, month::Int, day::Int, utc::T) where {T <: AbstractFloat}

Returns the corresponding `TDBEpoch` object. 
"""
function parse_date(year::Int, month::Int, day::Int, utc::T) where {T <: AbstractFloat}
    # Hour + minutes (in hours)
    hr, min_ = intndec(24 * utc)
    # Minutes + seconds (in minutes)
    min, sec_ = intndec(60 * min_)
    # Seconds + fraction of seconds 
    sec, frac = intndec(60 * sec_)
    # Return TDBEpoch
    return from_utc(year, month, day, hr, min, sec, frac; scale = TDB)
end

@doc raw"""
    ra(hrs::Int, min::Int, sec::T) where {T <: AbstractFloat}

Returns the right ascension in radians. 
"""
function ra(hrs::Int, min::Int, sec::T) where {T <: AbstractFloat}
    # Convert hours minutes seconds to degrees
    return deg2rad((15 * hrs) + (0.25 * min) + (1//240 * sec))
end

@doc raw"""
    dec(sgn::String, deg::Int, min::Int, sec::T) where {T <: AbstractFloat}

Returns the declination in radians. 
"""
function dec(sgn::SubString{String}, deg::Int, min::Int, sec::T) where {T <: AbstractFloat}
    # Convert degrees minutes seconds to degrees
    if sgn == "+"
        return deg2rad(deg + (min / 60) + (sec / 3600))
    elseif sgn == "-"
        return -deg2rad(deg + (min / 60) + (sec / 3600))
    end
end

@doc raw"""
    MPCRadec(m::RegexMatch)

Converts `m` to a `MPCRadec` object. 
"""
function MPCRadec(m::RegexMatch)
    date = parse_date(
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
    file_to_MPCRadec(file::String)  

Returns a vector of `MPCRadec` objects, each one corresponding to a line in MPC formatted `file`.
"""
function file_to_MPCRadec(file::String)
    # Read lines of mpc formatted file 
    lines = readlines(file)
    # Apply regular expressions
    matches = match.(mpc_line_regex, lines)
    # Eliminate nothings
    filter!(!isnothing, matches)
    # Convert matches to MPCRadec
    radecs = MPCRadec.(matches)
    # Sort observations by date
    sort!(radecs, by = x -> x.date)
    
    return radecs
end

@doc raw"""
    text_to_MPCRadec(f::Function, text::String)

Returns the matches of `TO.mpc_line_regex` in `text` and filters them by `f`. 
"""
function text_to_MPCRadec(f::Function, text::String)
    
    # Vector of observations 
    radecs = Vector{MPCRadec}(undef, 0)
    # Iterate over the matches 
    for m in eachmatch(mpc_line_regex, text)
        if f(m)
            push!(radecs, MPCRadec(m))
        end
    end
    # Sort observations by date
    sort!(radecs, by = x -> x.date)
    
    return radecs
end

@doc raw"""
    j2000_days(obs::MPCRadec)    

Returns the days since J2000 epoch. 
"""
j2000_days(obs::MPCRadec) = value(j2000(obs.date))

# Print method for MPCRadec
# number(α, δ) at TDBEpoch
function show(io::IO, m::MPCRadec{T}) where {T <: AbstractFloat} 
    print(io, m.number, "(", m.α, ", ", m.δ, ") at ", m.date)
end

@doc raw"""
    get_raw_text(url::String = "https://minorplanetcenter.net/mpec/K20/K20YA9.html")

Returns the raw html text of webpage `url`.
"""
function get_raw_text(url::String = "https://minorplanetcenter.net/mpec/K20/K20YA9.html")
    resp = get(url)
    text = String(resp.body)
    return text
end

# MPC main page url 
const mpc_url = "https://minorplanetcenter.net"

# Regex for next circular url 
const next_circular_regex = r"<a href=\"(?P<next>.*)\"><img src=\"/iau/figs/RArrow.gif\""

@doc raw"""
    iterate_mpc_circulars(
        number::Int;
        url1::String = "https://minorplanetcenter.net/mpec/K20/K20YA9.html",
        url2::String = "https://minorplanetcenter.net/mpec/K21/K21JL0.html",
        max_iter::Int = 1_000
    )


Iterates MPC circulars from `url1` to `url2` and returns the matches of `TO.mpc_line_regex` 
filtered by f. 
"""
function iterate_mpc_circulars(
    f::Function; 
    url1::String = "https://minorplanetcenter.net/mpec/K20/K20YA9.html",
    url2::String = "https://minorplanetcenter.net/mpec/K21/K21JL0.html",
    max_iter::Int = 1_000
)

    # Vector of observations
    obs = Vector{MPCRadec}(undef, 0)
    
    # Number of urls checked 
    n = 0
    # First url 
    u = url1

    while n < max_iter
        n += 1
        # Raw html text of webpage u 
        text = get_raw_text(u)
        # Observations found in text 
        obs_ = text_to_MPCRadec(f, text)
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
    # Sort observations by date 
    sort!(obs, by = x -> x.date)

    return obs
end