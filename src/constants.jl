# Module abbreviation
const TO = TaylorOrbits

# Earth/Sun mass ratio [dimensionless]
const μ_ES = 1/333_000

# Gauss gravitational constant [AU^(3/2) / day / M_S^(1/2)]
const k_gauss = 0.017_202_098_95

# Earth's position and velocity indices
const earthdofs = nbodyind(27, ea)
# Sun's position and velocity indices
const sundofs = nbodyind(27, su)
# Moon's position and velocity indices
const moondofs = nbodyind(27, mo)