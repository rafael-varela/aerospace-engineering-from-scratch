# P04.05  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact (closed form) and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31). Use μ = 3.986004418 × 10¹⁴ m³/s² and R = 6371 km for the Earth unless a problem says otherwise.

## 1. (Warm-up: a circular orbit)

A small Earth-observation satellite flies a circular orbit at altitude h = (500 + 10s) km.

(a) Compute the local gravitational acceleration g(h), the orbital speed, the period in minutes and the specific orbital energy.
(b) Compute the escape speed at that radius and the Δv needed to escape from the circular orbit with one burn along the velocity.
(c) Check (a) two ways: v²/r must equal g(h), and the kinetic energy per unit mass must equal half the magnitude of the potential energy per unit mass.

Answer form: five numbers with units, two checks.

## 2. (Derivation: vis-viva, eccentricity, the third law, and the leapfrog's angular momentum)

(a) Using only conservation of energy and of angular momentum between periapsis r_p and apoapsis r_a, show that ε = −μ/(r_p + r_a) = −μ/2a and hence v² = μ(2/r − 1/a).
(b) Show that h² = μa(1 − e²) and invert it, with (a), to e = √(1 + 2εh²/μ²). What does e = 1 correspond to?
(c) Kepler's second law says the areal velocity is h/2. The area of an ellipse is πab with b = a√(1 − e²). Derive T = 2π√(a³/μ) for an elliptical orbit.
(d) For a central force a(r) = f(|r|) r, show that both the kick v ← v + (h/2) a(r) and the drift r ← r + h v leave r × v unchanged. Conclude that the leapfrog conserves angular momentum exactly. Does the same argument work for RK4? Why not?

Answer form: four short derivations and one sentence.

## 3. (Unfamiliar scenario: Mars)

For Mars, μ = 4.282837 × 10¹³ m³/s², mean radius 3389.5 km, sidereal rotation period 88 642.66 s.

(a) Find the radius, altitude and speed of an areostationary orbit (the Martian equivalent of geostationary).
(b) Find the surface gravity and the escape speed from the surface. Compare with the Earth's 9.82 m/s² and 11.2 km/s.
(c) The moon Phobos orbits at a = 9376 km. Compute its period. Is it shorter or longer than a Martian day, and what does that mean for the direction in which an observer on the surface sees Phobos rise?

Answer form: six numbers with units, one sentence.

## 4. (Numerical, uses your code: step size for an energy tolerance)

Take a = 8000 km and e = 0.05 + s/200.

(a) Give the perigee and apogee altitudes and the perigee and apogee speeds (vis-viva).
(b) Propagate ten orbits with `propagate_leapfrog` at 250, 500, 1000, 2000 and 4000 steps per orbit. Tabulate the one-orbit position error against `propagate_solve_ivp`, the maximum relative energy error over ten orbits, and the maximum relative angular-momentum error. Report the observed order of the position error.
(c) From the energy error at 1000 steps per orbit and the h² law, predict the number of steps per orbit for a maximum relative energy error of 10⁻⁶ over ten orbits. Run it and check.
(d) With 2000 steps per orbit over five orbits, measure the period with `period_from_trajectory` and compare with Kepler. Is the simulated orbit early or late?

Answer form: four numbers, a 5-row table, one order, one step count with its check, one period difference.

## 5. (Numerical: a transfer to a navigation-satellite orbit)

GPS satellites orbit at a radius of 26 560 km. A spacecraft starts in a circular parking orbit at (300 + s) km altitude.

(a) Compute Δv₁, Δv₂, the total and the transfer time of a coplanar Hohmann transfer, and the eccentricity of the transfer ellipse.
(b) Check Δv₁ independently: put the post-burn state into `orbital_elements` and confirm the apoapsis radius is 26 560 km.
(c) Real GPS orbits are inclined 55° to the equator. Without computing it, explain why the real transfer from a parking orbit of a different inclination costs more, and at which point of the transfer the plane change is cheapest.

Answer form: five numbers, one check, two sentences.

## 6. (Judgment: when the two-body model stops being enough)

For each case, say whether the point-mass two-body model with the leapfrog at 1000 steps per orbit is adequate, what fails first if not, and what you would add or change. One or two sentences each.

(a) Predicting where the ISS (400 km) will be in 24 hours, to within 1 km.
(b) Estimating the Δv of a LEO-to-GEO transfer for a feasibility study.
(c) Propagating a geostationary satellite for ten years to plan station-keeping.
(d) Checking that a comet on an e = 0.99 orbit returns at the predicted time.
(e) Computing the period of a satellite of 1000 kg against one of 1 kg in the same orbit.
(f) Long-term statistics of a debris cloud over 10⁵ orbits, where only the energy distribution matters.

Answer form: six short paragraphs.
