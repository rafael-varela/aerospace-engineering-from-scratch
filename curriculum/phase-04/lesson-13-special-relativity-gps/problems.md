# P04.13  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: an airliner's clock)

An airliner cruises for 10 h at (200 + 2s) m/s at 11 km altitude. Treat the ground clock as at rest on a non-rotating sphere of radius R = 6378.137 km.

(a) Compute the time the airliner's clock loses to the velocity term, in ns, with `time_dilation_rate`.
(b) Compute the time it gains to the gravitational term, in ns, from (GM/c²)(1/R − 1/(R + h)), and compare with the local form gh/c² with g = 9.80665 m/s².
(c) Which term wins, and could a clock with a stability of 10⁻¹³ over 10 h see the net?

Answer form: three numbers in ns, one sentence.

## 2. (Derivation: the Lorentz transformation and what follows)

(a) Starting from x′ = A(x − vt) and x = A(x′ + vt′), and requiring that a light ray x = ct is x′ = ct′ in S′, show A = γ and derive t′ = γ(t − vx/c²).
(b) Show that c²Δt′² − Δx′² = c²Δt² − Δx².
(c) Apply a boost v₁ and then a boost v₂ to an event and show the result is a single boost with w = (v₁ + v₂)/(1 + v₁v₂/c²). Show that |w| < c whenever |v₁|, |v₂| < c.
(d) From f = f₀/(γ(1 − β cos θ)), derive the head-on, receding and transverse limits, and show the product of the head-on and receding frequencies is f₀².

Answer form: four derivations.

## 3. (Unfamiliar scenario: a navigation constellation at Mars)

Mars: GM = 4.2828 × 10¹³ m³/s², R = 3389.5 km. Ignore Mars's rotation and oblateness. A navigation satellite is on a circular orbit at altitude (300 + 10s) km.

(a) Compute the velocity, gravitational and net terms (µs/day) relative to a clock on the Martian surface, with `gps_clock_budget` and the right arguments.
(b) Find the crossover altitude at Mars with `crossover_radius`, and explain why its ratio to R is the same as at Earth.
(c) An areostationary orbit (period 88 775 s) is proposed instead. Compute its radius and the net clock term. Which orbit's clock must be slowed before launch, and which sped up?

Answer form: three numbers per orbit with units, one radius, two sentences.

## 4. (Numerical, uses your code: the periodic terms)

A GPS satellite has e = 0.01 + s/1000.

(a) Compute the amplitude of the eccentricity correction in ns and in metres of range, and its maximum rate 2√(GMa) e n/c² (n the mean motion). What fraction of the factory offset is that rate?
(b) For an equatorial receiver, tabulate the Sagnac term every 15° of position angle from the western to the eastern horizon (use the geometry of `main()`).
(c) For each of the two terms, say whether it grows with time if left uncorrected, and its worst-case range error.

Answer form: three numbers and a fraction, a 13-row table, two sentences.

## 5. (Derivation and check: why the factory offset uses −3GM/2a)

(a) For a Kepler orbit, use dt = (1 − e cos E)/n dE and r = a(1 − e cos E) to show that the time average of 1/r over one orbit is exactly 1/a.
(b) Using vis-viva, show that the time average of −(v²/2 + GM/r)/c² is −3GM/(2ac²) for any e, and that the remainder integrates to (4) of the lesson.
(c) Check (b) numerically for e = 0.1 by integrating the instantaneous rate over E with `scipy.integrate.quad`. Report the instantaneous rate at perigee and apogee and the mean, in µs/day relative to a clock far away.

Answer form: two derivations, three numbers.

## 6. (Judgment: which corrections does a service need?)

Three services use the GPS constellation: (i) a surveying service needing 1 m, (ii) a car navigation service needing 10 m, (iii) a maritime service needing 1 km. Clock parameters are uploaded to the satellites once a day. For each service and each term (velocity, gravity, geoid terms, eccentricity at e = 0.02, Sagnac), say whether it must be corrected, and justify with the time after which the omitted term alone exceeds the requirement, or its bounded maximum. Then say which single term you would never omit, and why.

Answer form: a 3 × 5 table of yes/no with numbers, one sentence.
