# P04.09  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a charged particle in the Earth's field)

A proton moves at v = (1 + s/10) × 10⁶ m/s perpendicular to a uniform field of 30 μT (m_p = 1.6726 × 10⁻²⁷ kg, e = 1.6022 × 10⁻¹⁹ C).

(a) Find the magnitude of the force, the gyroradius and the period of the circular motion. Check one result with `gyroradius` and `cyclotron_frequency`.
(b) Repeat for an electron (m_e = 9.109 × 10⁻³¹ kg) at the same speed. Which way does each particle circle when B points out of the page?
(c) A 1 m length of wire carrying 1 A lies across the same field. What force does it feel, and in which direction relative to the wire and the field?

Answer form: six numbers with units, two directions, one number and a direction.

## 2. (Derivation: the loop axis field, the solenoid and the torque on a loop)

(a) From the Biot–Savart law (2), derive the field (3) on the axis of a circular loop. Show which components cancel and why. Then show that for |z| ≫ R it becomes μ₀m/(2π|z|³) with m = IπR², and find the first correction in R/z.
(b) Integrate (3) over the windings of a solenoid of n turns per metre from −L/2 to L/2 to obtain (4). Show that the centre field tends to μ₀nI and the end field to half of it as L/R → ∞.
(c) A rectangular loop of sides a and b, current I, can rotate about an axis through its centre parallel to the sides of length a; that axis is perpendicular to a uniform B, and the loop's normal n̂ makes angle θ with B. Compute the force on each of the four sides, show the net force is zero, and show the torque is IabB sin θ in the direction of n̂ × B. Conclude τ = m × B with m = IA n̂.

Answer form: three derivations, one correction term.

## 3. (Unfamiliar scenario: an electrodynamic tether)

A satellite in a 5600 s circular orbit trails a straight conducting tether of length L = (500 + 10s) m, perpendicular both to its velocity 7600 m/s and to the local field of 30 μT. The circuit closes through the ionosphere; take its total resistance as 100 Ω.

(a) Compute the motional emf, the current and the electrical power.
(b) Compute the force on the tether and show by Lenz's law that it opposes the satellite's motion.
(c) For a 100 kg satellite, what is the deceleration and the speed lost per orbit? Is this a useful deorbit device?
(d) Show that the mechanical power lost (force × speed) equals the electrical power. What would change if a battery drove current the other way?

Answer form: three numbers, a force with a direction argument, two numbers and a sentence, an identity and a sentence.

## 4. (Numerical, uses your code: size a torquer for a 1U CubeSat)

A 1U CubeSat (I = 0.002 kg·m² about each axis) in a 400 km polar orbit needs a coil moment of m_req = 0.1 + s/100 A·m² from a square coil of side 9 cm. The coil current may not dissipate more than 0.5 W.

(a) With `field_along_orbit`, find the minimum, maximum and mean |B| around the orbit, and the largest and smallest torque the coil can give with m ⊥ B.
(b) With the `Torquer` dataclass, tabulate current, resistance, power, inductance estimate and L/R for N = 100, 150, 200, 300 and wire diameters 0.20, 0.25, 0.30 mm, each set to deliver m_req. Pick a design that meets the power limit and justify the choice in one sentence (mass, voltage, time constant).
(c) With `bdot_detumble` at the mean field, k = 2 × 10⁴ A·m²·s/T and m_max = m_req, detumble from 5°/s. Is the coil saturated at the start? Report T_d, the time to 0.1°/s, and the lower bound Iω₀/(m_max B).

Answer form: three fields and two torques, a 12-row table and a choice, three numbers and a yes/no.

## 5. (Numerical and derivation: a Helmholtz pair for magnetometer calibration)

Two coaxial loops of radius R = 0.1 m, 1 A each in the same sense, are separated by a distance d along their common axis.

(a) Using `loop_points(R, n, center=...)` and `wire_field`, compute B_z on the axis at the midpoint and at z = 0.1R and 0.2R from it, for d = 0.8R, R and 1.2R. Report the relative deviation from the midpoint value.
(b) Show from (3) that d = R makes d²B_z/dz² vanish at the midpoint, so the first non-zero deviation is fourth order in z, and find its coefficient. Compare with (a).
(c) Find the field at the midpoint for d = R in closed form and compare with your sum. How would you use the pair to test a satellite magnetometer?

Answer form: a 3 × 3 table, a derivation with a coefficient, a number with units and two sentences.

## 6. (Judgment: where a magnetorquer is the wrong tool)

Answer each in two or three sentences with a number where one is asked for.

(a) A team proposes magnetorquers as the only attitude actuator for a geostationary satellite. Estimate the field at GEO with the dipole model and give your verdict.
(b) A magnetometer is mounted 10 cm from the axis of the 0.154 A·m² torquer coil of the lesson. Estimate, with the dipole far-field formula, the field the coil makes at the magnetometer while it is on, compare with the Earth's field, and say what flight software must do about it.
(c) Explain why a satellite with three orthogonal torquers still cannot produce an arbitrary torque at a given instant, and why it can over an orbit.
(d) The B-dot law needs dB/dt in the body frame. Name two things other than the tumble that make the measured B change, and say which one limits how low the residual rate can go.

Answer form: four short paragraphs, two numbers.
