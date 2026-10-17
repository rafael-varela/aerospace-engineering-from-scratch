# P07.02  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact (closed form), which are numerical, and which are typical, illustrative inputs. Everything here is computed: **no pressure vessel is built, filled or vented at home in this course.**

Several problems use your personal seed **s**: the day of the month you were born (1–31). Use γ = 1.4 and M = 28.0134 g/mol for N₂, γ = 5/3 and M = 4.002602 g/mol for He, R_U = 8.314 462 618 J/(mol·K).

## 1. (Warm-up: closed systems)

(a) One kilogram of N₂ at 1 MPa and T₁ = 250 + 5s K expands to three times its volume along p vⁿ = constant with n = 1, 1.2 and 1.4. For each path find p₂, T₂, w, Δu and q. Which quantities are the same for every path, and why?
(b) A sealed 10 L helium bottle at 20 MPa and 300 K sits in the sun and warms by s kelvin. Find the mass, the final pressure and the heat absorbed. Show that the heat is p₀V ΔT/((γ − 1)T₀).

Answer form: a 3-row table and a sentence; three numbers and a one-line proof.

## 2. (Derivation: the blowdown from first principles)

(a) Starting from the control-volume energy equation, derive m cv dT/dt = −ṁ R T + Q̇ for a rigid, well-mixed tank with one exit. State each assumption where you use it.
(b) From h₀ = h + V²/2 and the isentropic relations, derive the choked mass flow ṁ = A Γ p₀ √(γ/(R T₀)) with Γ = (2/(γ + 1))^((γ+1)/(2(γ−1))), and the critical pressure ratio.
(c) For adiabatic choked flow, show that the remaining gas follows p ∝ ρ^γ and derive p(t)/p₀ = (1 + (γ − 1)t/(2τ))^(−2γ/(γ−1)) with τ = V/(Cd A Γ a₀). Show that τ = m₀/ṁ₀.
(d) Explain why `energy_residual` is at round-off for any Runge–Kutta step size, while the agreement with (c) is not.

Answer form: four derivations.

## 3. (Unfamiliar scenario: filling a tank)

The opposite of a blowdown: an evacuated rigid tank is filled adiabatically from a supply line at constant T_line (16.Unified notes, 2.3.3).

(a) Show that, filled from vacuum to any pressure, the gas in the tank ends at T = γ T_line. Evaluate for N₂ and He from a 300 K line.
(b) A 50 L tank is filled from vacuum with N₂ to 30 MPa this way. Find the mass, and the pressure after the tank cools back to 300 K. How much more mass would an isothermal fill to 30 MPa hold?
(c) The same tank, starting at 10 MPa and 300 K, is topped up adiabatically to 30 MPa from the 300 K line. Find the final temperature and mass. (Hint: m₂u₂ − m₁u₁ = h_line(m₂ − m₁).)
(d) Why do real filling procedures fill slowly or in stages?

Answer form: two temperatures, three numbers, two numbers, two sentences.

## 4. (Numerical, uses your code: your orifice)

Use the engineering-question tank (50 L, 30 MPa, 300 K, N₂, back pressure 0.5 MPa, Cd = 0.8) with an orifice of diameter d = 1.5 + 0.05s mm.

(a) Compute ṁ₀, τ and the closed-form adiabatic and isothermal times to 3 MPa.
(b) Run `simulate_blowdown` in both modes with `p_stop=3e6`. Report the event times, final temperatures, stranded masses and energy residuals, and compare with (a).
(c) Compute the recovery pressure after the adiabatic case warms back to 300 K.
(d) By what factor did the orifice change τ relative to d = 2 mm? Is that the factor you expected?

Answer form: four numbers, a table, one number, a ratio and a sentence.

## 5. (Numerical, uses your code: which pressurant, and how much does the wall help?)

For the engineering-question tank and 2 mm orifice, run nitrogen and helium with hA = 0, 10 and 100 W/K to 3 MPa.

(a) Tabulate time, final temperature and stranded mass for the six cases.
(b) Compute the initial thermal time constant m₀ cv/(hA) for each gas and conductance, and find the hA at which it equals the blowdown time constant τ.
(c) Using (b), explain why helium needs a larger hA than nitrogen to approach its isothermal limit.
(d) Helium's stored mass is a seventh of nitrogen's. Why do launch vehicles nevertheless prefer helium as a pressurant? (Think about what the stage has to carry, and about condensation in cryogenic propellant tanks.)

Answer form: a table, six time constants and two crossover values, two paragraphs.

## 6. (Judgment: is this model good enough?)

For each case say which assumption of the lesson's model fails first, what the consequence is, and what you would do instead. Two sentences each.

(a) Sizing a 30 MPa COPV's nitrogen load to ±1 % for a flight vehicle.
(b) Predicting the gas temperature at the outlet of a bottle during a 5 s emergency vent.
(c) A slow, hour-long draw of helium from a bottle mounted on an aluminium structure.
(d) The last 20 % of an adiabatic nitrogen blowdown to a 0.5 MPa back pressure.
(e) Using the adiabatic result as a conservative (worst-case) temperature for a seal specification.

Answer form: five short paragraphs.
