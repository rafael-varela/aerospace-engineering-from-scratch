# P10.18  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact (closed form) and which are numerical. Unless a problem says otherwise, use the lesson's reference rotor and assumptions (κ = 1.15, C_d0 = 0.010, a = 5.73 /rad, θ_tw = −8°, f = 1.0 m², K = 4.65).

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: hover high and heavy)

The helicopter weighs 1000 + 10s kg.

(a) At sea level (ρ = 1.225 kg/m³) find v_h, the ideal power, the hover power from (5) and the figure of merit.
(b) Repeat at 2000 m (ρ = 1.0066 kg/m³, ISA). Which part of the power grows, which shrinks, and why?
(c) C_T/σ also grows with altitude. Why does that matter for the blades?

Answer form: four numbers per altitude; two sentences; one sentence.

## 2. (Derivation: the blade element and the edgewise disk)

(a) Derive dC_T = (σa/2)(θr² − λr)dr from the lift of a blade element with small angles, and integrate it with θ = θ₀ + θ_tw r and uniform λ to obtain (3).
(b) Show that with uniform inflow the blade element's induced power ∫λ dC_T equals the momentum-theory C_T^{3/2}/√2 in hover.
(c) Derive the profile-power factor K = 3 by averaging (r + μ sin ψ)³ over the disk. Why does counting the radial flow raise it?
(d) Show that Glauert's inflow (6) tends to λ_h as μ → 0 and to C_T/(2μ) when μ ≫ λ, and give the dimensional form of the second limit.

Answer form: four derivations.

## 3. (Unfamiliar scenario: a quieter rotor)

Rotor noise rises steeply with tip speed, so a noise-sensitive operator asks to slow the rotor to 195 m/s or 180 m/s at the same mass.

(a) For each tip speed, compute the hover power, its induced and profile parts, the figure of merit, C_T/σ, the mean blade lift coefficient 6C_T/σ and the pitch at 0.75R.
(b) Explain why the induced power does not change and the profile power falls.
(c) What limits how far the tip speed can be reduced? Compute the advance ratio at 60 m/s for each tip speed and comment on the retreating blade.
(d) Compute the minimum forward-flight power for each tip speed.

Answer form: a table; two sentences; three numbers and a paragraph; three numbers.

## 4. (Numerical, uses your code: the power curve of your helicopter)

Take f = 0.8 + s/50 m².

(a) Find the minimum-power speed and power, the maximum-range speed and power, and the autorotation estimates (least descent rate and best glide ratio).
(b) The engine and transmission deliver 160 kW at 90 % efficiency to the rotor. Find the largest rate of climb (at which speed?) and the maximum level speed. Check the climb result with `forward_flight_power(V, rotor, climb=V_c)`.
(c) Plot the power curve and draw the tangent from the origin. Explain in one sentence why the maximum-range speed is the tangent point.

Answer form: five numbers; three numbers and a check; a plot and a sentence.

## 5. (Numerical, uses your code: the profile-power factor)

(a) Compute K(μ) = (I(μ) − 1)/μ² from `profile_power_integral` with and without radial flow for μ = 0.05, 0.1, 0.2, 0.3, 0.4 on a 400-point grid.
(b) Without radial flow K is slightly above 3. Where does the excess come from? (Hint: compare |u_T|³ with u_T³.)
(c) Recompute the minimum-power and maximum-range speeds and the power at 60 m/s with K = 3, 4.65 and 5. Where on the curve does K matter?

Answer form: a table; a sentence; a table and a sentence.

## 6. (Judgment: where does the model fail?)

(a) Vary κ (1.10, 1.15, 1.25) and C_d0 (0.008, 0.010, 0.012). Tabulate the hover power, FM and the least autorotative descent rate. Which input should a flight-test programme pin down first, and for which question?
(b) At 40, 60 and 70 m/s compute μ, the radius inside which the retreating blade is in reverse flow, the retreating-tip speed and the advancing-tip Mach number (a = 340.3 m/s). Which phenomena does the model leave out at 70 m/s?
(c) Why can this model not predict a vertical autorotation (descent at zero forward speed)? What happens to the flow through the rotor?
(d) The energy estimate gives a least descent rate of 6.9 m/s. List two effects that make the real value larger.

Answer form: a table and a sentence; a table and two phenomena; a paragraph; two effects.
