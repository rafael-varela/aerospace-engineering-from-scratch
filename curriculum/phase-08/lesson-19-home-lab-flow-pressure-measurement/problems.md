# P08.19  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and its uncertainty, and label every number as predicted, measured, SYNTHETIC, typical or a manufacturer's figure.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

Water at 20 °C: ρ = 998.2 kg/m³, μ = 1.002 × 10⁻³ Pa·s, σ = 0.073 N/m (typical). Air: ρ from the room, μ = 1.81 × 10⁻⁵ Pa·s.

## 1. (Warm-up: your own bottle, predicted)

Your bottle has an inner diameter D = 80 + s mm (± 0.5 mm) and a hole d = 3.5 + s/20 mm (± 0.1 mm), 30 mm above the base.

(a) How much water does it hold up to 180 mm above the hole centre? Is that within the lab's limit?
(b) With C_d = 0.61 ± 0.03, predict the time to drain from 180 mm to 30 mm, with its budget term by term. Which term dominates?
(c) Compute the jet's Reynolds number at 180 mm and at 30 mm.

Answer form: a volume and a yes or no; a time ± u and four terms; two numbers.

## 2. (Derivation: Torricelli, exactly, and the tilted manometer)

(a) Derive the draining law √h = √h₀ − kt with k = C_d(a/A)√(g/2) from Bernoulli's equation and mass conservation. Then keep the surface velocity and show that k becomes k/√(1 − (C_d a/A)²). How large is the correction for the lesson's bottle?
(b) Explain why C_d = C_c C_v, and why a rounded entry raises C_d.
(c) For a U-tube lying on a board tilted at θ, show that the along-tube separation L of the menisci gives Δh = L sin θ, and derive the relative uncertainty of V including the angle's uncertainty u_θ.
(d) Show that the pressure needed to push water out through a hole of diameter d is of order 2σ/d to 4σ/d, and compute the corresponding head for the lesson's 4 mm hole.

Answer form: four derivations and two numbers.

## 3. (Numerical, uses your code: reduce the SYNTHETIC data, or your own)

(a) Fit C_d for each run of `SYNTHETIC_seed819_drain.csv` (or of your measurement), with u_h = 0.5 mm, u_t = 0.08 s, d = 4.0 ± 0.1 mm, D = 90.0 ± 0.5 mm. Give the budget and χ²/dof.
(b) Compute E_n of the mean against 0.61 ± 0.03 and of run 1's 200 → 30 mm time against the prediction.
(c) The generator's hidden truth is C_d = 0.62 with d = 4.05 mm and D = 90.4 mm. Which of your results "found" it? What should the report say?
(d) Reduce `SYNTHETIC_seed819_manometer.csv`: airspeed ± u for every row; the calibration check of the tilted tube against the upright one; which readings are below resolution.

Answer form: a table; two E_n values; a number and a sentence; a table, an E_n and a list. Label everything SYNTHETIC if you used the supplied files.

## 4. (Unfamiliar scenario: a wind-tunnel contraction and an aeroplane's pitot)

(a) A small open-circuit wind tunnel has a contraction of area ratio CR = 4 + s/5 between its settling chamber and its test section. A water manometer tilted at 15° connects a wall tap in the settling chamber to a wall tap in the test section and reads L = 40 + s mm along the tube. The room is at 101.0 kPa and 21 °C. Derive V_ts = √(2Δp/(ρ_a(1 − 1/CR²))) and compute the test-section speed with its uncertainty (u_L = 0.6 mm, u_θ = 0.5°, ρ_a to 1 %). How wrong would you be if you ignored the settling-chamber velocity?
(b) The light aeroplane's pitot-static system at 60 m/s TAS and 3000 m (ISA) is connected to a water U-tube for a ground check. What Δh would it read in flight? Why is a water U-tube a good ground-test instrument for the airspeed indicator but not something you would fly? (This is a paper exercise: testing a real aircraft's pitot-static system is done only by certified maintenance personnel with approved test equipment; do not connect anything to an aircraft's ports yourself.)

Answer form: a derivation, a speed ± u and a percentage; a height and two sentences.

## 5. (Judgment: reading a 3 m/s fan to 5 %)

The fan blows 3 m/s at the straw (ρ_a = 1.196 kg/m³); each meniscus is read to 0.3 mm and the board's angle to 0.5°.

(a) Compute u_V/V upright and tilted at 30°, 10° and 5.3°. Then at 5.3° and 10° with the angle measured to 0.1°.
(b) An alternative is to average many upright readings. How many would you need for 10 %, assuming independent readings? Why is the assumption doubtful?
(c) Which single change would you make to reach 5 %, and what would you have to measure better to make it work?

Answer form: six relative uncertainties; a number and a reason; a recommendation.

## 6. (Judgment: what does the Reynolds number tell you here?)

(a) For the 4 mm hole, Re_jet is 3000–8000. For a 1 mm hole in the same bottle it would be 760–2000. Do you expect the same C_d? What would you measure to find out?
(b) The straw is 6 mm across. At 3 m/s its Reynolds number is about 1200; a 1 mm tube would be about 200. Which would you trust as a pitot probe, and why?
(c) The fan's own Reynolds number is about 7 × 10⁴. Does that matter for your measurement? What property of the fan's jet does matter, and how did your procedure deal with it?

Answer form: an expectation and a measurement; a choice and a reason; a sentence and a property.
