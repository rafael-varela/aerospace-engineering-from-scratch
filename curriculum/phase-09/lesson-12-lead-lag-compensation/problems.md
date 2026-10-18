# P09.12  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact, which come from your code, and which are simulations of a representative model.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a lead by hand)

A lead (Ts + 1)/(αTs + 1) has α = 0.1 + s/100 and must give its maximum phase at 40 rad/s.

(a) Compute φ_max, the gain at 40 rad/s (ratio and dB), the zero and pole frequencies, and the high-frequency gain. (b) Check φ_max and the gain with `lead_max_phase` and with `freq_response` at 40 rad/s. (c) If the same lead were centred at 400 rad/s instead, which of your answers would change?

Answer form: six numbers, two checks, one sentence.

## 2. (Derivation: lead, lag and the margin-damping rule)

(a) Differentiate φ(ω) = arctan(ωT) − arctan(αωT) and show that the maximum is at ω = 1/(T√α), with sin φ_max = (1 − α)/(1 + α) and gain 1/√α there. (b) Show that the lag (s + z)/(s + z/β) has DC gain β and high-frequency gain 1, and that its phase at ω = 10z lies between −5.71° and 0° for every β > 1. (c) Show that the lead-lag loop Kc C_lead,1 C_lag K/(s(τs + 1)) has Kv = Kc β K. (d) For L = ω_n²/(s(s + 2ζω_n)), find the crossover frequency and show PM = arctan[2ζ/√(√(1 + 4ζ⁴) − 2ζ²)]. Tabulate PM and the closed-loop overshoot for ζ = 0.3, 0.5, 0.55, 0.7 and compare with PM ≈ 100ζ.

Answer form: four derivations and a four-row table.

## 3. (Numerical, uses your code: your heavier surface)

The servo of the lesson drives a heavier control surface: τ = 0.05 + s/1000 s, K = 20 s⁻¹. Targets: ω_c = 30 rad/s, PM = 55°, Kv ≥ 200 s⁻¹.

(a) Give the plant phase at 30 rad/s and the lead alone (Kc, α, zero, pole, φ_max); measure its crossover, PM, GM and Kv with `margins`. (b) Design the lead-lag with z = ω_c/10 and give β, the lag pole, the new α and Kc; measure the margins and Kv. (c) Tabulate overshoot, rise time, settling time and ramp error for both designs. (d) Compare with the lesson's τ = 0.05 s: which parameters moved most, and why?

Answer form: one phase and five design numbers with four measured ones; five design numbers with measured ones; an eight-entry table; two sentences.

## 4. (Unfamiliar scenario: a camera gimbal with a motor pole)

A UAV camera-gimbal axis has G = 10/(s(0.1s + 1)(0.02s + 1)): a mechanical lag and an electrical (or amplifier) lag.

(a) With C = 1, measure ω_c, PM, the phase crossover and GM. Check GM by hand (find where the phase is −180°). (b) Design a lead for ω_c = 15 rad/s and PM = 50°; measure all four margins of the compensated loop, and the closed-loop overshoot and settling. What happened to the gain margin, and why could that not happen in the lesson's servo? (c) Try ω_c = 40 rad/s with the same PM. What does `design_lead` say, and how much phase would one stage have to supply?

Answer form: four numbers and a hand check; four margins, two step metrics and two sentences; one message and one number.

## 5. (Judgment: what lead and lag cost)

(a) For lead designs on the lesson's servo at ω_c = 30 rad/s with PM = 40°, 50°, 55°, 60° and 70°, tabulate α, Kc, Kv and the high-frequency gain Kc/α. A potentiometer feedback signal carries 2 mV rms of noise per 1 V full scale. Which design passes the most noise to the motor command, by what factor relative to PM = 40°, and why? (b) The lead-lag of the lesson settles in 471 ms against 147.5 ms for the lead alone. Explain the slow tail in terms of the lag's pole-zero pair, and say when you would accept it. (c) When is a PI (a lag with β → ∞) better than a lag, and what does it change in the system type? (d) The servo is to be implemented on a microcontroller at 1 kHz. Compare the phase cost of the half-sample hold delay at 30 rad/s with the design's phase and delay margins, and say whether the continuous design can be used.

Answer form: a five-row table with one factor and one sentence; three short paragraphs.

## 6. (Hard: two lead stages)

For the gimbal of Problem 4 the target is ω_c = 40 rad/s with PM = 50°.

(a) Compute the total lead needed and design two identical cascaded stages, each supplying half, both centred at 40 rad/s: give α, the zero and pole, and the common gain Kc. (b) Measure the crossover, PM, phase crossover and GM of the compensated loop. Report the closed-loop overshoot, rise time and settling time. (c) The high-frequency gain of two stages is Kc/α². Compare it with what a single stage would need at ω_c = 30 rad/s (which one stage can just do) and with two stages there. Is "use two stages" a free lunch? (d) Explain why the gain margin of the two-stage loop is smaller than that of the uncompensated loop, using the phase crossover frequency.

Answer form: one total and four design numbers; four margins and three step metrics; a three-entry comparison and one sentence; two sentences.
