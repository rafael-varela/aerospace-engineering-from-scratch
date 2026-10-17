# P09.10  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical. Sketch every Nyquist plot by hand, including the image of any indentation, before computing it.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: your servo gain)

Take the lesson's servo with K = 10 + 2s s⁻¹.

(a) Give the gain margin (factor and dB) from formula (2) and the phase crossover frequency.
(b) Find the phase margin, the gain crossover and the delay margin with `margins`.
(c) Count N, P and Z with `nyquist_count`, and give the peak sensitivity M_s = max 1/\|1 + L(jω)\|.
(d) Does the design meet "GM ≥ 6 dB and PM ≥ 45°"?

Answer form: two numbers and a frequency; three numbers; three integers and a number; yes or no.

## 2. (Derivation: the criterion and the servo's margin)

(a) From the argument principle, derive Z = N + P for unity negative feedback, stating the direction of the contour and of N.
(b) Show that on the indentation around a k-fold pole at s = 0, L(s) traces an arc of angle kπ at large radius, clockwise. Sketch it for k = 1 and k = 2.
(c) For L = K/(s(τ_m s + 1)(τ_f s + 1)), derive ω_pc = 1/√(τ_mτ_f) and GM = (τ_m + τ_f)/(Kτ_mτ_f), and show it equals the Routh limit.
(d) Show that the delay margin is PM/ω_gc when there is one gain crossover. What changes when there are several?

Answer form: four derivations, two sketches.

## 3. (Unfamiliar scenario: an unstable airframe with a computer)

The lesson's unstable airframe at K = 20, L = 20(s + 3)/((s − 2)(s + 10)), is flown by a digital computer that adds a delay.

(a) Find its phase margin, gain crossover and delay margin.
(b) Count Z for delays of 0.5, 0.9, 1.1 and 1.5 times the delay margin, by Nyquist with the exact delay and by an order-8 Padé closed loop.
(c) How does N change as the delay passes the delay margin? Explain why the loop needs N = −1.
(d) Why is the airframe's gain margin "infinite" in the usual sense but its gain margin downwards only −9.5 dB? Which one matters?

Answer form: three numbers; a table; a sentence; a short argument.

## 4. (Numerical, uses your code: the servo sweep)

(a) For K = 5, 10, 20, 40, 60, 100 s⁻¹ tabulate GM (dB), PM, ω_gc, DM and M_s.
(b) Plot PM and DM against K on log axes. Which falls faster, and why?
(c) Find, to two figures, the largest K that meets GM ≥ 6 dB and PM ≥ 45°. Which requirement binds?

Answer form: a table, a plot with one sentence, a gain and a requirement.

## 5. (The satellite: where to put the crossover)

For the satellite loop with K_p = 0.1, 0.3162, 1 and 3 N m/rad:

(a) Tabulate PM, ω_gc and DM.
(b) Show analytically that the phase margin is largest when ω_gc = 1/√(T_dτ), and give its value.
(c) The attitude computer runs at 2 Hz and takes one full period to compute (a 0.5 s delay). Which K_p values keep at least 30° of phase margin with that delay?

Answer form: a table, a derivation with one number, a list.

## 6. (Judgment: what the margins miss)

(a) The servo's filter lag doubles to τ_f = 0.02 s. Recompute K_crit and, at K = 20, GM, PM and DM.
(b) Add a 20 ms computation delay to the original servo at K = 20. Recompute GM, PM and ω_pc. Why does a pure delay reduce the gain margin, which is about gain?
(c) The servo at K = 20 has GM = 15.6 dB, PM = 43.2° and M_s = 1.77. What does M_s tell you that the two margins do not?
(d) Name one physical effect in a real control-surface actuator that none of these linear margins can capture.

Answer form: four numbers; three numbers and a sentence; two sentences; one sentence.
