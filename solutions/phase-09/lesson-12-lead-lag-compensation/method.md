# P09.12  Method

The route to each answer, without the final numbers. Pair with `hints.md`; `worked-solutions.md` has the numbers for s = 15.

## 1. A lead by hand

φ_max = arcsin((1 − α)/(1 + α)); gain 1/√α and 20 log₁₀(1/√α) dB; T = 1/(ω_m√α), zero 1/T = ω_m√α, pole 1/(αT) = ω_m/√α; high-frequency gain 1/α. Check with `lead_max_phase(alpha)` and `abs(freq_response(*lead_compensator(1, alpha, T), [40.0]))`. Centring at 400 rad/s scales the zero and pole by 10 and leaves φ_max and both gains unchanged.

## 2. Derivations

(a) T/(1 + ω²T²) = αT/(1 + α²ω²T²) ⇒ ω²T² = 1/α. At that point tan φ = (1 − α)/(2√α), and then sin φ = tan φ/√(1 + tan²φ) = (1 − α)/(1 + α). The gain is √((1 + 1/α)/(1 + α)) = 1/√α. (b) Evaluate the lag at s → 0 and s → ∞. Its phase at 10z is arctan 10 − arctan 10β, which decreases with β toward arctan 10 − 90°. (c) lim_{s→0} s L(s): the lead tends to Kc, the lag to β, and s G tends to K. (d) |L(jω)|² = ω_n⁴/(ω²(ω² + 4ζ²ω_n²)) = 1 ⇒ ω_c²/ω_n² = √(1 + 4ζ⁴) − 2ζ². arg L = −90° − arctan(ω_c/(2ζω_n)) ⇒ PM = arctan(2ζω_n/ω_c). Overshoot = e^{−πζ/√(1 − ζ²)}.

## 3. Heavier surface

Use `unwrapped_phase_deg(ng, dg, [30.])` for the plant phase. The lead comes from `design_lead`, and the lead-lag from `design_lead_lag(..., lag_ratio=10)`. Measure every design with `margins` and `velocity_constant`. For the steps, use `unity_feedback(*loop)` with `tf_response` and `step_metrics`; for the ramp error, take the 10 s value of t − y with the ramp held at its midpoint. Compare the results with the lesson's design. The extra phase lag of the heavier surface must come from the lead, so α falls. Kc and β also move, but by less.

## 4. Gimbal

(a) Read `margins(num, den)` with C = 1. By hand, the phase is −180° where 0.1ω · 0.02ω = 1, i.e. ω = √500; the gain margin is 1/|G| there. (b) Design with `design_lead(num, den, 15, 50)` and measure all four margins. The lead raises |L| at the phase crossover, so the GM falls. The second-order servo has no phase crossover, so its GM is infinite and cannot fall. (c) `design_lead(num, den, 40, 50)` raises `ValueError`. The phase that would be needed is in its message, and equals PM − (180° + arg G(j40)).

## 5. Judgment

(a) Loop `design_lead` over the five PMs and tabulate Kc/α. The motor-command noise is (Kc/α) × noise, so compare the PM 70° design with the PM 40° one. (b) The lag's pole at z/β is very slow. The closed loop keeps a pole near it and a zero near z, and the small residue of that pair produces a slow creep of the last few percent. Accept it when the ramp-following or disturbance requirement matters more than the 2 % settling time. (c) A PI raises the system type by one (zero ramp error for type 1 → type 2 plant), at the price of a pole at the origin and more phase lag at low frequency. Use it when a steady error to a ramp or to a constant disturbance must be zero, not just small. (d) Compute ω_c dt/2 in degrees and compare it with the PM, and compare dt/2 with the delay margin.

## 6. Two stages

Per-stage α from `lead_alpha(total/2)`, T = 1/(ω_c√α). Find Kc from |Kc C₁(jω_c)² G(jω_c)| = 1 and measure with `margins`. The step comes from `tf_response`. For (c), compute Kc/α² for two stages and Kc/α for one stage at 30 rad/s (with `design_lead(num, den, 30, 50)`), and repeat the two-stage design at 30 rad/s. (d) The lead's phase also raises the phase near the old phase crossover, which pushes the −180° crossing to a higher frequency. The lead's gain there is larger than at ω_c, so |L| at the new crossing is larger and the GM is smaller.
