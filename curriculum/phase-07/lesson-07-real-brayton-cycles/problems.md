# P07.07  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact, which are closed forms and which are numerical. Draw the station diagram and mark every heat and work arrow before writing any equation.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

Gas data: constant cp with γ = 1.4 and R = 287.05 J/(kg·K), or the NASA dry-air model (`AIR`), as each problem says. Component values are typical, illustrative values; the lesson's set is `REAL` (η_c = 0.85, η_t = 0.88, combustor loss 4 %), `IC` (intercooler to T₁, 2 % loss) and `REC` (ε = 0.85, losses 3 % cold, 4 % hot).

## 1. (Warm-up: the real simple cycle by hand)

A simple cycle has π = 8 + s/5, turbine inlet temperature 1400 K, inlet 288.15 K, compressor η_c = 0.85, turbine η_t = 0.88 and no pressure losses.

(a) With constant cp, compute T₂, T₇ (the turbine exit), w_c, w_t, w_net, q_in and η by hand, using equation (1) of the lesson. Compare with the ideal cycle at the same π.
(b) Check (a) with `cycle_performance(constant_cp_gas(1.4), CycleSpec(...))`, then repeat with NASA air.
(c) At what value of x = π^((γ−1)/γ) does the real cycle's work vanish? Compare with the ideal cycle's.

Answer form: two temperatures, five energies and an efficiency, with the ideal comparison; two efficiencies; a value and a comparison.

## 2. (Derivation: the recuperated and intercooled cycles)

(a) Derive the ideal recuperated efficiency η = 1 − (T₁/T₃)π^((γ−1)/γ), and show that it equals the simple cycle's at π = (T₃/T₁)^(γ/(2(γ−1))).
(b) Show that, for two compressors of the same isentropic efficiency with intercooling back to T₁ and no pressure loss, the total work is least when each has pressure ratio √π (constant cp).
(c) Derive the ideal intercooled-recuperated efficiency, equation (6) of the lesson.
(d) Starting from the first law for the cycle and the entropy generated in each component (heat from a reservoir at T_H, rejected to the atmosphere at T₀), derive w_net = q_in(1 − T₀/T_H) − T₀Σs_gen.

Answer form: four derivations.

## 3. (Unfamiliar scenario: should an APU have a recuperator?)

A ground APU (typical, illustrative) has π = 4, turbine inlet temperature 1250 K and the lesson's `REAL` components.

(a) Compute η and w_net without a recuperator, and with the lesson's recuperator at ε = 0.6, 0.75, 0.85 and 0.95 (its 3 % and 4 % losses each time). Why does w_net not change between the recuperated cases?
(b) By how much does the recuperator at ε = 0.85 cut the fuel (heat) needed per unit of shaft work?
(c) Give one reason the answer favours a recuperator on a ground APU but not necessarily on an aircraft APU.

Answer form: a table of five rows and a sentence; a percentage; a reason.

## 4. (Numerical, uses your code: the design study at your own temperature)

Take a turbine inlet temperature of 1300 + 5s K and the lesson's component set.

(a) For each of the four configurations, sweep π from 2 to 40 in steps of 0.5 with `sweep` and find the π of maximum efficiency and of maximum work, and the efficiency at π = 12.
(b) Which configuration and π would you choose for (i) a turboprop core that must be small and light, (ii) a ground generator set that runs continuously? Justify with your numbers.

Answer form: a four-row table; two choices with reasons.

## 5. (Numerical: where the work goes)

For the intercooled-recuperated design point of the lesson (π = 12, 1400 K, NASA air):

(a) Use `entropy_generation` to list s_gen and T₀s_gen (T₀ = 288.15 K) for each component, and rank them.
(b) Verify w_net = q_in(1 − T₀/T_H) − T₀Σs_gen with T_H = 1400 K.
(c) Why do the "heat-transfer" components (combustor, intercooler, exhaust, recuperator) generate more entropy in total than the turbomachines? Which of them could a better design reduce?

Answer form: a ranked list of seven; a check; a paragraph.

## 6. (Judgment: recuperators on high-pressure-ratio aero engines)

The core of a large turbofan runs at π ≈ 40.

(a) With the lesson's component set at 1400 K, compute η at π = 25, 30 and 40 for the simple, recuperated and intercooled-recuperated cycles. Explain why the intercooled-recuperated cycle still beats the simple cycle at π = 40 while the recuperated cycle does not.
(b) Double every pressure loss of the intercooled-recuperated cycle at π = 12. What happens to η and w_net?
(c) List three costs of heat exchangers on an aircraft that the model does not include, and say for each whether it argues for or against fitting them.

Answer form: a 3 × 3 table and an explanation; two numbers; a list.
