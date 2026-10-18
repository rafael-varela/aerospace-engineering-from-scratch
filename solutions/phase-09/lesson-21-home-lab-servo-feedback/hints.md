# P09.21  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. `zeta_from_overshoot`, then (3) of the lesson turned round: KpK = 1/(4τζ²). The first-sample duty is the proportional term plus the compensation; the linear model holds while that stays below the limit.
2. (a) Partial fractions of K/(s²(τs + 1)), then shift by L. (b) Match τs² + s + KpK to τ(s² + 2ζω_n s + ω_n²). (c) Put the PI zero on the plant pole and look at what is left of the loop. (d) The loop ω_c/(s(θs + 1)) is the plant of (b) with τ → θ and K Kp → ω_c; the critically damped step is 1 − e^{−x}(1 + x) with x = ω_n t.
3. The continuous design is Problem 1 with K = 1. The delay is the one-sample write plus half a sample of hold; the phase it costs at ω_c is ω_c θ in radians. For (d), `design_p_sampled` with `delay_samples=1`.
4. `analyse_dataset` returns everything per run in `res["motor"]["runs"]` and `res["servo"]["runs"]`; for (c) loop over seeds and collect `run["meas"]["overshoot_pct"]` and `run["sim_metrics"]["overshoot_pct"]` for the first motor run.
5. Bisection on Kp: rise time falls as Kp rises. Call `simulate_closed_loop` with `ki = kp/tau`, `plant_delay=L`, `quant`, the ±60° limit, then `step_metrics`.
6. For each pattern ask what in the loop is missing from the linear continuous model (delay, dead zone, a limit, a slew rate, an integrator meeting a dead band, a changing supply), and which single change of condition would make that effect larger or smaller.
