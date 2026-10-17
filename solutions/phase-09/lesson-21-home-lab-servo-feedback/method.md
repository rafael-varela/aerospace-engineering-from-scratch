# P09.21  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
(a) K = 400 + 10s, τ = (30 + s)/1000 s, M_p = 5 + s/2 %; ζ = `zeta_from_overshoot(M_p)`; Kp = 1/(4τζ²K) (`design_p_position`); ω_n = √(KpK/τ). (b) `predict_closed_loop([K], [τ, 1, 0], Kp, 0, 20, 1.5)["metrics"]`, or `second_order_metrics(ζ, ω_n)`; compare `settling_time` with 4/(ζω_n). (c) u₀ = Kp × 20 + 0.09; largest linear step = (0.6 − 0.09)/Kp. (d) ±(0.10 − 0.09)/Kp.

## 2.
(a) K/(s²(τs + 1)) = K[1/s² − τ/s + τ²/(τs + 1)]; inverse transform; multiply by u; replace t by t − L for t ≥ L; the exponential term vanishes for large t, leaving Ku(t − L − τ); subtract. (b) Closed loop KpK/(τs² + s + KpK); divide by τ; read ω_n² = KpK/τ and 2ζω_n = 1/τ; invert M_p(ζ) and solve for Kp. (c) C G = Kp(s + 1/τ)/s × K/(τ(s + 1/τ)) = KKp/(τs); T = 1/(1 + τs/(KKp)); for 1 − e^{−t/τ_cl} solve for the 10 %, 90 % and 98 % times. (d) L(s) = ω_c/(s(θs + 1)); T = ω_c/(θs² + s + ω_c); read ζ, ω_n; at θω_c = 1/4, ζ = 1 and ω_n = 2ω_c; find the 10 % and 90 % roots x₁₀, x₉₀ of 1 − e^{−x}(1 + x) by bisection; compare (x₉₀ − x₁₀)/(2ω_c) with ln 9/ω_c.

## 3.
(a) τ = (15 + s/2)/1000; `design_p_position(1, τ, 5)`; `predict_closed_loop([1], [τ, 1, 0], Kp, 0, 10, 1)`. (b) θ = 1.5 Ts; ω_c from |L(jω)| = 1 with L = Kp/(jω(τjω + 1)) (a dense frequency grid); PM = 90° − atan(τω_c), minus ω_cθ in degrees; lag model: `predict_closed_loop` with den = conv([τ, 1, 0], [θ, 1]); check: `simulate_closed_loop(..., Ts, delay_samples=1)` then `step_metrics`. (c) Ts = 2.5 ms. (d) `design_p_sampled([1], [τ, 1, 0], 10, 0.01, 5, kp_lo, kp_hi, delay_samples=1)`; compare its rise time with (a).

## 4.
(a) `write_synthetic_csvs(d, s)`, `res = analyse_dataset(d)`; `res["motor"]["K"]`, `["tau"]`, `["dead_zone"]`, each `fits[i]["ls"]["p"]` and `["se"]`; `res["servo"]["K"]`, `["tau"]`, `["delay"]`. (b) For each run: `cmp_pred` (continuous against synthetic) and `cmp_sim` (sampled model against synthetic). (c) Loop over ten seeds; mean and standard deviation (ddof = 1) of both overshoots; compare the mean difference with the run-to-run standard deviation. (d) Label.

## 5.
(a) Identify K, τ, L from the shipped files (`analyse_dataset`); define rise(Kp) from `simulate_closed_loop([K], [τ, 1], Kp, Kp/τ, 20, 0.02, 1.2, u_min=-60, u_max=60, plant_delay=L, quant=0.066)`; bisect on Kp for 0.20 s; round Kp and Ki; re-simulate with the rounded gains. (b) `predict_closed_loop` with the rounded gains. (c) `design_pi_first_order(K, τ, 0.20)` and simulate those gains on the sampled model.

## 6.
Map each pattern to an effect: integrator against a dead band (hunting; test: Ki = 0 or an error dead band), delay (test: change the loop rate), saturation (test: the logged duty at its limit), asymmetric friction (test: steps and open-loop slopes in both directions), slew limit (test: the maximum rate against step size), dead zone (test: an open-loop duty ramp to the breakaway duty), supply sag (test: the pack voltage and a repeated open-loop step).
