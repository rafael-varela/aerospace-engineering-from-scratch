# P09.08  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe scripts that call it, given to the figures shown. Simulations are of the linear model, not of an aircraft.

## 1. Your auto-throttle gains (s = 15: K_p = 500 N per m/s, K_i = 50 N/m)

**Assumptions.** The lesson's plant: c = 39.27 N per m/s, lags 0.3 s and 1.0 s; unity feedback of true airspeed.

(a) **L(0) = 500/39.27 = 12.73.** Error after a 2 m/s command: 2/(1 + 12.73) = **0.1457 m/s** (7.28 %). Speed lost after a 200 N drag increase: 200/(39.27 + 500) = **0.371 m/s**.

(b) **Stable**: K_crit = 4989.7 N per m/s, about 10 times K_p.

(c) **Type 1**, **K_v = 50/39.27 = 1.273 s⁻¹**, ramp error 0.5/1.273 = **0.393 m/s**. Closed-loop poles **−3.503, −0.374 ± 0.475j, −0.118 s⁻¹**: stable.

**Check.** The final-value theorem on 200·G S/s returns 0.37087 m/s, equal to 200/(c + K_p).

**Common mistakes.** Using K_p itself instead of L(0) = K_p/c (the units differ); forgetting to check stability before quoting an error.

## 2. Errors from the loop structure

(a) V = G(CE + D) and E = R − V, so E + GCE = R − GD and **E = SR − GSD**, S = 1/(1 + CG). T + S = L/(1 + L) + 1/(1 + L) = **1**.

(b) With R = A/s^{q+1}: sE = A sⁿ D₀/(s^q (sⁿD₀ + K₀N₀)). As s → 0, D₀, N₀ → 1:
- n = 0: A/(s^q(1 + K₀)): q = 0 gives **A/(1 + K_p)**; q ≥ 1 diverges (**∞**).
- n = 1: A s/(s^q(s + K₀)): q = 0 gives **0**; q = 1 gives **A/K_v** (K_v = K₀); q = 2 gives **∞**.
- n = 2: A s²/(s^q(s² + K₀)): q = 0, 1 give **0**; q = 2 gives **A/K_a**.

(c) The error from a disturbance alone is e = −v = −GS d. With d = −F/s (a drag increase), e_ss = lim s·GS·F/s = G(0)F/(1 + C(0)G(0)). With G(0) = 1/c: **e_ss = F/(c + C(0))**. If C has an integrator, C(0) = ∞ and **e_ss = 0**.

(d) L = 1/s² under unity feedback: T = 1/(s² + 1), a step gives E = s/(s² + 1), and lim sE = lim s²/(s² + 1) = 0, but e(t) = cos t oscillates for ever. The theorem's condition fails: sE has poles at ±j, on the imaginary axis. Likewise the P loop with K_p = 6000 (Problem 4d).

**Common mistakes.** Applying the theorem without checking the poles of sE; putting the disturbance after the plant and getting the wrong transfer function.

## 3. A satellite's attitude loop is type 2

(a) **L = (K_d s + K_p)/(J s²) = (s + 0.1)/(10 s²)**: **type 2**, K_p = K_v = ∞, **K_a = 0.1/10 = 0.01 s⁻²**. Closed loop 10s² + s + 0.1: poles **−0.05 ± 0.0866j s⁻¹** (ω_n = 0.1 rad/s, ζ = 0.5).

(b) A parabola of acceleration 10⁻⁴ rad/s² (θ_c = ½αt², i.e. A = α in the table): error **α/K_a = 10⁻⁴/0.01 = 0.01 rad (0.57°)**.

(c) θ/τ_d = 1/(J s² + K_d s + K_p), so **θ_ss = τ_d/K_p = 10⁻⁵/0.1 = 10⁻⁴ rad (0.0057°)**. The type counts integrators in the *loop*, but the disturbance enters after the controller and before the plant's integrators. To reject it the integrator must be in the controller, upstream of where the disturbance enters.

(d) PID: closed loop 10s³ + s² + 0.1s + 0.002. Routh: K_d K_p = 0.1 > J K_i = 0.02, so **stable**, with poles −0.0377 ± 0.0820j and −0.0245 s⁻¹. The loop is **type 3**. θ/τ_d = s/(10s³ + s² + 0.1s + 0.002) → **0**. A 3000 s simulation gives 1.2 × 10⁻¹⁸ rad at the end, peaking at 1.01 × 10⁻⁴ rad on the way; the PD loop stays at 1.000 × 10⁻⁴ rad.

**Check.** The PD disturbance error 10⁻⁴ rad is the torque balance K_p θ = τ_d: a spring holding the disturbance.

**Common mistakes.** Concluding from "type 2" that every constant disturbance is rejected.

## 4. Gain sweeps and the theorem's limits

(a)

| K_p | step error | least ζ |
|---|---|---|
| 100 | 28.20 % | 1.000 |
| 400 | 8.94 % | 0.722 |
| 1000 | 3.78 % | 0.387 |
| 2000 | 1.93 % | 0.194 |
| 4000 | 0.97 % | 0.041 |

(b) K_p = 400:

| K_i | ramp error (0.5 m/s²) | least ζ | slowest real pole |
|---|---|---|---|
| 10 | 1.964 m/s | 0.726 | −0.024 s⁻¹ |
| 40 | 0.491 | 0.734 | −0.127 |
| 100 | 0.196 | 0.453 | −0.182 |
| 200 | 0.098 | 0.146 | −0.070 |

The ramp error falls as 1/K_i, but past K_i ≈ 100 the oscillatory pair loses its damping and the slowest real pole slows down again.

(c) Linear simulation, dt = 0.1, 0.05, 0.025 s: 0.51591, 0.50341, 0.49716 m/s against the theorem's 0.49091 m/s: differences **0.02500, 0.01250, 0.00625 m/s = A dt/2** exactly. The zero-order hold feeds a staircase that lags the ramp by half a step on average; the theorem is exact.

(d) 1/(1 + K/c) = 0.0065 ⇒ **K ≈ 6000 N per m/s**, above K_crit = 4990. `final_value` raises; the simulated error after 60 s is 19.6 m/s and growing.

**Common mistakes.** Blaming the theorem for the ramp offset; tuning K_i for ramp error alone.

## 5. A climb is a disturbance

(a) F = 1100 × 9.80665 × sin 2° = **376.5 N**. With no auto-throttle the aircraft loses F/c = **9.59 m/s**, with time constant m/c = **28.0 s** (about 84 s to settle within 5 %).

(b) Linear simulation:

| controller | final loss | peak loss (time) | within 0.1 m/s of final after |
|---|---|---|---|
| K_p = 400 | 0.857 m/s | 0.889 m/s (7.7 s) | 4.5 s |
| PI 400/40 | 0 | 0.782 m/s (5.9 s) | 22.7 s |

(c) In steady state the engine must supply the full **376.5 N** extra. The PI result assumes it can. If the engine is near full power (a climb often starts there), the throttle saturates, the speed decays as in (a), and the integrator keeps growing while the error persists: windup, and an overshoot when the climb ends.

**Check.** The P loop's final loss F/(c + K_p) = 376.5/439.27 = 0.857 m/s.

**Common mistakes.** Using sin γ with γ in degrees; reporting the PI loop's zero final error without its 0.78 m/s transient.

## 6. What the steady-state answer leaves out

(a) PI loop: \|S\| = **0.0102, 0.2295, 0.9411, 1.405, 1.279, 1.027** at 0.01, 0.1, 0.3, 0.627, 1 and 3 rad/s. Slow drag changes (periods of minutes) are rejected almost completely; gusts with periods of 5–20 s are made worse by up to 40 %; fast gusts are untouched. Steady-state error describes only the left end of this curve.

(b) 1100 × 0.5 = 550 N for the acceleration plus up to c × 10 = 393 N for the extra drag at +10 m/s: about **940 N** by the end of the ramp. If the engine cannot supply it, the throttle saturates, the error grows beyond the 0.49 m/s predicted, and the integrator winds up.

(c) Trust the low-frequency results whose ingredients are well known (the zero step and drag error of PI, the ordering of the gains, the 1/K_i law), and treat with caution anything that depends on c, on the lags and on the transient (the 18 % overshoot, the 23 s settling, the critical gain). In a real aircraft a speed change also changes lift and height through the phugoid, and the altitude loop's reaction changes the speed. The flight test checks the dynamic numbers; the model only predicts them.

**Common mistakes.** Treating a small steady-state error as proof of good disturbance rejection at all frequencies.

## Check-yourself answers

1. T = CG/(1 + CG), S = 1/(1 + CG), disturbance path G/(1 + CG); S + T = 1.
2. The number of free integrators in L. Type n gives zero error to inputs of polynomial order below n, a finite error to order n, and an infinite one above.
3. 1/(1 + 9) = **10 %**.
4. 0.5/2 = **0.25 m/s**.
5. All poles of sE(s) must have negative real parts (the closed loop stable, and no persistent oscillation). Ignoring it returns a finite number, for example 0.65 % for K_p = 6000, while the true error grows without bound.
6. The integrator keeps adding thrust while any error remains, so the only steady state is the one where the extra thrust equals the extra drag and the error is zero; a proportional controller needs a nonzero error to produce the extra thrust.
7. Raising the gain reduces the damping of the closed-loop poles until the loop oscillates and then goes unstable (here above 4990 N per m/s), so the error is limited by stability before it reaches zero.
