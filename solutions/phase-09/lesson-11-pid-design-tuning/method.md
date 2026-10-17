# P09.11  Method

The route to each answer, without the final numbers. Pair with `hints.md`; `worked-solutions.md` has the numbers for s = 15.

## 1. The filtered PID

C(s) = [(Kp Tf + Kd)s² + (Kp + Ki Tf)s + Ki]/[s(Tf s + 1)]. Evaluate with `num, den = pid_tf(...)`, `H = np.polyval(num, 1j*w)/np.polyval(den, 1j*w)`, then `abs` and `np.degrees(np.angle(...))`. High-frequency limit: Kp + Kd/Tf (the proportional and the saturated derivative add in phase). Compare with |Kp + Ki/(jω) + Kd jω| at ω = 100. The integral term dominates at 1 rad/s, so the phase approaches −90° there. The filter caps the gain, so the phase returns to 0° at 10⁴ rad/s. Without the filter the gain grows like Kd ω, and gyro noise at hundreds of rad/s reaches the servo amplified.

## 2. Derivations

(a) u = C_r r − C y with C_r = b Kp + Ki/s (+ D if the derivative acts on e) and C = Kp + Ki/s + D. y = P(C_r r − C y) ⇒ T = P C_r/(1 + PC). The closed-loop poles are the roots of den_P den_C + num_P num_C, and only C appears there. (b) T = (P C_r)/(1 + PC). As s → 0, C ~ Ki/s → ∞ and C_r/C → 1, so T → PC/(1 + PC) → 1. (c) Expand χ(s) = den_P(s)(Tf s² + s) + num_P(s)[(Kp Tf + Kd)s² + (Kp + Ki Tf)s + Ki]. At fixed s it equals χ₀(s) + Kp χ_p(s) + Ki χ_i(s) + Kd χ_d(s). The equations are Re χ(λ) = 0, Im χ(λ) = 0 and χ(μ) = 0. (d) With Tf = 0 and P = β/(s(s + a)): χ = s²(s + a) + β(Kd s² + Kp s + Ki) = s³ + (a + βKd)s² + βKp s + βKi. Match coefficients. (e) Steady state of (1) with constant δ: θ̇ = q = 0, α̇ = Z_α α + q = 0 ⇒ α = 0, q̇ = 0 ⇒ M_δ(δ + d) = 0 ⇒ δ = −d. A PD has δ = Kp(b r − θ) in steady state (the derivative of a constant is zero). With b = 1 the error is e = r − θ = −d/Kp.

## 3. Your servo

`p = dict(pitch_plant(), w_servo=10 + s)`, `num, den = pitch_plant_tf(p)`, `Ku, wu, Pu = ultimate_gain(num, den)`, `ziegler_nichols(Ku, Pu)['PID']`. `pid_place_three_poles(num, den, [-4+4j, -4-4j, -1.5], 0.02)`. Get the poles with `np.roots(closed_loop_tf(num, den, gains)[1])`. For each design and b, run `simulate_pitch_loop(p, gains, np.radians(10), 8.0, 1e-3, b=b)` and pass θ to `step_metrics`; take the peak elevator as `np.abs(sim['u']).max()`. Decide by peak elevator against 15° and by overshoot.

## 4. Roll

(a) φ = p/s and p(s − L_p) = L_δ δ_a. (b) The solve as in 2(c). The unplaced root sits near −1/Tf. (c) u = C_r r − C y and y = P u ⇒ u(1 + PC) = C_r r. Multiply through by den_P den_C to get polynomials. (d) `tf_step` on (num_P num_Cr, χ) for φ and on (num_Cr den_P, χ) for δ_a, each times 30° in radians. Compare the peak of |δ_a| with 20°.

## 5. Windup and Tt

Run the six configurations for a 10 s horizon and use `step_metrics`. Back-calculation with small Tt drains the integrator almost as fast as it fills, so it behaves like clamping. With large Tt it hardly drains, so it behaves like no anti-windup. In between, the overshoot rises smoothly with Tt. The default Tt = √(Ti Td) with Ti = Kp/Ki and Td = Kd/Kp.

## 6. Judgment

(a) GM: `ultimate_gain(np.polymul(num_P, num_C), np.polymul(den_P, den_C))[0]`, and in dB 20 log₁₀. PM: dense log sweep of L(jω), unwrapped phase, the first crossing of |L| = 1. (b) Instability in flight risks loss of control and structural loads. The usual alternative is a small-amplitude frequency sweep or a multisine injected at the servo command in closed loop. The broken-loop response is computed from the measured signals, and the margins are read off it without crossing the stability boundary, usually at altitude, with a safety pilot and a disengage switch. (c) Rank the omitted effects by how soon each would bind in a 20° pull-up: servo rate limit, airspeed loss and the phugoid, nonlinear lift near α = 14°. (d) The high-frequency gain is Kp + Kd/Tf. Tf of two samples is barely resolved by the backward-Euler filter, whose pole Tf/(Tf + dt) differs from e^{−dt/Tf}.

## 7. Back-calculation

(a) With u = δ_max and v = P + I + D: İ = Ki e + (δ_max − P − D − I)/Tt = Ki e − (I − I*)/Tt with I* = δ_max − P − D. (b) The fixed point is I = I* + Ki Tt e, approached with time constant Tt. Then v = P + D + I → δ_max + Ki Tt e: the unsaturated output parks a margin Ki Tt e above the limit. (c) Tt → 0: v → δ_max. Tt → ∞: the correction term vanishes. (d) Loop `pid_step` with constant (r, y), record I. Compare the last value with I* + Ki Tt e, which for y = 0 and D = 0 is δ_max − Kp b r + Ki Tt e. Fit log|I − I_final| against t. (e) In the closed loop the error and the P and D terms change on a 0.1–0.5 s time scale, the same as Tt, so the quasi-steady result is only a rough guide.
