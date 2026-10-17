# P09.19  Worked solutions

Worked for **s = 15**. Every number below is produced by `solution.py` (the functions named in each problem) or by `scipy.linalg.solve_continuous_lyapunov`; simulations use `scipy.integrate.solve_ivp`. Units: SI, angles in radians unless stated.

## 1. Lyapunov by hand

**Assumptions.** k = 4 + 15/5 = **7 s⁻²**, c = 1 s⁻¹, Q = I.

**Derivation.** (1,1): −2kp₁₂ + 1 = 0 ⇒ p₁₂ = 1/14 = **0.071429**. (2,2): 2p₁₂ − 2cp₂₂ + 1 = 0 ⇒ p₂₂ = (1 + 1/7)/2 = **0.571429**. (1,2): p₁₁ = cp₁₂ + kp₂₂ = 0.071429 + 4.0 = **4.071429**.

**(b)** Eigenvalues of P: 0.56997 and 4.07289, both positive. Along trajectories V̇ = xᵀ(AᵀP + PA)x = −xᵀx < 0.

**(c)** `lyapunov_solve` and `solve_continuous_lyapunov(A.T, -I)` both return [[4.07143, 0.07143], [0.07143, 0.57143]].

**(d)** c = 0.1: P = [[40.007, 0.071], [0.071, 5.714]]; c = 0.01: [[400.0, 0.071], [0.071, 57.14]]. P grows like 1/c and does not exist at c = 0: the undamped oscillator conserves energy, so no quadratic V can fall at the rate |x|².

**Common mistakes.** Solving PA + AᵀP = +Q (sign); forgetting that the (1,2) and (2,1) equations are the same.

## 2. Derivations

**(a)** For Hurwitz A, Z(t) = e^{Aᵀt}Qe^{At} → 0 and Ż = AᵀZ + ZA. Integrating, −Q = Aᵀ∫Z + ∫Z A, so **P = ∫₀^∞ e^{Aᵀt}Qe^{At}dt solves (1)**; xᵀPx = ∫‖Q^{1/2}e^{At}x‖²dt > 0 for x ≠ 0. Conversely, if P > 0, Q > 0 and Av = λv, then v*(AᵀP + PA)v = 2 Re(λ)v*Pv = −v*Qv < 0 forces Re λ < 0.

**(b)** Let β = arcsin(δ/A). By quarter-wave symmetry b₁ = (4/π)∫₀^{π/2} y(θ) sin θ dθ with y = kA sin θ for θ < β and kδ after. ∫₀^β kA sin²θ dθ = (kA/2)(β − sin β cos β); ∫_β^{π/2} kδ sin θ dθ = kδ cos β = kA sin β cos β. Sum: (kA/2)(β + sin β cos β). So **N = b₁/A = (2k/π)[β + sin β cos β] = (2k/π)[arcsin(δ/A) + (δ/A)√(1 − (δ/A)²)]**.

**(c)** With β = arcsin(h/A), y = +M on [β, π + β) and −M elsewhere. b₁ = (1/π)∫ y sin θ dθ = (4M/π)cos β; a₁ = (1/π)∫ y cos θ dθ = −(4M/π)sin β. **N = (4M/(πA))(cos β − j sin β) = (4M/(πA))[√(1 − (h/A)²) − jh/A]**, and −1/N = −(πA/(4M))(cos β + j sin β) = −(π/(4M))(√(A² − h²) + jh): **Im(−1/N) = −πh/(4M)** for every A.

**(d)** θ̈ = (mgl sin θ − bθ̇ + u)/(ml²) = v ⇔ **u = ml²v − mgl sin θ + bθ̇**. It needs m, l, b and g exactly, θ and θ̇ measured exactly and without delay, and as much torque as the formula asks.

## 3. Thruster limit cycle

**Assumptions.** K = 2 + 15/10 = **3.5 rad/s²**, τ = 0.05 s, relay ±1.

**(a)** arg G(jω) = −90° − arctan(0.05ω) > −180° for every ω (the minimum on a grid to 10⁴ rad/s is −179.89°): the curve never meets the negative real axis, so no cycle is predicted. The simulation does not settle either: it **chatters**, the switching intervals shrinking until `simulate_relay_loop` stops at more than 2000 switches by t = 1.04 s. An ideal relay on this plant produces infinitely fast switching (a sliding motion), which the describing function cannot describe.

**(b)** h = 0.01 rad: predicted **A = 0.02875 rad (1.65°) at ω = 53.91 rad/s**.

**(c)** Simulated (after 4 s): **A = 0.02785 rad at 54.23 rad/s**: prediction 3.2 % high in amplitude, 0.6 % low in frequency. Each thruster fires once per cycle: 54.23/π = **17.3 firings per second** in total.

**(d)** A smaller hysteresis shrinks the cycle but speeds it up: h = 0.005 rad gives 0.0180 rad at 69.0 rad/s (22.0 firings/s); h = 0.02 rad gives 0.0462 rad at 41.7 rad/s (13.3 firings/s). Going from h = 0.01 to 0.005 rad cuts the predicted amplitude by 38 % (0.0288 to 0.0180 rad) and raises the predicted firing rate by 28 % (17.2 to 22.0 per second), which wears valves and wastes gas on minimum impulse bits; real roll-control systems add a dead zone (no firing near zero error) to trade amplitude for gas.

## 4. Your loops (s = 15)

**(a)** m = 1.5 + 15/20 = 2.25; gain 2.25/0.24 = 9.375 V/rad, limit 1 V. Predicted **A = 0.2990 rad at 50.00 rad/s**; simulated **0.3075 rad at 49.07 rad/s**: −2.8 % and +1.9 %.

**(b)**

| Q | c* | at-rest reach | random starts inside that converge |
|---|---|---|---|
| diag(1, 0.15) | 0.0812 | 23.60° | 100 of 100 |
| diag(1, 1) | 0.2525 | 22.13° | 100 of 100 |

(The lesson's Q = diag(1, 0.1) gives 0.0689 and 23.57°.) The level c* is not comparable across Q (V changes with Q); the reach is, and it barely moves.

## 5. Judgment

**(a)** Do not trust it (i) for a plant with a lightly damped resonance near 3ω or 5ω, where the harmonics are amplified instead of filtered, and (ii) for a plant like the double integrator or K/(s(τs + 1)) with an ideal relay, where no intersection exists and the real behaviour is chattering (Problem 3). Also for asymmetric nonlinearities (a bias shifts the cycle) and for several nonlinearities in one loop.

**(b)** Not wrong: it is a proof for the region it states. It is conservative because a quadratic V has ellipse-shaped level sets, while the true region is bounded by the stable manifolds of the unstable equilibria near ±37.7°, and the scan stops at the first point where V̇ ≥ 0. A non-quadratic V (Extension 1), or several Lyapunov functions combined, enlarges it.

**(c)** The command saturates at 6 N·m, the gravity cancellation is incomplete, and the closed loop is nonlinear again; from beyond 37.7° at rest no controller can recover (Problem 6). Feedback linearisation is only as good as the actuator authority and the model; in practice it is used inside a margin of both.

**(d)** Near upright with true mass m_t = 1.1m: m_tl²θ̈ ≈ m_tglθ − bθ̇ + ml²v − mglθ + bθ̇, so θ̈ = (m/m_t)v + (1 − m/m_t)(g/l)θ. With v = −16θ − 5.6θ̇: θ̈ + 5.09θ̇ + (14.55 − 1.78)θ = 0, all coefficients positive (Routh, second order): **stable**, slightly slower. The uncancelled gravity reduces the stiffness but not enough to destabilise for a 10 % error.

## 6. Mismatch and torque limit

**(a)** From 60° with the nominal controller: true mass 1.1m gives a largest deviation of **4.50°** from the linear closed form; 1.2m gives **9.26°**. Both still converge (final state below 3 × 10⁻⁴).

**(b)** With the 6 N·m clip, bisection from rest: feedback-linearising controller **37.70°**, linear controller **37.70°** (bisection resolution 0.02°).

**(c)** At rest the gravity torque is mgl sin θ₀ = 9.81 sin θ₀ N·m. For θ₀ > arcsin(6/9.81) = 37.71° it exceeds the motor's 6 N·m, θ̈ > 0 whatever u is, and the stand falls. Below that limit both controllers saturate initially and then recover. The recovery limit is set by physics, not by the control law; feedback linearisation cannot create torque.

## Check-yourself answers (lesson.md)

1. AᵀP + PA + Q = 0. With Q > 0 and P > 0, V = xᵀPx is a Lyapunov function and A is Hurwitz (asymptotically stable).
2. V̇ = −xᵀQx + 2xᵀP·(higher-order terms); the higher-order terms are dominated only near 0, so V̇ < 0 is guaranteed only in a neighbourhood.
3. V̇ < 0 at every point of the set except 0, so trajectories cannot leave a sublevel set (V cannot increase), and they converge to 0.
4. N(A) = (b₁ + ja₁)/A, the first-harmonic gain for input A sin ωt; a limit cycle needs G(jω)N(A) = −1.
5. N = 4M/(πA); A = 4M|G|/π = 4 × 2 × 0.24/π = 0.611 rad.
6. The saturation's effective gain N(A) falls as the amplitude grows; the oscillation grows only until N(A)|G(jω₁₈₀)| = 1, where the loop is exactly critically stable for that amplitude.
7. u = ml²v − mgl sin θ + bθ̇; it needs an exact model (m, l, b) and enough torque (and exact state measurements).
8. A plant with a resonance near a harmonic of the cycle, an ideal relay on a plant whose phase never reaches −180° (chattering), or an asymmetric nonlinearity.
