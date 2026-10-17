# P09.16  Worked solutions

Worked for **s = 15**. Every number below is produced by `solution.py` (the functions named in each problem) or by `scipy.linalg.solve_continuous_are`; the "independent check" line says how it was verified by a different route. Units: SI; angles in radians inside the code, degrees where stated.

## 1. Scalar regulator

**Assumptions.** a = 15/10 = 1.5 s⁻¹ (unstable), b = 1, q = 1, r = 0.1.

**Derivation.** 2aP − P²/r + q = 0 ⇒ P = r(a + √(a² + q/r)) = 0.1(1.5 + √(2.25 + 10)) = 0.1(1.5 + 3.5) = **0.5**. K = bP/r = **5.0**. Pole a − bK = 1.5 − 5 = **−3.5 s⁻¹** = −√(a² + b²q/r) = −√12.25.

**(b) Limits.** Expensive control, r → ∞: the pole → −|a| = −1.5 s⁻¹ and K → 2a/b = 3 (r = 100: K = 3.0033, pole −1.5033; r = 10⁴: pole −1.500033). The pole does not go to 0 because the cost of an unstable state is infinite; the cheapest stabiliser moves the pole to its mirror image, which costs the least control energy. Cheap control, r → 0: pole ≈ −b√(q/r) → −∞ and K ≈ √(q/r) → ∞: the regulator spends freely and the response becomes arbitrarily fast (until actuator limits, absent from the model, bite).

**(c) Independent check.** `care_hamiltonian` → [[0.5]]; `care_newton_kleinman` → ([[0.5]], 5 iterations); `solve_continuous_are` → [[0.5]].

**Common mistakes.** Taking the negative root (P = 0.1(1.5 − 3.5) = −0.2, a "solution" whose closed-loop pole is +3.5 s⁻¹); writing the pole as −√(q/r) and forgetting a²; claiming expensive control turns the loop off (it cannot, for an unstable plant).

## 2. Hamiltonian and margins

**(a)** With ℋ = xᵀQx + uᵀRu + λᵀ(Ax + Bu): ∂ℋ/∂u = 2Ru + Bᵀλ = 0 ⇒ **u = −½R⁻¹Bᵀλ**. The costate equation is λ̇ = −∂ℋ/∂x = −2Qx − Aᵀλ. With μ = λ/2: ẋ = Ax − BR⁻¹Bᵀμ, μ̇ = −Qx − Aᵀμ, i.e. **d/dt(x, μ) = H(x, μ)** with H = [[A, −BR⁻¹Bᵀ], [−Q, −Aᵀ]].

**(b)** If (x, μ) = (X₁c, X₂c) with X₁ invertible, μ = X₂X₁⁻¹x = Px. Differentiate: μ̇ = Pẋ, so −Qx − AᵀPx = P(A − BR⁻¹BᵀP)x for all x, i.e. **AᵀP + PA − PBR⁻¹BᵀP + Q = 0**. The control is u = −R⁻¹BᵀPx = −Kx. On the stable subspace x(t) decays, so the boundary condition μ(∞) = 0 of the infinite-horizon problem holds.

**(c)** J₂ = [[0, I], [−I, 0]], J₂⁻¹ = −J₂. J₂H = [[−Q, −Aᵀ], [−A, BR⁻¹Bᵀ]], and (J₂H)J₂⁻¹ = [[−Aᵀ, Q], [BR⁻¹Bᵀ, A]]. Meanwhile −Hᵀ = [[−Aᵀ, Q], [BR⁻¹Bᵀ, A]] (Q and BR⁻¹Bᵀ symmetric). So **J₂HJ₂⁻¹ = −Hᵀ**: H is similar to −Hᵀ, which has the eigenvalues of −H. Eigenvalues come in pairs (λ, −λ).

**(d)** From the CARE, −AᵀP − PA = Q − KᵀRK (using PB = KᵀR). Add and subtract jωP: P(jωI − A) + (−jωI − Aᵀ)P = Q − KᵀRK. Multiply on the left by Bᵀ(−jωI − Aᵀ)⁻¹ = G(−jω)ᵀ and on the right by G(jω) = (jωI − A)⁻¹B, and use BᵀP = RK:
RK G(jω) + G(−jω)ᵀKᵀR = G(−jω)ᵀQG(jω) − G(−jω)ᵀKᵀRKG(jω).
With L = KG this rearranges to **(1 + L(−jω))ᵀR(1 + L(jω)) = R + G(−jω)ᵀQG(jω)**. For one input, R > 0 is a scalar and G(−jω)ᵀQG(jω) = G(jω)ᴴQG(jω) ≥ 0, so **|1 + L(jω)| ≥ 1**. The disk |1 + L| < 1 contains (−2, 0) on the real axis, so real-axis crossings are at −x, x ≥ 2, and scaling by k changes the encirclements only when kx = 1, i.e. k ≤ ½: **GM ∈ [½, ∞)**. At |L| = 1, L lies on the unit circle about 0 outside the unit circle about −1; those circles intersect at e^{±j120°}, so ∠L ≥ −120° and **PM ≥ 60°**.

**Independent check.** `test_kalman_return_difference_and_guaranteed_margins` verifies |1 + L| ≥ 1 on 3000 frequencies, PM ≥ 60° and stability for k from 0.51 to 1000 on four random systems.

**Common mistakes.** Forgetting the factor ½ between λ and Px; deriving the margins for multi-input loops without saying they hold loop-at-a-time; concluding PM ≥ 60° from |1 + L| ≥ 1 at frequencies other than the crossover.

## 3. Roll attitude

**Assumptions.** Decoupled roll, L_p = −2 s⁻¹, L_δ = 15 + 15 = **30 s⁻²**; perfect measurements; no aileron servo dynamics.

**(a)** φ_max = 10° = 0.17453 rad, p_max = 30°/s = 0.52360 rad/s, δ_max = 15° = 0.26180 rad: **Q = diag(32.828, 3.6476), R = 14.590**.

**(b)** `lqr` → **K = [1.5000, 0.52869]** (rad of aileron per rad, per rad/s). Closed-loop poles **−3.0354 and −14.825 s⁻¹**: both real (ζ = 1 for each), natural frequencies 3.035 and 14.83 rad/s.

**(c)** Peak aileron **15.00°** at t = 0, exactly K₁ × 10° = 1.5 × 10° = 15°: Bryson's rule with φ_max equal to the initial error puts the initial command at the limit, as for the satellite.

**(d)** **PM = 87.1° at ω_c = 15.98 rad/s**; min |1 + L| = 1.0001 ≥ 1. 2 % settling of φ: 1.365 s.

**(e)** R/10: K = [4.7434, 1.6128], poles −3.0034 and −47.381 s⁻¹, peak aileron 47.43°, PM 88.9° at 48.43 rad/s, settling 1.325 s.

| R | K₁ | slow pole (s⁻¹) | peak δ (°) | PM (°) | ω_c (rad/s) | settling (s) |
|---|---|---|---|---|---|---|
| Bryson | 1.500 | −3.035 | 15.00 | 87.1 | 15.98 | 1.365 |
| Bryson/10 | 4.743 | −3.003 | 47.43 | 88.9 | 48.43 | 1.325 |

Ten times cheaper control bought 0.04 s (3 %) of settling and cost 3.2 times the aileron, far beyond the 15° limit. The slow pole is pinned near −p_max/φ_max = −3 s⁻¹ (the lesson's −√(q₁/q₂) argument), so more authority cannot speed up the bank capture.

**Independent check.** For a 2 × 2 roll system the hand formulas (with L_p ≠ 0) are messier than (9), but the stable eigenvalues of H equal the poles, and `solve_continuous_are` gives the same P (test 2).

**Common mistakes.** Bryson weights in degrees (they then differ by (180/π)² and the design changes); reading the 47° deflection as acceptable because "LQR is optimal".

## 4. Your satellite (s = 15: J = 20 kg·m²)

**(a)** **K = [0.22918 N·m/rad, 3.79734 N·m·s/rad]**; poles **−0.09493 ± 0.04946j rad/s**; hand formulas k₁ = √(q₁/r) = 0.22918 (independent of J), k₂ = √(q₂/r + 2·20·k₁) = 3.79734: identical.

**(b)** 5° error: peak torque **0.0200 N·m**, 2 % settling **42.6 s**, peak wheel-momentum change **0.0743 N·m·s**.

**(c)** 15° error:

| | time with |u| at or above 0.02 N·m (s) | settling (s) | peak rate (°/s) | peak momentum change (N·m·s) |
|---|---|---|---|---|
| torque-limited | 8.45 (at the limit) | 46.15 | 0.608 | 0.2121 |
| unlimited | 4.80 (above the limit, peak 0.06 N·m) | 42.6 | 0.639 | 0.2229 |

Saturation lengthens settling by 3.6 s. More important: in both cases the momentum change exceeds the 0.2 N·m·s capacity, so a 15° slew of this heavier spacecraft would saturate the wheel. The design must slow the slew (a smaller ω_max) or the wheel must be larger.

**(d)** `care_newton_kleinman` takes **15 iterations** from the Bass gain; max relative difference from the Hamiltonian P **4.6 × 10⁻¹⁶**.

**Independent check.** Closed form (9) equals both solvers; the 5° peak torque equals k₁θ₀ = 0.22918 × 0.087266 = 0.0200 N·m.

**Common mistakes.** Expecting k₁ to scale with J; checking the torque limit and not the momentum limit.

## 5. Judgment

**(a)** Three violations on a real spacecraft: the rate and angle are *estimated* (gyro plus star tracker through a filter), not measured; the measurement and the torque command are delayed (attitude determination and wheel-drive electronics); and the plant has dynamics the model omits (wheel motor and current loop, flexible appendages, propellant slosh).

**(b)** With an observer the loop at the plant input is no longer K(sI − A)⁻¹B, and the return-difference inequality no longer applies: Doyle (1978) showed by a counterexample that an LQG regulator can have arbitrarily small margins. Check the margins of the *actual* loop (controller plus observer) with the methods of `P09.10` and `P09.18`, or design for recovery (loop-transfer recovery).

**(c)** PM = 77.2° = 1.347 rad at ω_c = 0.3216 rad/s: **delay margin 4.19 s**. A 0.25 s latency costs 0.25 × 0.3216 rad = 4.6° of phase at crossover, leaving about 73°: not a concern for this slow loop.

**(d)** At 0.8 rad/s the rigid-body loop gain is |L(j0.8)| = 0.394, only a factor 2.5 below 1. A lightly damped array mode multiplies the loop gain near its frequency by roughly 1/(2ζ), so with ζ = 0.005 (typical of a bare structure) a factor of 100: the loop would cross 1 again near 0.8 rad/s with uncontrolled phase. The rigid design is not trustworthy without a notch filter, a lower crossover, or the mode in the model.

## 6. Tightness and the cost curve

**(a)** From Problem 1, K = (a + √(a² + b²q/r))/b → 2a/b as r → ∞ (for a = 1.5: r = 100 gives 3.0033, r = 10⁴ gives 3.000033, r = 10⁸ gives 3.0000000033). With gain kK the closed-loop pole is a − bkK → a − 2ka = a(1 − 2k): negative iff k > ½. Any lower gain margin bound larger than ½ would be contradicted by this family, so **½ is tight**.

**(b)** Satellite (J = 10 kg·m², Bryson weights), x₀ = (5°, 0):

| k | 0.25 | 0.5 | 0.8 | 1 | 1.25 | 2 | 4 |
|---|---|---|---|---|---|---|---|
| J(kK) | 17.2716 | 14.4817 | 13.7642 | 13.6845 | 13.7642 | 14.4817 | 17.2716 |

Minimum at k = 1 (as it must be: J(K) = x₀ᵀPx₀ = 13.6845), and J(kK) = J(K/k) to every printed digit.

**(c)**

| change | J(0.5K) | J(2K) |
|---|---|---|
| Q = diag(q₁, 0), new K | 10.5094 | 10.5094 |
| Q = diag(q₁, 10q₂), new K | 33.3045 | 33.3045 |
| original Q, x₀ = (0, 0.1°/s) | 0.32638 | 0.30551 |

The weights do not break the symmetry; the initial state does. For a pure attitude error the cost curve is symmetric in ln k; with an initial rate, too little gain costs more than too much. A flat, symmetric minimum is good news for a designer unsure of the wheel's torque constant (which scales the effective gain): a ±25 % error in that constant costs only 0.6 % more (13.7642 against 13.6845).

**Common mistakes.** Recomputing K when only k should change in (b); comparing costs from different initial states.

## Check-yourself answers (lesson.md)

1. J = ∫(xᵀQx + uᵀRu) dt; AᵀP + PA − PBR⁻¹BᵀP + Q = 0; K = R⁻¹BᵀP. For the satellite K = [N·m/rad, N·m·s/rad].
2. Q ⪰ 0 so the state cost is never negative (otherwise the "optimum" could be to let errors grow); R ≻ 0 so every control has a price and R⁻¹ exists. With R = 0 the problem has no minimiser: infinite gain is "free", K and P blow up (the cheap-control limit).
3. H = [[A, −BR⁻¹Bᵀ], [−Q, −Aᵀ]]; J₂HJ₂⁻¹ = −Hᵀ makes it similar to −Hᵀ, so eigenvalues come in pairs ±λ; the stable half are the closed-loop poles eig(A − BK).
4. Each iteration solves the Lyapunov equation (A − BK_k)ᵀP + P(A − BK_k) + Q + K_kᵀRK_k = 0, which has a positive semidefinite solution equal to the cost of K_k only if A − BK_k is Hurwitz; a destabilising start gives a meaningless (or non-stabilising) P.
5. Q_ii = 1/x_i,max², R_jj = 1/u_j,max²; θ_max = 2° = 0.034907 rad and ω_max = 0.2°/s give Q = diag(820.70, 82070), R = 1/0.01² = 10000.
6. k₁ = √(q₁/r), k₂ = √(q₂/r + 2Jk₁); only k₂ depends on J.
7. Gain margin [½, ∞) and phase margin ≥ 60°, from the Kalman inequality |1 + L(jω)| ≥ 1 at the plant input.
8. Any two of: an observer (or filter) estimating the state; a time delay; unmodelled actuator or structural dynamics; the loop broken anywhere other than the plant input; several inputs considered simultaneously.
