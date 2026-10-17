# P09.13  Worked solutions

Worked for **s = 15**. Every number below is produced by `solution.py` (the functions named in each problem) or by the probe scripts that call it; simulations are of representative models. Units: SI.

## 1. Second-order sensor (s = 15: ω_n = 3.5 rad/s, ζ = 0.3)

**(a)** Controllable: **A = [[−2.1, −12.25], [1, 0]], B = (1, 0)ᵀ, C = [0, 12.25], D = 0**. Observable: **A = [[−2.1, 1], [−12.25, 0]], B = (0, 12.25)ᵀ, C = [1, 0], D = 0**.

**(b)** G(jω_n) = ω_n²/(ω_n² − ω_n² + 2jζω_n²) = 1/(2jζ) = **−1.6667j**: magnitude **1/(2ζ) = 1.667**, phase −90°. `ss_eval` gives −1.6667j for both realisations (real part below 10⁻¹⁶). At ω_n the spring and mass terms cancel and only the damping limits the response: the resonance amplification of `P02.03`.

**(c)** The physical state z = (y, ẏ) relates to the companion state x = (ẇ, w) by y = ω_n²w and ẏ = ω_n²ẇ, so z = T⁻¹x with **T⁻¹ = [[0, 12.25], [12.25, 0]]**, T = [[0, 1/12.25], [1/12.25, 0]]. `similarity(A, B, C, D, T)` returns A_p = [[0, 1], [−12.25, −2.1]], B_p = (0, 12.25)ᵀ, C_p = [1, 0]: exactly the physical model ẏ = ẏ, ÿ = −ω_n²y − 2ζω_nẏ + ω_n²u.

**Common mistakes.** Putting y in the first companion state (it is the highest derivative of w); forgetting that T and T⁻¹ swap roles between x = Tz and z = T⁻¹x.

## 2. Derivations

**(a)** C T(sI − T⁻¹AT)⁻¹T⁻¹B + D = C T[T⁻¹(sI − A)T]⁻¹T⁻¹B + D = C T T⁻¹(sI − A)⁻¹T T⁻¹B + D = C(sI − A)⁻¹B + D. The eigenvalues are the roots of det(sI − T⁻¹AT) = det T⁻¹ · det(sI − A) · det T = det(sI − A).

**(b)** With w⁽ⁿ⁾ = u − a₁w⁽ⁿ⁻¹⁾ − ⋯ − a_nw and y = b₀w⁽ⁿ⁾ + b₁w⁽ⁿ⁻¹⁾ + ⋯ + b_nw, substitution gives y = b₀u + (b₁ − b₀a₁)w⁽ⁿ⁻¹⁾ + ⋯ + (b_n − b₀a_n)w. With x = (w⁽ⁿ⁻¹⁾, …, w) that is (3). **D = b₀** is the part of the output that responds instantly to u, which exists only when deg num = deg den.

**(c)** For a scalar G, G = Gᵀ = [C(sI − A)⁻¹B]ᵀ + D = Bᵀ(sI − Aᵀ)⁻¹Cᵀ + D. That is the TF of (Aᵀ, Cᵀ, Bᵀ, D). For MIMO, Gᵀ is a different transfer matrix (outputs and inputs swapped), so the dual is a different system.

**(d)** (sI − A)(N₁sⁿ⁻¹ + ⋯ + N_n) = (sⁿ + c₁sⁿ⁻¹ + ⋯ + c_n)I. The coefficient of sⁿ gives N₁ = I; of s^{n−k}, N_{k+1} − AN_k = c_kI, i.e. **N_{k+1} = AN_k + c_kI**; of s⁰, −AN_n = c_nI. Unrolling the recursion, N_n = A^{n−1} + c₁A^{n−2} + ⋯ + c_{n−1}I, so −AN_n − c_nI = −(Aⁿ + c₁A^{n−1} + ⋯ + c_nI) = 0: **Cayley–Hamilton**. The trace step (c_k = −tr(AN_k)/k) comes from tr adj(sI − A) = (d/ds) det(sI − A).

**(e)** From §2.1: **θ/δ = M_δ(s − Z_α)/(sΔ) = 12(s + 1.2)/(s(s² + 3.7s + 9))**, type 1, relative degree 2 (rolls off at −40 dB/decade). **h/δ = 432/(s²(s² + 3.7s + 9))**, type 2, relative degree 4 (−80 dB/decade).

## 3. Hidden mode

**(a)** **A = [[−3, −2], [1, 0]], C = [1, 1]**, eigenvalues **−2 and −1**.

**(b)** Largest difference between the 2-state and 1-state step responses: **6.1 × 10⁻¹⁵**; y(5 s) = **0.49998** (→ 1/2).

**(c)** From x₀ = (1, 1): y(0) = **2.0**, y(1) = **0.27067**, y(5) = **9.08 × 10⁻⁵**. These are exactly 2e^{−2t} (2e^{−2} = 0.27067, 2e^{−10} = 9.08 × 10⁻⁵): the output shows only the mode at −2. The mode at −1 is in the state but invisible at the output.

**(d)** rank [C; CA] = **1**, rank [B, AB] = **2**. The realisation is controllable but not observable: the pole-zero cancellation has hidden the mode at −1 from the output (`P09.14`).

## 4. Flexible array (s = 15: J₂ = 2.5 kg m²)

**(a)** J₁θ̈₁ = u − k(θ₁ − θ₂) − c(θ̇₁ − θ̇₂) and J₂θ̈₂ = k(θ₁ − θ₂) + c(θ̇₁ − θ̇₂). With x = (θ₁, θ₂, ω₁, ω₂): A rows [0, 0, 1, 0], [0, 0, 0, 1], [−k/J₁, k/J₁, −c/J₁, c/J₁], [k/J₂, −k/J₂, c/J₂, −c/J₂]; B = (0, 0, 1/3, 0)ᵀ, C = [1, 0, 0, 0], D = 0.

**(b)** `ss_to_tf`: **θ₁/u = (0.33333 s² + 0.0026667 s + 0.26667)/(s⁴ + 0.014667 s³ + 1.46667 s²)** = (1/3)(s² + 0.008s + 0.8)/(s²(s² + 0.01467s + 1.4667)). Poles: **0, 0** (the rigid body; computed as ±8 × 10⁻⁹j) and **−0.00733 ± 1.21104j** (the flexible mode). Zeros: **−0.004 ± 0.89442j** (the antiresonance). Checks: √(k(J₁ + J₂)/(J₁J₂)) = √(2 × 5.5/7.5) = **1.2111** rad/s and √(k/J₂) = √0.8 = **0.8944** rad/s. At the antiresonance the array swings against the spring so that the hub, where the sensor and actuator sit, stands still.

**(c)** θ₁(30 s) = **0.81905 rad** for a 0.01 N m step. SciPy differs by **1.3 × 10⁻¹⁴ rad**. The rigid-body result 0.01 × 30²/(2 × 5.5) = **0.81818 rad** differs by 0.9 mrad: the flexible mode, barely damped, oscillates on top.

**(d)** No. cond(V) = **1.9 × 10⁸**: the double eigenvalue at 0 has a single eigenvector (angle and rate of the rigid body form a chain of integrators). A is defective and has no eigenvector basis. The modal form exists only for the flexible part, or with a Jordan block for the rigid body.

**Common mistakes.** Forgetting the reaction on the array (the spring acts on both bodies with opposite signs); taking the zero frequency as the resonance.

## 5. Judgment

**(a)** The lesson's study: rounding to five figures leaves a 6th-order companion form accurate to 5.5 × 10⁻¹³, but at 8th order the poles move by 0.64 and at 10th order by 2.75, with spurious complex pairs (imaginary parts up to 2.42). A 10th-order flexible model in companion form at five figures has unreliable poles and could show a false instability. Store it in modal form (real 2 × 2 blocks [[σ, ω], [−ω, σ]] with mode shapes), or in physical coordinates, where each number moves one pole by about its own error.

**(b)** A TF is enough for: the frequency response and the margins of a loop closed around h; the step response from trim. The state model is needed for: the angle of attack during the manoeuvre (stall margin), a response from a disturbed initial state (a gust upset), and any state-feedback or observer design (`P09.15`), and also for elevator plus throttle in one model.

**(c)** The model predicts an altitude growing faster than t² (a loop at 1.6°/s) for as long as the elevator is held. In 30 s the pitch angle would pass 45°, outside the small-angle model, and the real aircraft would have slowed down long before. The missing state is the **airspeed perturbation u**; the derivatives **X_u, Z_u** (and M_u) couple it in and create the phugoid mode.

**(d)** **B_g = (Z_α, M_α, 0, 0)ᵀ = (−1.2, −6, 0, 0)ᵀ**. `ss_to_tf` gives h/w_g = (36s² + 90s)/(s²(s² + 3.7s + 9)) = **36(s + 2.5)/(s(s² + 3.7s + 9))**. By hand: α + w_g = w_g s(s − M_q)/Δ and γ = −Z_α(α + w_g)/s, so h/w_g = −VZ_α(s − M_q)/(sΔ). The poles are the same. The numerator gains a zero at s = M_q = −2.5, which cancels one integrator: a gust changes the flight path directly, not only through pitch.

## 6. Faddeev–LeVerrier accuracy

**(a)** A = Q diag(1, …, n)Qᵀ, Q orthogonal from seed 1315. Largest relative coefficient error:

| n | 5 | 10 | 15 | 20 |
|---|---|---|---|---|
| Faddeev–LeVerrier | 2.5 × 10⁻¹⁵ | 7.0 × 10⁻¹⁴ | 6.1 × 10⁻¹² | 2.1 × 10⁻⁸ |
| `np.poly` | 2.7 × 10⁻¹⁵ | 3.7 × 10⁻¹⁵ | 3.3 × 10⁻¹⁵ | 2.3 × 10⁻¹⁵ |

`np.poly` computes the eigenvalues first (backward stable, error 10⁻¹³ here) and multiplies out. Faddeev–LeVerrier accumulates products Aᵏ whose norms grow like nᵏ. The traces then cancel catastrophically, and the error grows by one and a half to three and a half decades every five orders.

**(b)** With Q = I every N_k is diagonal with integer entries and every trace is an integer smaller than 2⁵³. The division by k is exact because c_k is an integer, so the arithmetic is exact: the error is **0** for n = 5, 10, 15. At n = 20 the coefficients exceed 2⁵³ ≈ 9 × 10¹⁵ (20! = 2.4 × 10¹⁸) and rounding appears: **1.2 × 10⁻⁹** (`np.poly` 1.3 × 10⁻¹⁶).

**(c)** Roots of the (exact up to n = 15) Faddeev–LeVerrier coefficients of diag(1, …, n), by `np.roots`:

| n | 5 | 10 | 15 | 20 |
|---|---|---|---|---|
| largest root error | 2.2 × 10⁻¹³ | 2.8 × 10⁻⁹ | 9.2 × 10⁻⁶ | 1.3 × 10⁻¹ |

The coefficients are exact, but the roots of a polynomial with clustered roots are violently sensitive to coefficient perturbations, including the rounding inside the root finder. This is Wilkinson's polynomial: the condition number of a root with respect to the coefficients grows factorially with n.

**(d)** For a high-order model, never go through the characteristic polynomial if you can avoid it. Keep the model in state space (eigenvalues from A directly) or in factored pole-zero form. Use `ss_to_tf` only for low orders (up to about 6) or for display.

## Check-yourself answers (lesson.md)

1. A = [[−3, −2], [1, 0]], B = (1, 0)ᵀ, C = [0, 1], D = 0 (state (ẏ, y)).
2. Any invertible T gives x = Tz with (T⁻¹AT, T⁻¹B, CT, D), and the TF is unchanged. There are infinitely many T.
3. CT(sI − T⁻¹AT)⁻¹T⁻¹B = CT T⁻¹(sI − A)⁻¹T T⁻¹B = C(sI − A)⁻¹B.
4. The real part σ (decay rate, so ζω_n = −σ) and the damped frequency ω of the mode.
5. Two at s = 0: θ = ∫q (attitude kinematics) and h = ∫V(θ − α) (altitude from climb rate).
6. With airspeed held constant, the steady short period under a fixed elevator has α = 4/3 δ and q = 1.6 δ per second, which never decays. The missing physics is the airspeed state (phugoid): a real aircraft would slow down.
7. Every square matrix satisfies its own characteristic polynomial: Aⁿ + c₁Aⁿ⁻¹ + ⋯ + c_nI = 0. In Faddeev–LeVerrier it is the last step, AN_n + c_nI = 0, a residual that should be zero to rounding.
8. Its poles depend on coefficients spanning many decades, so tiny coefficient errors move them by large amounts. In the lesson, rounding to five figures moved the poles of a 10th-order model by up to 2.75 and created spurious oscillatory pairs.
