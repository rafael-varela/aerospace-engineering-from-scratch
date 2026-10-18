# P09.15  Worked solutions

Worked for **s = 15**. Every number below is produced by `solution.py` (the functions named in each problem) or by the probe scripts that call it; simulations are of representative models and every noise record is SYNTHETIC. Units: SI.

## 1. One axis (s = 15: J = 3.5 kg m²)

**(a)** k₁ = Jω_n² = 3.5 × 0.04 = **0.14 N m/rad**, k₂ = 2ζω_nJ = 2 × 0.8 × 0.2 × 3.5 = **1.12 N m s/rad**, **N = k₁ = 0.14**.

**(b)** l₁ = 0.8 + 1.0 = **1.8 s⁻¹**, l₂ = 0.8 × 1.0 = **0.8 s⁻²**.

**(c)** `ackermann` → (0.14, 1.12), `observer_gain` → (1.8, 0.8), `reference_gain` → 0.14. The loop eigenvalues are **−0.16 ± 0.12j, −0.8, −1.0**: the regulator pair (−ζω_n ± jω_n√(1 − ζ²)) and the two observer poles.

**(d)** Nr = 0.14 × 0.087266 = **12.2 mN m** at t = 0. That is above a 10 mN m wheel: this design would saturate on a 5° step. Slow the poles or shape the reference.

**Common mistakes.** Forgetting J in k₁ and k₂; treating the observer poles as affecting the step torque (with a perfect estimate they do not).

## 2. Derivations

**(a)** In controllable canonical form A = [[−a₁, …, −a_n], [1, 0, …], …] and B = e₁, so BK = e₁K changes only the first row: A − BK has first row −(a₁ + k₁, …, a_n + k_n). Its characteristic polynomial is sⁿ + (a₁ + k₁)sⁿ⁻¹ + ⋯ + (a_n + k_n), so k_i = a_{d,i} − a_i places any desired polynomial.

**(b)** With A_K = A − BK, A_K² = A² − ABK − BKA_K, and in general A_Kᵏ = Aᵏ − Σ_{j=0}^{k−1} A^{k−1−j}BKA_Kʲ. Cayley–Hamilton for A_K: 0 = φ_d(A_K) = φ_d(A) − [B, AB, …, A^{n−1}B]·M, where M is a stack of row blocks whose last block is KA_K⁰ = K (the coefficient of A^{n−1}B comes only from the leading term a₀ = 1 with j = 0). Multiplying by e_nᵀW_c⁻¹ picks out that last block: e_nᵀW_c⁻¹φ_d(A) = K. That is (2), and it needs W_c invertible: controllability.

**(c)** eig(A − LC) = eig((A − LC)ᵀ) = eig(Aᵀ − CᵀLᵀ). This is state feedback for the pair (Aᵀ, Cᵀ) with gain Lᵀ, so **Lᵀ = ackermann(Aᵀ, Cᵀ, poles)**. It exists when (Aᵀ, Cᵀ) is controllable, i.e. (A, C) is observable.

**(d)** With x̂ = x − e: ẋ = Ax + B(Nr − K(x − e)) = (A − BK)x + BKe + BNr. ė = ẋ − x̂̇ = Ax + Bu − Ax̂ − Bu − LC(x − x̂) = (A − LC)e. This is (5), block upper triangular, so the eigenvalues are eig(A − BK) ∪ eig(A − LC). (i) e is driven by neither r nor K. (ii) The observer modes are unreachable from r, so they cancel in the r → y transfer function, which equals the state-feedback one. (iii) e(0) = 0 ⇒ e ≡ 0 ⇒ the state-feedback trajectory.

**(e)** At steady state 0 = (A − BK)x_ss + BNr ⇒ y_ss = −C(A − BK)⁻¹BNr = r ⇒ **N = −1/(C(A − BK)⁻¹B)**. For the satellite, A − BK = [[0, 1], [−k₁/J, −k₂/J]]. Solving (A − BK)z = B = (0, 1/J)ᵀ gives z₂ = 0 and −k₁z₁/J = 1/J, so z₁ = −1/k₁, Cz = −1/k₁ and **N = k₁**.

## 3. Pitch (s = 15: a = 4.5)

**(a)** Poles −4.5 ± 4.5j, −6, −20: **K = (−10.5097, 1.1237, 16.875, 0.565)** (on α, q, θ, δ), **N = 16.875** (equal to the θ gain). eig(A − BK) returns the four targets.

**(b)** Observer poles −9 ± 9j, −12, −40: **L = (−481.821, 397.69, 46.3, 143.262)**. The loop eigenvalues are **−40, −20, −12, −9 ± 9j, −6, −4.5 ± 4.5j**: the union.

**(c)** 5° step with x̂(0) = 0: overshoot **114.9 %**, peak elevator command **84.37°**, peak angle of attack **8.46°**. θ(4 s) = 5.0000°.

**(d)** 84° is more than five times the ±15° limit of `P09.11`; this linear result would not happen on the aircraft. The peak is u(0) = Nr = 16.875 × 0.0873 = 1.473 rad: the reference enters as a step multiplied by the large θ gain. The large overshoot (with poles of ζ = 0.707) comes from the zeros of the r → θ path, as in `P09.11` §2.5. Fixes: shape the reference (a ramp or a first-order prefilter) and/or choose slower poles; and keep anti-windup in the implementation.

## 4. Thrust-vector-controlled rocket

**(a)** Open-loop poles **±2 s⁻¹** (eigenvalues of [[0, 1], [4, 0]]): one unstable.

**(b)** The closed-loop polynomial is s² + bk₂s + (bk₁ − a), so **k₁ = (ω_n² + a)/b, k₂ = 2ζω_n/b**. For ζ = 0.7, ω_n = 4: **K = (2.0, 0.56)** (= (16 + 4)/10 and 5.6/10). For the double pole at −2 (ω_n = 2, ζ = 1): **K = (0.8, 0.4)**. Because ω_n² > 0, bk₁ > a: the gimbal's restoring moment per radian must exceed the aerodynamic destabilising moment before any damping can help.

**(c)** C = [1, 0], poles −10, −12: **L = (22, 124)** (l₁ = 22, l₂ = 120 + a = 124).

**(d)** With C = [0, 1], [C; CA] = [[0, 1], [4, 0]] has **rank 2**: observable. `observer_gain` gives **L = (31, 22)** for the same poles. The angular acceleration is ω̇ = aθ + bδ. From the measured rate history and the known gimbal command the observer infers θ = (ω̇ − bδ)/a. The satellite has a = 0, so its rate carries no angle information.

## 5. Judgment

**(a)** True loop, model J = 3 kg m²: with true **J = 3.6**, eigenvalues **−0.5518, −0.3857, −0.0519 ± 0.0660j**; with true **J = 2.4**, **−0.6443, −0.1324, −0.1323 ± 0.0941j**. Both are stable, and **θ_ss/r = 1.0000** in both. In steady state the satellite is at rest, so the wheel torque is zero. A plant at rest with zero torque behaves the same whatever J is, so the observer's prediction matches it, e = 0, and the proportional term balances at θ = r. The inertia error changes the transient, not the final value.

**(b)** Separation follows from the cancellation in ė = (A − LC)e, which uses the plant's own A and B inside the observer. With a model error the (x, e) matrix is no longer block triangular, and the poles mix (part (a): none of the four is at its placed value).

**(c)** Twenty times the regulator's ω_n = 0.1 rad/s is ω₀ = 2 s⁻¹. By the ω₀^{3/2} law, the rate noise would be 3.394 × (2/1.2)^{1.5} = **7.3 arcsec/s** (predicted), against 0.424 for the lesson's ω₀ = 0.3, a factor of 17. The SYNTHETIC torque noise rose from 0.586 to 7.78 μN m between ω₀ = 0.15 and 1.2 and would rise further. The regulator, with a 0.1 rad/s bandwidth, cannot use an estimate that converges in 0.5 s rather than in 5 s: the extra speed buys only noise and sensitivity to model error.

**(d)** The innovation y − Cx̂. If the model is right it is zero-mean and white, with the tracker's noise level. A bias, a correlation or a spectral line (a flexible mode) says the model is wrong.

## 6. Complex observer poles

**(a)** **L = (0.6, 0.18)** (2σ' = 0.6 and σ'² + ω² = 0.09 + 0.09).

**(b)** Fitted rate of ‖e‖: **0.29427 s⁻¹ over 20–60 s** and **0.29822 s⁻¹ over 5–60 s**. ‖e(t)‖ is e^{−0.3t} times a periodic factor of period π/0.3 = 10.5 s. A least-squares line through its logarithm picks up part of that factor's slope at the window ends, and the bias depends on how many half-periods the window holds.

**(c)** The modal-norm fit gives **0.300000 s⁻¹**. The 5 local maxima of ‖e‖ between 5 and 60 s give **0.29991 s⁻¹**. The modal coordinates of a complex pair have |z_k(t)| = |z_k(0)|e^{σt} exactly, so the modal fit is unbiased. The peaks fall on nearly the same phase of the periodic factor in every cycle, so the peak fit is nearly unbiased.

**(d)** e(t) = c v e^{λt} + c̄ v̄ e^{λ̄t} = 2Re(c v e^{(σ+jω)t}) = e^{σt}(p cos ωt + q sin ωt) with real constant vectors p and q. So ‖e‖ = e^{σt}√(|p|²cos²ωt + 2pᵀq cos ωt sin ωt + |q|²sin²ωt), and the square root is periodic with period π/ω. log‖e‖ = σt + g(t) with g periodic. Over a whole number of half-periods g has no net trend, and the fit recovers σ up to the end effects of the least-squares weighting. Over any other window, g's partial slope biases the result, by 2 % here.

## Check-yourself answers (lesson.md)

1. K = e_nᵀW_c⁻¹φ_d(A), with W_c = [B, AB, …, A^{n−1}B] and φ_d the desired characteristic polynomial. It requires (A, B) controllable (W_c invertible).
2. k₁ = Jω_n², k₂ = 2ζω_nJ.
3. With one input, the n gain entries are fixed by the n polynomial coefficients. With m inputs there are nm entries and only n conditions, and the freedom is used to shape the eigenvectors (robustness), as `place_poles` does.
4. L = ackermann(Aᵀ, Cᵀ, poles)ᵀ (duality). It requires (A, C) observable.
5. The loop with u = Nr − Kx̂ has the eigenvalues eig(A − BK) ∪ eig(A − LC). Proof: change to (x, e = x − x̂), which makes the matrix block upper triangular.
6. ė = (A − LC)e: subtracting the observer from the plant cancels the input Bu exactly, so neither r nor K appears.
7. At the slowest observer pole, 0.3 s⁻¹, once the −0.6 mode has decayed. Measure it by fitting log‖x − x̂‖ against t over a window that starts after the fast mode is gone and ends before rounding noise: the lesson's fit over 20–60 s gives 0.29998 s⁻¹.
8. Rate-estimate noise grows like ω₀^{3/2} (each doubling about ×2.8), it reaches the wheel as torque jitter, and the model error grows in importance. The regulator cannot use estimates much faster than its own bandwidth.
