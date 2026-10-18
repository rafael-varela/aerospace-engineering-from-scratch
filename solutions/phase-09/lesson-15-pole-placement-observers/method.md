# P09.15  Method

The route to each answer, without the final numbers. Pair with `hints.md`; `worked-solutions.md` has the numbers for s = 15.

## 1. One axis

(a) k₁ = Jω_n², k₂ = 2ζω_nJ, N = k₁. (b) Observer poles −p₁, −p₂: l₁ = p₁ + p₂, l₂ = p₁p₂. (c) `ackermann(A, B, second_order_poles(0.8, 0.2))`, `observer_gain(A, C, [-0.8, -1.0])`, `reference_gain`, and `np.linalg.eigvals(observer_closed_loop(...)[0])`. (d) Nr with r = 5° in radians.

## 2. Derivations

(a) In canonical form A has first row −a and B = e₁, so A − BK has first row −(a + K): the closed-loop coefficients are a_i + k_i. To place the poles, choose k_i = a_{d,i} − a_i. (b) Cayley–Hamilton gives φ_d(A − BK) = 0. Expanding the powers of A − BK, φ_d(A − BK) = φ_d(A) − W_c M for some stacked matrix M whose last block row is K. Multiplying by e_nᵀW_c⁻¹ and using e_nᵀW_c⁻¹W_c M = K gives K = e_nᵀW_c⁻¹φ_d(A). (c) eig(M) = eig(Mᵀ), and (A − LC)ᵀ = Aᵀ − CᵀLᵀ is a state-feedback problem for (Aᵀ, Cᵀ) with gain Lᵀ. (d) With x̂ = x − e: ẋ = Ax + B(Nr − Kx + Ke) and ė = ẋ − x̂̇ = (A − LC)e. The matrix is block triangular. The three consequences follow from the zero lower-left block and the zero lower input block. (e) x_ss = −(A − BK)⁻¹BNr, so y_ss = −C(A − BK)⁻¹BNr = r. For the satellite, (A − BK)⁻¹B gives C(A − BK)⁻¹B = −1/k₁.

## 3. Pitch

Use `char_poly` to check that the targets are a valid set. Get K from `ackermann`, N from `reference_gain`, L from `observer_gain(A, C, 2*reg)`, and the eight eigenvalues from `observer_closed_loop`. Simulate with r = 5° for 4000 steps at 1 ms. The overshoot comes from x[:, 2], the peak elevator command from |u| and the peak α from |x[:, 0]|. The peak command is u(0) = Nr, with N equal to the gain on θ. Fixes: slower poles, reference shaping (a ramp or a prefilter), or putting the reference through the observer (feeding r as a model-following command) rather than as an N step; and in any case anti-windup if the limit is reached.

## 4. Rocket

(a) The eigenvalues of [[0, 1], [a, 0]] are ±√a. (b) The polynomial is s² + bk₂s + (bk₁ − a), so k₁ = (ω_n² + a)/b and k₂ = 2ζω_n/b, and for the double pole ω_n = 2, ζ = 1. Since ω_n² > 0, k₁ > a/b always: the feedback must first overcome the aerodynamic instability. (c) `observer_gain(A, [[1, 0]], [-10, -12])`. (d) [C; CA] with C = [0, 1] gives [[0, 1], [a, 0]], of rank 2 when a ≠ 0. Then `observer_gain` with that C. Physically ω̇ = aθ + bδ: the measured angular acceleration minus the known control term reveals θ. A satellite (a = 0) has no such term.

## 5. Judgment

(a) Form the 4 × 4 true-loop matrix for each J and take its eigenvalues. Solve for the steady state with the reference input; θ/r comes out exactly 1 because at steady state ω = 0, the observer error is zero (the plant is at rest, so model and plant agree on a stationary state), and N = k₁ cancels the θ feedback. (b) Separation needs the observer to use the same A and B as the plant. (c) Use ω₀^{3/2}: twenty times faster than ω_n = 0.1 is ω₀ = 2. Extrapolate from 3.394 arcsec/s at ω₀ = 1.2 by (2/1.2)^{1.5}, and point to the torque noise. (d) The innovation sequence y − Cx̂: its mean, whiteness and size against the sensor noise.

## 6. Complex observer poles

(a) `observer_gain(A, C, [-0.3+0.3j, -0.3-0.3j])`. (b) `fit_decay_rate` over both windows. The norm of e oscillates at the damped frequency, so a straight-line fit of log‖e‖ depends on where the window cuts the oscillation. (c) z = V⁻¹e from `np.linalg.eig(A - L @ C)`, and fit log‖z‖. For the peaks use `argrelmax` on ‖e‖ restricted to 5–60 s and fit log‖e‖ at those times. The modal norm is exactly e^{σt}|z(0)|, so its fit is unbiased; the peak fit is nearly so. (d) e(t) = 2Re(c v e^{(σ+jω)t}) = e^{σt}(a cos ωt + b sin ωt), with constant real vectors a and b. Its norm is e^{σt} times √(quadratic form in cos and sin), which is periodic with period π/ω. The least-squares slope over a window that is not a whole number of half-periods picks up part of that oscillation's slope.
