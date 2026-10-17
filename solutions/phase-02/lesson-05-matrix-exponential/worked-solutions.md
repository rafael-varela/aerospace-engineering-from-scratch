# P02.05  Worked solutions

Worked for **s = 15**. Every number below is produced by `solution.py` (the functions named in each problem) or by `scipy.linalg.expm`; the "independent check" line says how it was verified by a different route. Units: SI; angles in radians inside the code.

## 1. Exponentials by hand

**Derivation.** (a) N² = 0 ⇒ e^{Nt} = I + Nt = [[1, t], [0, 1]]. (b) ω = 1.5 rad/s. A² = −ω²I, so Σ A^k t^k/k! = (Σ(−1)^k(ωt)^{2k}/(2k)!) I + (Σ(−1)^k(ωt)^{2k+1}/(2k+1)!) A/ω = cos(ωt) I + sin(ωt) A/ω = [[cos ωt, sin ωt], [−sin ωt, cos ωt]]: orthogonal, determinant cos² + sin² = 1 = e^{tr A·t}. At t = 1: cos 1.5 = 0.0707, sin 1.5 = 0.9975, so **e^{A} = [[0.0707, 0.9975], [−0.9975, 0.0707]]**. (c) e^{Dt} = diag(e^{−t}, e^{−3t}).

**Result.** `expm_taylor(2N, 5)` = [[1, 2], [0, 1]] exactly; `expm_scaling_squaring(A)` matches the rotation to 2 × 10⁻¹⁶; `expm_eig(A)` matches to 2 × 10⁻¹⁶ (eigenvalues ±1.5i, two eigenvectors) and `expm_eig(D)` is trivially diag(e^{−1}, e^{−3}); **`expm_eig(N)` raises ValueError**: N has the double eigenvalue 0 with a single eigenvector (1, 0), so V is singular and e^{Nt} = I + Nt contains the t that no diagonal e^{Λt} can supply.

**Independent check.** The rotation by 1.5 rad applied to (1, 0) gives (cos 1.5, −sin 1.5): a clockwise rotation, consistent with the sign pattern of A (ẋ₁ = ωx₂, ẋ₂ = −ωx₁).

**Common mistakes.** Writing e^{At} = [[e⁰, e^{ωt}], [e^{−ωt}, e⁰]] (exponentiating entry by entry); forgetting the 1/ω in the odd terms; expecting `expm_eig` to work for N because "0 is an eigenvalue twice".

## 2. Properties of Φ

**(a)** Term-by-term: d/dt Σ_{k≥0} A^k t^k/k! = Σ_{k≥1} A^k t^{k−1}/(k−1)! = A Σ_{j≥0} A^j t^j/j! = AΦ; since A^k = A^{k−1}A the factor can equally be taken out on the right: ΦA. Φ(0) = I because only the k = 0 term is non-zero at t = 0.
**(b)** X(t) = Φ(t)Φ(s): Ẋ = AΦ(t)Φ(s) = AX, X(0) = Φ(s). Y(t) = Φ(t + s): Ẏ = AY, Y(0) = Φ(s). The matrix IVP Ẏ = AY, Y(0) given has a unique solution (`P02.01`, applied column by column), so X ≡ Y. With s = −t: Φ(t)Φ(−t) = Φ(0) = I, hence Φ⁻¹(t) = Φ(−t) and Φ(t) is never singular.
**(c)** Jacobi: d(det Φ)/dt = det Φ · tr(Φ⁻¹Φ̇) = det Φ · tr(Φ⁻¹AΦ) = det Φ · tr A (trace is invariant under similarity). Then det Φ = e^{tr A · t} det Φ(0) = e^{tr A · t}. (For diagonalisable A directly: det(V e^{Λt} V⁻¹) = Π_k e^{λ_k t} = e^{(Σλ_k)t} and Σλ_k = tr A.) Check: det Φ_CW(P/2) = 1.000000000000 in `main()` since tr A_CW = 0.
**(d)** e^A e^B = Σ_{j,k} A^jB^k/(j!k!). If AB = BA, (A + B)^m = Σ_j C(m, j) A^j B^{m−j}, so the coefficient of A^jB^k in e^{A+B} is C(j+k, j)/(j+k)! = 1/(j!k!): the series agree. Counterexample: A = [[0, 1], [0, 0]], B = [[0, 0], [1, 0]], AB = diag(1, 0) ≠ BA = diag(0, 1). e^A = I + A, e^B = I + B, **e^A e^B = [[2, 1], [1, 1]]**; (A + B)² = I so e^{A+B} = cosh(1) I + sinh(1)(A + B) = **[[1.5431, 1.1752], [1.1752, 1.5431]]**. Different (and e^{A+B} is symmetric while e^Ae^B is not).
**(e)** M = [[A, B], [0, 0]]: M² = [[A², AB], [0, 0]] and by induction M^k = [[A^k, A^{k−1}B], [0, 0]] for k ≥ 1. e^{Mh} = I + Σ_{k≥1} M^k h^k/k! has upper-left I + Σ A^k h^k/k! = e^{Ah}, upper-right Σ_{k≥1} A^{k−1}h^k/k! · B = (∫₀^h Σ A^{k−1}s^{k−1}/(k−1)! ds) B = (∫₀^h e^{As} ds) B, lower-right I. Check: `discretize(OSC, B, 0.05)` reproduces `expm(OSC·0.05)` to 10⁻¹⁰ and the quadrature of the integral to 10⁻⁷ (test 7).

**Common mistakes.** Assuming Φ(t)Φ(s) = Φ(s)Φ(t) needs proof (it follows from (b): both equal Φ(t + s)); using Jacobi's formula with tr(Φ̇Φ⁻¹) in the wrong order (it does not matter for the trace, but say why); claiming e^{A+B} = e^Ae^B whenever the two "look small".

## 3. Nutation

**Assumptions.** Linearised rigid-body transverse dynamics about a steady spin; λ = 0.5 + 15/50 = **0.8 rad/s**, d = 0.02 s⁻¹.

**Derivation.** (a) A = [[−d, λ], [−λ, −d]] = −dI + λJ with J = [[0, 1], [−1, 0]]. I commutes with J, so by 2(d) e^{At} = e^{−dt} e^{λJt}, and e^{λJt} = R(λt) = [[cos λt, sin λt], [−sin λt, cos λt]] from 1(b). (b) Nutation period 2π/0.8 = **7.854 s**; the amplitude √(ω₁² + ω₂²) = e^{−dt} × initial halves at t = ln 2/0.02 = **34.66 s**, which is 34.66/7.854 = **4.41 cycles**. (c) e^{−0.2} = 0.8187, cos 8 = −0.1455, sin 8 = 0.9894: **e^{10A} = [[−0.1191, 0.8100], [−0.8100, −0.1191]]**, determinant 0.1191² + 0.8100² = **0.6703 = e^{−0.4}**. (d) 0.1 e^{−dt} = 0.001 ⇒ t = ln 100/0.02 = **230.3 s**, independent of λ (the rotation changes the direction of the transverse rate vector, not its length). `is_stable(A)` → **True** (eigenvalues −0.02 ± 0.8i).

**Independent check.** `state_transition(A, 10.0)` → [[−0.11913, 0.81002], [−0.81002, −0.11913]], det 0.670320 (= e^{−0.4} = 0.670320). The decay time from `transient_decay_time`-style reasoning: amplitude ratio 100 at rate 0.02 s⁻¹ ⇒ 230 s.

**Common mistakes.** Treating λ as a damping; forgetting that e^{−dI t} = e^{−dt}I is a scalar factor (which is why the split works); quoting the half-life in cycles without converting.

## 4. Your rendezvous

**Assumptions.** Clohessy–Wiltshire in-plane model, n = 1.131401 × 10⁻³ rad/s (P = 92.56 min), impulsive burns, (x₀, y₀) = (0, −(500 + 300)) = **(0, −800) m**, at rest.

**Result.** (a) After one orbit (and at every time) the chaser is still at **(0, −800) m**: with x₀ = ẋ₀ = ẏ₀ = 0 the only contribution is the y₀ column of Φ, which is (0, 1, 0, 0)ᵀ for all t. Same altitude, same speed, same orbit: a co-orbital chaser never drifts. (b) T = 0.4P = 2221.4 s = **37.02 min**. `rendezvous_delta_v` → **Δv₁ = (−0.3262, +0.0530) m/s**, **Δv₂ = (−0.3262, −0.0530) m/s**, total **0.661 m/s**. (c) `propagate` every T/8 (minutes, x m, y m): 0: (0, −800); 4.63: (−84.5, −758.0); 9.26: (−151.6, −668.1); 13.88: (−194.6, −543.2); 18.51: (−209.4, −400.0); 23.14: (−194.6, −256.8); 27.77: (−151.6, −131.9); 32.40: (−84.5, −42.0); 37.02: (0.0, 0.0). Maximum radial excursion **−209.4 m** at mid-transfer: the chaser drops *below* the station to a faster orbit in order to catch up, and the path is symmetric about the midpoint (the transfer is a half of a CW ellipse). (d) Totals: 0.25 orbit **1.231 m/s**, 0.5 orbit **0.453 m/s**, 0.75 orbit **0.183 m/s**. (e) `rendezvous_delta_v(n, r0, v0, P)` raises `ValueError`; det Φ_rv(T) = [8(1 − cos nT) − 3nT sin nT]/n² is 0 at nT = 2π: Φ_rv(P) = [[0, 0], [0, −6π/n]] has rank 1 (any velocity change returns x to its starting value after a full orbit), so no first burn can bring x to zero at that time.

**Independent check.** The arrival position after the planned burn is 1.9 × 10⁻¹⁰ m (rounding); `solve_ivp` of `cw_rhs` with the post-burn state reaches the origin to 10⁻⁵ m (test 8). The symmetry Δv₂ = (Δv₁,x, −Δv₁,y) follows from the time-reversal symmetry of the equations for a transfer between two states at rest at the same altitude. At T = 0.4P the determinant formula gives 10.04/n² = 7.84 × 10⁶ m², far from singular.

**Common mistakes.** Expecting the chaser in (a) to drift (only a *radial* offset or a velocity drifts); adding Δv components instead of norms for the total; concluding from (d) that "longer is always cheaper" (it is singular at one orbit and 1.25 orbits costs more than 0.75).

## 5. Judgment

**(a)** 50 states, 10 000 times, structural (hence diagonalisable with distinct damped frequencies): `modal_decomposition` once (O(50³)), then Φ(t_k) = V diag(e^{λt_k}) V⁻¹ at O(50²) each, or better, propagate the modal coordinates and never form Φ at all. For the 4-state CW matrix at one T: the matrix is defective, so the eigen method is out; use the closed form (exact) and `expm_scaling_squaring` or `scipy.linalg.expm` as the check.
**(b)** ‖At‖₁ = 4.8 × 60 = 288: the Taylor terms grow until k ≈ 288 and the largest is of order e^{288} ≈ 10¹²⁵; after 30 terms the partial sum (≈ 288³⁰/30! ≈ 10³⁰) is nowhere near converged and bears no relation to the answer. The true entries are of order e^{−60} ≈ **10⁻²⁶** (`expm` gives 8.8 × 10⁻²⁷ as the largest). Use scaling and squaring (10 squarings) or note that the answer is zero to any engineering purpose.
**(c)** For the 1 km, half-orbit (46 min, 0.52 m/s) transfer: burn execution error at 1 % of Δv, 5.2 mm/s × 2777 s ≈ **14 m**; differential drag of 10⁻⁶ m/s² (dissimilar area-to-mass at 400 km), ½ a T² ≈ **3.9 m**; linearisation, of order (separation/a) × separation ≈ (1000/6.78 × 10⁶) × 1000 ≈ **0.15 m** per orbit. Burn execution dominates; the remedy is a mid-course re-plan, not a better matrix.
**(d)** Between updates the controller's output is literally constant, so the zero-order-hold map (3) is the exact solution of the equations the plant actually experiences: no truncation error at any step length. RK4 with h = 10 s evaluates the right-hand side at the step midpoint and end as if the input varied smoothly, and its error is O(h⁴) in the dynamics plus an inconsistency in how it treats the hold; for the CW matrix at h = 10 s (‖Ah‖ ≈ 0.02) that error is small but not zero. RK4 is the right choice when the plant is nonlinear (`P02.06`), when the input varies continuously inside the step (a ramp or a sinusoid), or when A changes with time.

## 6. Exact versus Euler

**Assumptions.** k = 15² = **225 s⁻²**, c = 0.2 × 15 = **3 s⁻¹**, m = 1: ω_n = 15 rad/s, ζ = 3/(2 × 15) = 0.1, ω_d = 14.925 rad/s.

**Derivation.** (a) λ = −1.5 ± 14.925i s⁻¹. A_d = e^{Ah} = V e^{Λh} V⁻¹ (A is diagonalisable: distinct eigenvalues), so its eigenvalues are e^{λh} with modulus e^{−1.5h} < 1 for every h > 0. (b) A_E = I + hA has eigenvalues 1 + hλ; |1 + hλ|² = (1 − 1.5h)² + (14.925h)² = 1 − 3h + 225h² < 1 ⇔ **h < 3/225 = 2ζ/ω_n = 0.01333 s**. In general h_max = −2 Re λ/|λ|².

**Result.** `discretize(A, [[0], [1]], h)` and `I + hA`:

| h (s) | max |eig A_d| | e^{−1.5h} | max |eig A_E| |
|---|---|---|---|
| 0.01 | 0.98511 | 0.98511 | 0.99624 |
| 0.05 | 0.92774 | 0.92774 | 1.18849 |
| 0.1 | 0.86071 | 0.86071 | 1.71756 |
| 0.2 | 0.74082 | 0.74082 | 3.06594 |

The exact discretisation is stable at every step and reproduces e^{−ζω_n h} to 12 figures; forward Euler is stable only at h = 0.01 s (below 0.0133 s) and *amplifies* a decaying oscillation by 3.07 per step at h = 0.2 s. A 0.2 s step is 1/2 of the 0.42 s period, which is coarse for *resolving* the oscillation but perfectly fine for the ZOH map, which is exact at any step.

(d) ‖A·1 s‖₁ = max column sum = max(225, 1 + 3) = **225**. Bound ‖A‖^N/N! · e^{‖A‖}: N = 10 → 4.8 × 10¹¹⁴; 20 → 2.4 × 10¹²⁶; 30 → 7.2 × 10¹³⁵. Actual `expm_taylor` error against `expm`: **1.24 × 10⁶, 8.2 × 10⁵, 2.8 × 10³**, against true entries of order 2.4 (largest |entry| of e^{A}). The bound is astronomically loose because it bounds each term by its norm and ignores the cancellation that makes e^{A} small (the eigenvalues have |e^λ| = e^{−1.5} = 0.22); the actual error is nonetheless 10³–10⁶ against an answer of order 1 because the same cancellation, performed in floating point among terms of size up to 225^k/k! ≈ 10⁹⁶, leaves no significant figures, and 30 terms is also far short of the ≈ 225 needed for the terms even to start shrinking. `expm_scaling_squaring` uses s = ⌈log₂(225/0.5)⌉ = **9 squarings** and reaches 10⁻¹² of `expm`.

**Independent check.** The Euler limit by direct iteration: with h = 0.0133 the amplitude after 1000 steps is 0.9999 of the start; with h = 0.014 it has grown by a factor 1.4. e^{−1.5 × 0.05} = 0.927743 by calculator.

**Common mistakes.** Taking the Euler limit as h < 2/ω_n (that is the undamped-centre heuristic; the correct limit involves ζ and is much smaller here); reading |eig A_d| < 1 as "A_d is accurate" (it is exact, but a 0.2 s sample cannot *show* the oscillation); comparing the Taylor error with the bound and concluding the method "works".

## Check-yourself answers (lesson.md)

1. e^{A} = Σ_{k≥0} A^k/k!; ‖A^k/k!‖ ≤ ‖A‖^k/k! and Σ‖A‖^k/k! = e^{‖A‖} converges, so the matrix series converges absolutely for every A.
2. Φ(0) = I; Φ(t)Φ(s) = Φ(t + s); Φ⁻¹(t) = Φ(−t); dΦ/dt = AΦ; det Φ = e^{tr A·t}. Invertibility (the third, or the fifth: det never zero) means two different states cannot be mapped to the same state, so trajectories never cross or merge.
3. The series for e^{A+B} contains (A + B)² = A² + AB + BA + B², while e^Ae^B contains A² + 2AB + B² at that order; they agree only if AB = BA. Condition: A and B commute.
4. The eigen method: A_CW has eigenvalue 0 twice with one eigenvector (defective), so V is singular. Physically the defect is the secular drift (terms linear in t), which a sum of pure exponentials cannot represent.
5. x_{k+1} = A_d x_k + B_d u_k with A_d = e^{Ah}, B_d = ∫₀^h e^{As} ds B; e^{Mh} = [[A_d, B_d], [0, I]] for M = [[A, B], [0, 0]].
6. Δy = −12π x₀ = −12π(−100) = +3770 m: ahead of the target (the lower orbit is faster).
7. The position at T is linear in the post-burn velocity: Φ_rv Δv₁ = −(Φ_rr r₀ + Φ_rv v₀), a 2×2 system. No solution when det Φ_rv = 0, at T = one orbital period (and its multiples).
8. All eigenvalues must have strictly negative real part. A defective eigenvalue with zero real part gives solutions growing like t (or t^k), so the system is unstable even though no exponential grows: Clohessy–Wiltshire drift.
