# Phase 01 Exam: worked solutions and pass thresholds

Attempt the exam first, in the time you set yourself. Then grade each question against this key and record the result per competency (1 reasoning, 2 analysis, 3 computation, 4 interpretation). **Pass threshold:** proficient (3 of 4) or better on every competency, where "proficient" means the reasoning is right and the result is right to within rounding, even if the presentation is rough. In Part B a correct final number reached by a wrong derivation is *not* proficient: the derivation is what the question tests. Every number below was reproduced with NumPy and SymPy before it was written down.

## Part A

**1.** A forward difference (f(x+h) − f(x))/h differs from f′(x) by two errors. The *truncation* error comes from the terms of the Taylor series that the quotient drops; for the forward difference it is ≈ (h/2) f″(x), proportional to h, and for the central difference (f(x+h) − f(x−h))/(2h) it is ≈ (h²/6) f‴(x), proportional to h². The *rounding* error comes from subtracting two nearly equal numbers that each carry a relative error of about ε = 2.2 × 10⁻¹⁶; the difference has absolute error ≈ ε|f|, and dividing by h gives ≈ ε|f|/h, which *grows* as h shrinks. The two balance where h ≈ √ε ≈ 10⁻⁸ for the forward difference and h ≈ ε^{1/3} ≈ 6 × 10⁻⁶ for the central difference; below that the answer gets worse. (Check: for sin at x = 1, the forward-difference error is 3 × 10⁻⁹ at h = 10⁻⁸ but 4 × 10⁻⁵ at h = 10⁻¹²; the central difference reaches 3 × 10⁻¹¹ at h = 10⁻⁶.) "Matches exactly" is therefore a statement about a tolerance that was not stated; the honest sentence is "agrees with the symbolic derivative to within 10⁻⁸ relative at h = 10⁻⁸, as the error model predicts". *Proficient:* both error sources named with their scaling in h, the order of magnitude √ε for the forward difference, and the conclusion that smaller h is not better.

**2.** e = 0: circle of radius p. 0 < e < 1: ellipse with the focus at the origin. e = 1: parabola (open, escape at zero excess speed). e > 1: hyperbola. The radius is smallest where the denominator is largest, at θ = 0 (periapsis), r_p = p/(1 + e). Constancy of h: d(r × v)/dt = v × v + r × a = 0 + r × a, and r × a = 0 because a is parallel to r for a central force; so h is constant. The rate at which r sweeps area is dA/dt = ½|r × v| = ½|h|, so equal areas are swept in equal times (Kepler's second law) for *any* central force, not only inverse square. *Proficient:* four shapes, periapsis at θ = 0, the two-line cross-product argument, and the areal-rate conclusion.

**3.** rank A = 7 < 8 means the 8 unknowns (member forces and reactions) are not all determined: the nullspace is one-dimensional and contains a non-zero set of internal forces and reactions that is in equilibrium with *zero* applied load (a state of self-stress). By rank–nullity for the square 8 × 8 matrix, the left nullspace is also one-dimensional: there is a load vector direction y with A^T y = 0 that no combination of member forces can balance, so the structure has a mechanism (a rigid-body motion of part of the truss that no member resists) in that direction. Counting m + r = 2j is necessary, not sufficient: the members are arranged so that one is redundant where another is missing. The designer should look for a joint whose members are collinear (no stiffness perpendicular to them), a panel without a diagonal, or reactions that all pass through one point (free rotation about it); then move a member or a support so that the rank becomes 8. *Proficient:* nullspace = self-stress, left nullspace = unsupportable load / mechanism, rank–nullity invoked, and a concrete geometric cause.

**4.** C^T C = I (orthogonal: it preserves lengths and angles) and det C = +1 (proper: it preserves handedness). det C = −1 is an orthogonal matrix that includes a reflection; it would turn a right-handed body frame into a left-handed one, which no rotation of a rigid body can do. Multiplying many rotations in floating point commits a rounding error of order ε at each step; the errors accumulate so that C^T C drifts from I and det C from 1, and the columns slowly stop being orthonormal. Restore them with the orthogonal factor of a polar decomposition via the SVD, C ← U V^T (lesson 18), or with Gram–Schmidt / QR on the columns (lesson 17). *Proficient:* both properties, the reflection interpretation, the accumulation argument, one named fix.

## Part B

**5.** (a) By definition, d(sin α)/dα = lim_{h→0} [sin(α + h) − sin α]/h. The angle-sum identity gives sin(α + h) = sin α cos h + cos α sin h, so the quotient is sin α (cos h − 1)/h + cos α (sin h)/h. With the two given limits, the first term → sin α · 0 and the second → cos α · 1, hence cos α. (SymPy: `limit((sin(a+h)-sin(a))/h, h, 0)` = cos(a).)
(b) dC_L/dα = 2π cos α, so at α = 0 the lift-curve slope is 2π per radian ≈ 6.28 per radian (0.110 per degree).
(c) Multiply by (1 + cos α)/(1 + cos α): (1 − cos²α)/(α²(1 + cos α)) = (sin α/α)² · 1/(1 + cos α) → 1² · ½ = ½. Check at α = 0.01: (1 − cos 0.01)/10⁻⁴ = 0.499996, which is ½ to five figures (the deviation is the α²/12 term of the series). *Proficient:* a definition-based derivation in (a), 2π in (b), ½ with a numeric check in (c).

**6.** L/D = C_L/(C_D0 + K C_L²). Quotient rule: d(L/D)/dC_L = [(C_D0 + K C_L²) − C_L · 2K C_L]/(C_D0 + K C_L²)² = (C_D0 − K C_L²)/(C_D0 + K C_L²)². Zero when K C_L² = C_D0, so C_L* = √(C_D0/K) = √0.4 = 0.632. The numerator changes sign from + to − as C_L passes C_L* (and the second derivative there is −39.5 in these units), so it is a maximum. At C_L* the two drag terms are equal by construction (K C_L*² = C_D0), so C_D = 2C_D0 = 0.040 and (L/D)_max = C_L*/(2C_D0) = 0.632/0.040 = 15.8; in closed form 1/(2√(C_D0 K)) = 1/(2 × 0.0316) = 15.8. Glide range from 1000 m in still air: 15.8 km. Airspeed at C_L* for 300 kg, 12 m², sea level: V = √(2mg/(ρ S C_L*)) = √(2 × 300 × 9.81/(1.225 × 12 × 0.632)) = √(5886/9.297) = 25.2 m/s. *Proficient:* the differentiation, a maximum argument, C_L* = 0.632, L/D = 15.8, range 15.8 km, the equal-split statement; the speed is a bonus check.

**7.** (a) Substitute y = (b/2) sin φ, dy = (b/2) cos φ dφ, so √(1 − (2y/b)²) = cos φ and the limits become ±π/2: ∫L′ dy = L₀ (b/2) ∫_{−π/2}^{π/2} cos²φ dφ = L₀ (b/2)(π/2) = π b L₀/4. For the root moment put u = 1 − (2y/b)², du = −(8y/b²) dy, so y dy = −(b²/8) du and M = L₀ (b²/8) ∫₀¹ √u du = L₀ (b²/8)(2/3) = L₀ b²/12. Numbers: total lift = π × 10 × 400/4 = 3142 N, so m = 3142/9.81 = 320 kg; root moment from one half-wing = 400 × 100/12 = 3333 N·m. (Cross-check: half the lift, 1571 N, acting at the centroid of the half-ellipse 2b/(3π) = 2.12 m from the root gives 1571 × 2.12 = 3333 N·m.)
(b) The region is 0 ≤ x ≤ c, 0 ≤ y ≤ s x/c. Area A = ∫₀^c ∫₀^{sx/c} dy dx = ∫₀^c (s x/c) dx = cs/2. First moments: ∫∫x dA = ∫₀^c (s x²/c) dx = c² s/3, so x̄ = 2c/3; ∫∫y dA = ∫₀^c ½ (s x/c)² dx = c s²/6, so ȳ = s/3. Numbers: A = 0.36 m², x̄ = 0.800 m, ȳ = 0.200 m, mass = 4 × 0.36 = 1.44 kg. *Proficient:* both substitutions carried through in (a) with 3142 N and 3333 N·m; correct limits of integration and (2c/3, s/3) in (b).

**8.** sin α = α − α³/6 + α⁵/120 − … . Taking the approximation sin α ≈ α as the Taylor polynomial of degree 2 (the α² coefficient is zero), the Lagrange remainder is R₂ = f‴(ξ) α³/3! with f‴ = −cos ξ, so |R₂| ≤ |α|³/6. At 10° = 0.17453 rad the bound is 0.17453³/6 = 8.86 × 10⁻⁴; the actual error is 0.174533 − 0.173648 = 8.85 × 10⁻⁴ (the bound is tight because cos ξ ≈ 1 over the interval). Relative error at 10°: 8.85 × 10⁻⁴/0.1736 = 0.51%. The bound reaches 1% relative error when α³/6 ≈ 0.01 α, i.e. α² = 0.06, α = 0.245 rad = 14.0°. For cos α ≈ 1 − α²/2 the next term is α⁴/24 with |f⁗| = |cos ξ| ≤ 1, so the error at 10° is at most 0.17453⁴/24 = 3.87 × 10⁻⁵ (actual 3.86 × 10⁻⁵). *Proficient:* the remainder written with the correct derivative bound, 8.86 × 10⁻⁴ vs 8.85 × 10⁻⁴, 0.5%, about 14°, and 3.9 × 10⁻⁵ for the cosine.

**9.** (a) ∂L/∂α = ½ ρ V² S C_Lα = 0.5 × 1.225 × 2500 × 16 × 5.7 = 139 650 N/rad = 2437 N/deg; ∂L/∂V = ρ V S C_Lα α = 1.225 × 50 × 16 × 5.7 × 0.069813 = 390.0 N per m/s. For reference L = 139 650 × 0.069813 = 9749 N.
(b) ln L = const + 2 ln V + ln α, so dL/L = dα/α + 2 dV/V. A 2% increase in α gives +2% of L (195 N); a 2% increase in V gives +4% of L (390 N). Speed matters twice as much, ratio 2.
(c) Linearised: 2% + 4% = 6.0% (585 N). Exact: (1.02)² × 1.02 − 1 = 1.061208 − 1 = 6.12% (597 N). The 0.12% difference is the second-order term (the cross term 2 dα dV/(αV) and dV²/V²); the linear model is adequate for sensitivities of this size and the gap grows with the square of the perturbation. *Proficient:* both partials correct with the per-degree conversion, the relative-differential derivation, ratio 2, and 6.0% vs 6.12% with the reason.

**10.** (a) Speed = √(u² + v²) = Γ √(x² + y²)/(2π(x² + y²)) = Γ/(2πr). At R = 0.5 m: 12/(2π × 0.5) = 3.82 m/s.
(b) Parametrise x = R cos t, y = R sin t, dx = −R sin t dt, dy = R cos t dt. Then u dx + v dy = [Γ/(2πR²)] (−R sin t)(−R sin t) dt + [Γ/(2πR²)] (R cos t)(R cos t) dt = (Γ/2π)(sin²t + cos²t) dt = (Γ/2π) dt, and the integral from 0 to 2π is Γ = 12 m²/s, independent of R.
(c) With r² = x² + y²: ∂v/∂x = (Γ/2π)(r² − 2x²)/r⁴ and −∂u/∂y = (Γ/2π)(r² − 2y²)/r⁴, so ∂v/∂x − ∂u/∂y = (Γ/2π)(2r² − 2x² − 2y²)/r⁴ = 0 for r ≠ 0. Green's theorem ∮(u dx + v dy) = ∬(∂v/∂x − ∂u/∂y) dA requires the field to be differentiable throughout the enclosed region; the circle encloses the origin, where the field is singular, so the theorem does not apply and there is no contradiction: all the vorticity is concentrated at the origin (a delta function of strength Γ). The square 1 ≤ x ≤ 2, 1 ≤ y ≤ 2 does not contain the origin, the field is smooth inside it, and Green's theorem gives circulation 0. (Numerical check: 2 × 10⁻¹⁶ for the square, 12.0 for a square enclosing the origin.) *Proficient:* the parametrised integral giving Γ independent of R, the curl computed to zero, the singularity named as the reason the theorem does not apply, and 0 for the off-centre square.

**11.** Augmented matrix and elimination (multipliers written as "row ← row + m × pivot row"):

```
[ 2   1  −1 |   8 ]
[−3  −1   2 | −11 ]   R2 ← R2 + 1.5 R1
[−2   1   2 |  −3 ]   R3 ← R3 + 1.0 R1

[ 2   1   −1  |  8 ]
[ 0   0.5  0.5 |  1 ]
[ 0   2    1   |  5 ]   R3 ← R3 − 4 R2

[ 2   1   −1  |  8 ]
[ 0   0.5  0.5 |  1 ]
[ 0   0   −1   |  1 ]
```

Pivots 2, 0.5, −1; determinant = product of pivots = −1 (no row exchanges, so no sign change). Back-substitution: −F₃ = 1 → F₃ = −1; 0.5F₂ + 0.5(−1) = 1 → F₂ = 3; 2F₁ + 3 + 1 = 8 → F₁ = 2. So (F₁, F₂, F₃) = (2, 3, −1) kN, the negative sign meaning F₃ acts opposite to the assumed direction (compression if tension was assumed). Check: 4 + 3 + 1 = 8 ✓; −6 − 3 − 2 = −11 ✓; −4 + 3 − 2 = −3 ✓. (NumPy: `solve` gives [2, 3, −1], `det` gives −1.0.) *Proficient:* the three multipliers, the pivots, det = −1, the solution, and the substitution check.

**12.** Characteristic polynomial: det(I − λ1) = (5 − λ)(2 − λ) − 4 = λ² − 7λ + 6 = (λ − 6)(λ − 1), so the principal moments are 6 and 1 kg·m². Checks: trace 5 + 2 = 7 = 6 + 1 ✓; det 10 − 4 = 6 = 6 × 1 ✓. Eigenvectors: for λ = 6, (5 − 6)x − 2y = 0 → x = −2y → (2, −1)/√5; for λ = 1, (5 − 1)x − 2y = 0 → y = 2x → (1, 2)/√5. Dot product 2 × 1 + (−1) × 2 = 0, orthogonal; guaranteed because I is real symmetric (spectral theorem), which is why principal axes are always a right-angled set. The axis of maximum inertia (6 kg·m²) is along (2, −1), at atan2(−1, 2) = −26.6° from the body x-axis; equivalently the formula tan 2θ = 2I_xy/(I_xx − I_yy) with I_xy = −2 gives θ = ½ atan2(−4, 3) = −26.6°. The minimum-inertia axis is at +63.4°. (NumPy `eigh`: eigenvalues [1, 6], and Q^T I Q with Q = [(1,2),(2,−1)]/√5 gives diag(1, 6).) *Proficient:* λ = 6 and 1, both checks, orthogonal unit eigenvectors with the symmetry reason, angle −26.6° (or 153.4°, the same axis).

**13.** Let x = C_L² = (0.04, 0.16, 0.36, 0.64) and y = C_D. Design matrix columns (1, x). Sums: n = 4, Σx = 1.20, Σx² = 0.5664, Σy = 0.140, Σxy = 0.04 × 0.023 + 0.16 × 0.027 + 0.36 × 0.039 + 0.64 × 0.051 = 0.00092 + 0.00432 + 0.01404 + 0.03264 = 0.05192. Normal equations:

```
[ 4     1.20   ] [C_D0]   [0.140  ]
[ 1.20  0.5664 ] [ K  ] = [0.05192]
```

Determinant 4 × 0.5664 − 1.44 = 0.8256. K = (4 × 0.05192 − 1.20 × 0.140)/0.8256 = (0.20768 − 0.168)/0.8256 = 0.03968/0.8256 = 0.04806; C_D0 = (0.140 − 0.04806 × 1.20)/4 = 0.02058. Fitted C_D: 0.02250, 0.02827, 0.03788, 0.05134. Residuals r = y − ŷ = (+0.00050, −0.00127, +0.00112, −0.00034), ‖r‖ = 0.0018. Orthogonality: Σr = 0.00050 − 0.00127 + 0.00112 − 0.00034 = 0.00001 ≈ 0 (rounding), Σ x r = 0.04 × 0.00050 − 0.16 × 0.00127 + 0.36 × 0.00112 − 0.64 × 0.00034 = 0.00002 − 0.000203 + 0.000403 − 0.000218 ≈ 0; with unrounded residuals both are zero to 10⁻¹⁷, which is what the normal equations A^T(b − Aθ) = 0 state. (L/D)_max = 1/(2√(C_D0 K)) = 1/(2√(0.02058 × 0.04806)) = 1/(2 × 0.03145) = 15.9, at C_L* = √(C_D0/K) = 0.654. (NumPy `lstsq` and QR both give θ = [0.020581, 0.048062].) *Proficient:* the four sums, the 2 × 2 solve, C_D0 = 0.0206 and K = 0.0481, residuals that sum to zero against both columns, (L/D)_max ≈ 15.9.

## Part C

**14.** A proficient answer looks like:

```python
import numpy as np
import pytest


def power_iteration(A: np.ndarray, n_iter: int = 200) -> tuple[float, np.ndarray]:
    """Dominant eigenpair of a real symmetric matrix by power iteration.

    Returns (lambda, v) with v a unit vector and lambda the Rayleigh
    quotient v^T A v at the final iterate.
    """
    A = np.asarray(A, dtype=float)
    v = np.ones(A.shape[0])
    v /= np.linalg.norm(v)
    for _ in range(n_iter):
        w = A @ v
        v = w / np.linalg.norm(w)
    lam = float(v @ A @ v)
    return lam, v


INERTIA = np.array([[6.0, -1.0, 0.0],
                    [-1.0, 4.0, -1.0],
                    [0.0, -1.0, 3.0]])


def test_matches_numpy_oracle():
    lam, v = power_iteration(INERTIA)
    w, V = np.linalg.eigh(INERTIA)          # ascending eigenvalues
    assert lam == pytest.approx(w[-1], rel=1e-8)
    assert abs(v @ V[:, -1]) == pytest.approx(1.0, abs=1e-8)   # same axis, either sign


def test_eigen_residual_is_small():
    lam, v = power_iteration(INERTIA)
    assert np.linalg.norm(INERTIA @ v - lam * v) < 1e-8
    assert np.linalg.norm(v) == pytest.approx(1.0)


def test_identity_limiting_case():
    lam, v = power_iteration(np.eye(3), n_iter=5)
    assert lam == pytest.approx(1.0)
    assert np.linalg.norm(v) == pytest.approx(1.0)
```

For the matrix above the iteration returns λ = 6.4605 with residual 1 × 10⁻¹⁶, and `eigh` gives eigenvalues 2.3004, 4.2391, 6.4605. Why residual and oracle rather than a hard-coded number: the residual ‖Av − λv‖ is the *definition* of an eigenpair, so it tests the property the function promises and stays valid if the matrix in the test is changed; the oracle comparison guards against a method that converges to the wrong eigenpair while still having a small residual (the identity case shows every vector is an eigenvector). Convergence rate: the error shrinks by the ratio |λ₂/λ₁| per iteration (here 4.239/6.461 = 0.656), so a dominant eigenvalue that is not well separated converges slowly. *Proficient:* a working function with normalisation and Rayleigh quotient, the three tests present and passing, the residual/oracle explanation, and the |λ₂/λ₁| ratio.

## Part D

**15.** (a) Load uncertainty: ‖δx‖/‖x‖ ≤ 5 × 10⁴ × 10⁻⁴ = 5, i.e. up to 500% in the worst case. Rounding: 5 × 10⁴ × 2.2 × 10⁻¹⁶ ≈ 1 × 10⁻¹¹, negligible. The data uncertainty dominates by eleven orders of magnitude; the solver's precision is irrelevant to the quality of the answer, and six significant figures in a quantity that might be wrong by a factor of five is not justified. The honest statement would report the displacements with an uncertainty obtained by propagating the load uncertainty through K (or by re-solving with perturbed loads).
(b) The bound is attained when the true load lies along the singular direction of K with the *largest* singular value (so ‖x‖ is as small as possible relative to ‖b‖) and the load error lies along the direction with the *smallest* singular value σ_min (so it is amplified by 1/σ_min). The smallest singular value is the stiffness of the softest deformation mode of the structure: a nearly free mechanism, a long slender member in bending, or a node that is almost unconstrained. A load error in that direction produces a large displacement for little force. For a typical load the actual error is far below the bound; the SVD tells you *which* loads are dangerous.
(c) Modelling causes: (i) mixed units or scales, for example rotational degrees of freedom in N·m/rad beside translational ones in N/m, or a model that mixes millimetres and metres, which inflates κ without changing the physics; a diagonal rescaling (equilibration) often removes most of it. (ii) A very stiff element used to approximate a rigid connection or a penalty constraint, whose stiffness is many orders larger than the structure's, which pushes σ_max up; replace it with a true constraint (eliminate the degree of freedom). Also a nearly unconstrained degree of freedom from a missing or badly placed support (which pushes σ_min toward zero). Check first the units and scaling of each degree-of-freedom type, then look at the singular vector of σ_min to see which nodes move in the soft mode. *Proficient:* 500% vs 10⁻¹¹ with the conclusion about significant figures, the σ_min direction identified as the soft mode, two modelling causes and the first check.
