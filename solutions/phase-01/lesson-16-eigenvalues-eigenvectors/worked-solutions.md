# P01.16  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for $s = 15$; numerical results are from `solution.py` and given to the figures shown.

## 1. A 2×2 tensor by hand

**Derivation** ($s = 15$: A = [[35, −4], [−4, 17.5]] kg m²). tr A = 52.5, det A = 35 × 17.5 − 16 = 596.5. Characteristic polynomial λ² − 52.5λ + 596.5 = 0; discriminant 52.5² − 4 × 596.5 = 2756.25 − 2386 = 370.25, √ = 19.242. λ = (52.5 ∓ 19.242)/2 = **16.629** and **35.871** kg m².

Eigenvector for λ₁ = 16.629: (A − λ₁I)**v** = 0 gives (35 − 16.629)v₁ − 4v₂ = 0, so **v**₁ ∝ (4, 18.371) ∝ (0.2128, 0.9771); angle atan2(0.9771, 0.2128) = 77.7° from the x axis. For λ₂ = 35.871: (35 − 35.871)v₁ − 4v₂ = 0 gives **v**₂ ∝ (−4, 0.871) → (−0.9771, 0.2128); angle 167.7° (equivalently −12.3°).

**Check.** **v**₁·**v**₂ = −0.2079 + 0.2079 = 0 ✓. λ₁ + λ₂ = 52.5 = tr ✓; λ₁λ₂ = 596.5 = det ✓. `jacobi_eigen` returns (16.6291, 35.8709) and `np.linalg.eigh` the same; the eigenvectors agree up to sign (numpy returns (−0.2128, −0.9771) for the first, which is the same line).

**Common mistakes.** Using a + d for the determinant; forgetting that an eigenvector's sign (and hence its angle ± 180°) is arbitrary; reporting the angle from atan instead of atan2.

## 2. Spectral theorem and triangle inequality

**Derivation.** (a) Let A**v** = λ**v**, A real symmetric, **v** ≠ 0 possibly complex. Multiply by **v*** (conjugate transpose): **v***A**v** = λ **v***v**. The scalar **v***A**v** satisfies (**v***A**v**)* = **v***A***v** = **v***A**v** (A real and symmetric gives A* = A), so it is real; **v***v** = Σ|v_i|² > 0 is real; hence λ is real.

(b) With A**v**₁ = λ₁**v**₁, A**v**₂ = λ₂**v**₂: **v**₂ᵀA**v**₁ = λ₁**v**₂ᵀ**v**₁, and also **v**₂ᵀA**v**₁ = (Aᵀ**v**₂)ᵀ**v**₁ = (A**v**₂)ᵀ**v**₁ = λ₂**v**₂ᵀ**v**₁. Subtracting, (λ₁ − λ₂)**v**₂ᵀ**v**₁ = 0, and λ₁ ≠ λ₂ forces **v**₂ᵀ**v**₁ = 0.

(c) **n**ᵀ[(**r**·**r**)I − **rr**ᵀ]**n** = (**r**·**r**)(**n**·**n**) − (**r**·**n**)² = |**r**|² − (**r**·**n**)². Decompose **r** = (**r**·**n**)**n** + **r**⊥; then |**r**|² = (**r**·**n**)² + |**r**⊥|², so the quadratic form equals |**r**⊥|². Summing with weights m_k gives Σ m_k|**r**_⊥,k|², the moment about the axis.

(d) I_xx + I_yy − I_zz = Σ m(y² + z²) + Σ m(x² + z²) − Σ m(x² + y²) = 2Σ m z² ≥ 0, with equality iff every mass has z = 0 (a flat body in the x–y plane). The proof used only the definition in an arbitrary frame; the principal frame is a frame, so I₁ + I₂ ≥ I₃ and its permutations hold for the principal moments.

**Check.** Satellite: 16.259 + 43.200 = 59.459 ≥ 52.515 ✓; smallest sum first.

**Common mistakes.** In (a), multiplying by **v**ᵀ instead of **v***, which does not give a positive real **v**ᵀ**v** for complex **v**; in (b), forgetting that the argument needs λ₁ ≠ λ₂.

## 3. Two-mass chain

**Assumptions.** Linear springs, no damping, small motion; $s = 15$: k = 2500 N/m, m₁ = 2 kg, m₂ = 1 kg.

**Derivation.** (a) With **x** = M^{−1/2}**y**: M M^{−1/2}**ÿ** + K M^{−1/2}**y** = 0; multiply by M^{−1/2}: **ÿ** + M^{−1/2}KM^{−1/2}**y** = 0. For **y** = **ŷ**e^{iωt}, S**ŷ** = ω²**ŷ** with S = M^{−1/2}KM^{−1/2}, symmetric because Sᵀ = M^{−1/2}KᵀM^{−1/2} = S. M⁻¹K has the same eigenvalues (it is M^{−1/2} S M^{1/2}, a similarity) but is not symmetric, so its eigenvectors are not orthogonal and Jacobi does not apply.

(b) M^{−1/2} = diag(0.70711, 1): S = [[2500, −1767.77], [−1767.77, 5000]] s⁻². `jacobi_eigen(S)` gives ω² = 1584.94 and 5915.06 s⁻², so ω = 39.811 and 76.910 rad/s, f = **6.336 Hz** and **12.241 Hz**.

(c) Mode shapes **x** = M^{−1/2}**ŷ**, normalised to x₁ = 1: mode 1 (6.34 Hz): (1, 0.7321); mode 2 (12.24 Hz): (1, −2.7321). The masses move in opposite directions in the higher mode, as always for a chain.

(d) K**x** − ω²M**x** evaluates to (0, 0) to 10⁻¹² for both modes (`np.allclose` true). M-orthogonality: **x**₁ᵀM**x**₂ = 2 × 1 × 1 + 1 × 0.7321 × (−2.7321) = 2 − 2.000 = 0 ✓.

**Check.** Limiting case: equal masses m₁ = m₂ = m give ω² = k/m and 3k/m with shapes (1, 1) and (1, −1); here the heavier m₁ lowers both frequencies and makes the light mass move more in the antisymmetric mode. Dimensions: k/m = N m⁻¹ kg⁻¹ = s⁻². ✓

**Common mistakes.** Taking the eigenvectors of S as the mode shapes without the M^{−1/2} factor (they are not M-orthogonal then... they are orthogonal, but they are **y**, not **x**); forgetting to divide ω by 2π.

## 4. A payload is added ($s = 15$: payload at (0.40, −0.30, 0.30) m)

**Result.** New total mass 128 kg, centre of mass (0.0219, 0.0063, 0.0428) m. Tensor about the new CM (kg m²):

$$
I = \begin{pmatrix} 45.195 & 0.158 & 2.276 \\ 0.158 & 17.589 & 0.676 \\ 2.276 & 0.676 & 53.004 \end{pmatrix},
$$

principal moments 17.575, 44.580, 53.632 kg m². (b) The major axis moved by 4.17° from its no-payload direction; it is now 15.15° from body z (was 19.13°). (c) `power_iteration` converges in 138 iterations; the ratio is 44.580/53.632 = 0.8312 and log(10⁻¹²)/log(0.8312) = 149: the law predicts within 8 %. (d) Closer: the payload's z-offset and its −y position partly cancel the antenna's product of inertia I_xz (2.888 → 2.276). The moment about body z is 53.004 kg m² against the major moment 53.632: 1.2 % below the maximum, consistent with a 15° tilt (cos² 15° = 0.933 weights toward the major moment, the rest toward the intermediate one).

**Check.** Trace 115.788 = sum of moments ✓. Triangle inequality 17.575 + 44.580 = 62.155 ≥ 53.632 ✓. The new tensor's I_zz is larger than before (53.004 against 51.513) because the payload adds 6 × (0.38² + 0.31²) ≈ 1.4 kg m² about z, minus the CM shift.

**Common mistakes.** Computing the new tensor about the *old* centre of mass; comparing eigenvectors without taking |**n**_old·**n**_new| (sign flips give 175° instead of 4°).

## 5. Spin stabilisation

(a) With **ω** = Ω**e**_z, **H** = I**ω** = Ω (2.888, 0.180, 51.513) kg m²/s: the angle between **H** and **ω** is arccos(51.513/√(2.888² + 0.180² + 51.513²)) = 3.22°. In torque-free flight **H** is fixed in inertial space, so **ω** and the body z axis cone around **H** at that angle: the satellite nutates with a 3.2° half-angle and its antenna sweeps a cone instead of pointing. (The 19° tilt of the principal axis produces only a 3.2° **H**–**ω** angle because the moments differ by less than 20 %.)

(b) Minor-axis spin is stable for a *rigid* body (the coefficient in the linearised Euler equation is positive for I₃ smallest) but any dissipation drives the body toward the minimum-energy state at fixed **H**, which is spin about the major axis. Explorer 1 (1958), spun about its long minor axis, went into a flat tumble within hours as its flexible whip antennas dissipated energy. A designer refuses because no satellite is perfectly rigid; the question is only the time scale.

(c) Two equal (or nearly equal) principal moments: `spin_axis_stability` raises `ValueError`. The linear rule has a zero coefficient there and says nothing; the axisymmetric-top solution of Euler's equations applies, with free precession at rate Ω(I₃ − I₁)/I₁ and the spin axis coning at a constant angle. Nearly equal moments also make the eigenvectors ill-conditioned: the "principal axes" in the degenerate plane are not meaningful.

(d) Tank at 16 kg: moments 16.259, 43.200, 52.515; relative major–intermediate gap (52.515 − 43.200)/52.515 = 17.7 %. Tank drained to 4 kg (replace the 16 in `satellite_point_masses` and recompute): the moments become 14.863, 42.706, 51.479 kg m², the major axis stays the major axis and the gap is 17.0 %. The tank sits close to the centre of mass, so it contributes little; the margin is set by the solar arrays (which fix the large I_xx and I_zz) and is robust to the propellant state.

**Common mistakes.** Equating the principal-axis tilt with the nutation angle; claiming the minor axis is "unstable" for a rigid body.

## 6. Damped system ($s = 15$: ζ = 0.15)

**Derivation.** (a) det(λI − A) = λ(λ + 2ζω) + ω² = λ² + 2ζωλ + ω². Roots λ = −ζω ± √(ζ²ω² − ω²) = −ζω ± iω√(1 − ζ²) for ζ < 1. With ω = 2, ζ = 0.15: λ = −0.3 ± 1.9774i s⁻¹.

(b) `linear_system_decay_rates(A)` returns (−0.3, −0.3) ✓. `power_iteration(A, max_iter=2000)` raises `RuntimeError`: both eigenvalues have modulus √(0.09 + 3.910) = 2 = ω, so the ratio |λ₂/λ₁| = 1 and the iterate rotates in the plane forever instead of settling (the real iterate traces the decaying oscillation's direction, which keeps turning).

(c) `expm(5A) @ (1, 0)` has norm 0.294 at t = 5 s, while e^{−ζω·5} = e^{−1.5} = 0.223: the ratio is 1.32. The justification for a constant above 1: **u**(t) = V e^{Λt} V⁻¹ **u**(0), so |**u**(t)| ≤ |V||V⁻¹| e^{−ζωt}|**u**(0)| = κ(V) e^{−ζωt}, where κ(V) is the condition number of the eigenvector matrix; with **u** = (x, ẋ) the eigenvectors are (1, λ) and (1, λ̄), far from orthogonal, and `np.linalg.cond` gives κ(V) = 2.04. The measured ratio 1.32 is below that bound, as it must be. (The lesson's test uses ζ = 0.1 and the looser constant 1.5.)

(d) `matrix_function_symmetric` calls `jacobi_eigen`, which raises `ValueError` on a non-symmetric input. Jacobi rotations are orthogonal similarities that preserve symmetry and drive a symmetric matrix to a *real* diagonal; A is not symmetric and its eigenvalues are complex, so no real orthogonal similarity can diagonalise it. `scipy.linalg.expm` (Padé approximation with scaling and squaring) or `numpy.linalg.eig` (complex eigenvectors) is the right tool.

**Check.** ζ → 0: λ = ±iω, pure oscillation, decay rate 0 ✓. ζ → 1: double root −ω, critical damping ✓.

**Common mistakes.** Expecting power iteration to "find the eigenvalue with the largest real part"; forgetting that the decay rate of **u** is governed by the real part of λ, not its modulus.

## Check-yourself answers

1. An eigenvector of A is a non-zero **v** with A**v** parallel to **v**; the eigenvalue λ is the factor A**v** = λ**v**. For the inertia tensor: a principal axis (spinning about it gives angular momentum along the spin) and the moment of inertia about it.
2. I₁ = trace, I₂ = sum of the principal 2×2 minors, I₃ = determinant; they are the coefficients of the characteristic polynomial, which is unchanged under A → CACᵀ (any change of frame), so they do not depend on the coordinates used to write A.
3. A**v** = λ**v** ⇒ **v***A**v** = λ**v***v**; the left side is real because A is real symmetric ((**v***A**v**)* = **v***A**v**), the right factor **v***v** > 0 is real; so λ is real.
4. **n**ᵀI**n** = Σ m[|**r**|² − (**r**·**n**)²] = Σ m|**r**⊥|². Its minimum over unit **n** is the smallest principal moment and its maximum the largest (attained on the corresponding principal axes).
5. Ratio r = |λ₂/λ₁| per step for the vector; steps ≈ log ε / log r = log(10⁻⁹)/log(0.9) = 197.
6. Torque-free spin about the axis of largest or smallest principal moment is stable, about the intermediate axis unstable. With energy dissipation only the largest-moment (major) axis is stable.
7. 16.26 + 43.20 = 59.46 ≥ 52.52, and the other two trivially; trace = 16.26 + 43.20 + 52.52 = 111.97 kg m²; determinant = 16.26 × 43.20 × 52.52 = 36 892 kg³ m⁶ with the rounded moments (36 885 with full precision: the product amplifies the rounding of each factor).
8. A state matrix is generally non-symmetric with complex eigenvalues; Jacobi requires symmetry and real orthogonal rotations. `eigvals` uses the QR algorithm, which handles the general case.
