# P06.07  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for **s = 15**; numerical results are from `solution.py` (and the probe scripts that call it) and are given to the figures shown. All "simulated" values are outputs of the model, not measurements.

## 1. An orbit as a 3-1-3 rotation (s = 15: Ω = 25°)

**Assumptions.** Passive DCM C_p/i; angles in radians in the code.

(a)

```
C = [[ 0.65363  0.64747  0.39185]
     [-0.68049  0.27622  0.67870]
     [ 0.33120 -0.71027  0.62115]]
```

orthogonality error **1.6 × 10⁻¹⁶**, det C = **1.000000000000000**.

(b) Row 3: **(0.33120, −0.71027, 0.62115)**. The formula gives sin 51.6° sin 25° = 0.33120, −sin 51.6° cos 25° = −0.71027, cos 51.6° = 0.62115. It is a row because v_p = C v_i: row i of C, dotted with a reference-frame vector, gives that vector's component along perifocal axis i, so row i *is* perifocal axis i in inertial components.

(c) 3-2-1: **(44.7288°, −23.0695°, 47.5351°)**; 3-1-3 back: **(25.0000°, 51.6000°, 30.0000°)**.

**Check.** The third row does not depend on ω, as it must not: the argument of periapsis turns the orbit within its own plane.

**Common mistakes.** Reading the normal from column 3 (that would be the inertial z axis in perifocal components); using R₃R₁R₃ in the wrong order.

## 2. Rates, determinants and the kinematic equation

(a) R₂(θ)e₃ = (−sin θ, 0, cos θ), so R₁(φ)R₂(θ)e₃ = (−sin θ, sin φ cos θ, cos φ cos θ); R₁(φ)e₂ = (0, cos φ, −sin φ); e₁. Hence **p = φ̇ − ψ̇ sin θ, q = θ̇ cos φ + ψ̇ sin φ cos θ, r = −θ̇ sin φ + ψ̇ cos φ cos θ**. Multiply q by cos φ and r by sin φ and subtract: **θ̇ = q cos φ − r sin φ**. Multiply q by sin φ and r by cos φ and add: q sin φ + r cos φ = ψ̇ cos θ, so **ψ̇ = (q sin φ + r cos φ)/cos θ**; then **φ̇ = p + ψ̇ sin θ = p + (q sin φ + r cos φ) tan θ**.

(b) det B = c₁ · (c₂ × c₃). The last two columns, c₃ = e_{s₃} and c₂ = R_{s₃}(a₃)e_{s₂}, are perpendicular unit vectors (s₂ ≠ s₃, and R_{s₃} leaves e_{s₃} alone), so n = c₂ × c₃ is a unit vector and det B is the cosine of the angle between c₁ and n. Multiplying all three columns by R_{s₃}(a₃)ᵀ changes neither dot nor cross products, so compare u = R_{s₂}(a₂)e_{s₁} with m = e_{s₂} × e_{s₃}. *Asymmetric:* e_{s₁} = ±m; the turn by a₂ about e_{s₂}, which is perpendicular to both, moves u away from m by a₂, so **|det B| = |cos a₂|**. *Symmetric:* e_{s₁} = e_{s₃}, which is perpendicular to m; the turn by a₂ tilts u towards ±m by a₂, so **|det B| = |sin a₂|**. For 3-2-1, det B = −cos θ, and at (30°, 10°, −5°) the code gives −0.984808.

(c) Differentiate CCᵀ = E: ĊCᵀ + CĊᵀ = 0, so W = ĊCᵀ satisfies W = −Wᵀ: it is skew-symmetric and can be written W = −[ω×]. Multiply on the right by C: **Ċ = −[ω×]C**. Then Ċᵀ = −Cᵀ[ω×]ᵀ = Cᵀ[ω×], so d(CᵀC)/dt = ĊᵀC + CᵀĊ = Cᵀ[ω×]C − Cᵀ[ω×]C = **0**: CᵀC keeps its initial value E.

(d) For 3-1-3 with (a₁, a₂, a₃) = (ψ, β, φ): R₁(β)e₃ = (0, sin β, cos β) and R₃(φ)(0, sin β, cos β) = (sin φ sin β, cos φ sin β, cos β); with β̇ = 0, ω_b = ψ̇ (sin β sin φ, sin β cos φ, cos β) + φ̇ e₃, which is **(6)**.

(e) With A = −[ω×] constant, k₁ = AC, k₂ = A(E + AΔt/2)C, k₃ = A(E + AΔt/2 + A²Δt²/4)C, k₄ = A(E + AΔt + A²Δt²/2 + A³Δt³/4)C; summing, C ← **(E + AΔt + (AΔt)²/2 + (AΔt)³/6 + (AΔt)⁴/24)C**, the series of exp(AΔt) cut after five terms. In the plane of rotation this multiplies by the complex number Σₙ₌₀⁴ (−ix)ⁿ/n! with x = |ω|Δt, whose squared magnitude is 1 − x⁶/72 + O(x⁸). So PᵀP − E has two entries of −x⁶/72 and **‖CᵀC − E‖ ≈ √2 x⁶/72 per step**. A probe gives 1.9617 × 10⁻⁸ at x = 0.1 against the estimate 1.964 × 10⁻⁸, and 3.068 × 10⁻¹⁰ at x = 0.05 against 3.07 × 10⁻¹⁰. Over a fixed time the defect accumulates to order Δt⁵, which Experiment 2 of the lesson observes.

**Common mistakes.** Using the active sign (+[ω×]); confusing ω of b relative to n with ω of n relative to b; forgetting that B's columns are rotated by the *later* turns only.

## 3. A rocket on the pad

(a) Rows are body axes in NED: **C = [[0, 0, −1], [0, 1, 0], [1, 0, 0]]**, equal to `dcm_from_euler([0, 90°, 0], "321")`.

(b) Singular (|det B| = 0): **1-2-3, 2-1-2, 2-3-2, 3-2-1**. The other eight have |det B| = 1.0 at this attitude.

(c) Through the pitch-over (3-2-1 pitch 90°, 60°, 30°, 0°): **1-3-2, 2-1-3, 2-3-1 and 3-1-2** keep |det B| = 1.0 throughout; 1-2-1, 1-3-1, 3-1-3 and 3-2-3 fall as 1.0, 0.866, 0.5, 0.0 and are singular when the vehicle is horizontal. Any of the four is a valid choice for this planar pitch-over; 2-1-3 is a common one.

(d) **No three-parameter set is non-singular everywhere**; a quaternion (four parameters, one constraint) is, so a vehicle whose attitude ranges widely carries a quaternion and computes angles only for display or telemetry.

**Check.** The four singular sequences are exactly those whose middle-axis turn is ±90° (asymmetric) or 0° (symmetric) at this attitude, as 2.3 predicts.

**Common mistakes.** Concluding that gimbal lock is avoided by choosing "a better sequence" for all attitudes; it is only moved.

## 4. Your own coning spacecraft (s = 15: ψ̇ = 0.35 rad/s, φ̇ = 3.0 rad/s, β = 17.5°)

(a) |ω| = **3.3355 rad/s**. `solve_ivp`, 60 s: largest misalignment **2.1 × 10⁻⁹ rad**, largest orthogonality error **1.3 × 10⁻⁹**.

(b)

| Δt (s) | misalignment (rad) | orthogonality error |
|---|---|---|
| 0.04 | 5.090 × 10⁻⁴ | 1.615 × 10⁻⁴ |
| 0.02 | 3.197 × 10⁻⁵ | 5.054 × 10⁻⁶ |
| 0.01 | 2.000 × 10⁻⁶ | 1.580 × 10⁻⁷ |
| 0.005 | 1.251 × 10⁻⁷ | 4.938 × 10⁻⁹ |

Misalignment orders **3.99, 4.00, 4.00**; orthogonality orders **5.0, 5.0, 5.0**.

(c) With renormalisation the orthogonality error is **4–13 × 10⁻¹⁶** at every step size, while the misalignment is unchanged to six figures (5.090 × 10⁻⁴ … 1.251 × 10⁻⁷ rad).

(d) 50 Hz is Δt = 0.02 s: **3.2 × 10⁻⁵ rad < 10⁻⁴ rad: yes**, with a margin of three.

**Check.** The orthogonality order is one more than the attitude order, as 2(e) predicts (a sixth-order defect per step).

**Common mistakes.** Measuring the error against the RK4 result at a finer step rather than the closed form (that hides any error common to both).

## 5. How fast must the computer run? (s = 15: δ = 0.20°)

(a) r/sin δ = 0.5/sin 0.2° = **143.2 rad/s**.

(b)

| rate (Hz) | Euler-angle RK4 (rad) | DCM RK4 (rad) |
|---|---|---|
| 50 | 2.02 × 10⁻³ | 2.6 × 10⁻¹⁰ |
| 100 | 8.74 × 10⁻⁵ | 1.6 × 10⁻¹¹ |
| 200 | 2.77 × 10⁻⁶ | 1.0 × 10⁻¹² |
| 400 | 1.74 × 10⁻⁷ | 6.3 × 10⁻¹⁴ |

(c) Orders **4.5, 5.0, 4.0**. The asymptotic RK4 order 4 appears only once ψ̇Δt is small (at 400 Hz, 143 × 0.0025 = 0.36); at 50 Hz ψ̇Δt ≈ 2.9 and the error is not in its asymptotic range, so the apparent order is erratic.

(d) **100 Hz** for the Euler angles (8.7 × 10⁻⁵ rad); the DCM meets 10⁻⁴ rad at **50 Hz** with seven orders of magnitude to spare, and would do so at much lower rates, limited only by the body rates themselves (0.5 rad/s).

**Common mistakes.** Running RK4 straight through the rate switches; the kink then dominates both errors (of order 10⁻³ rad at 100 Hz in a probe run), which hides the difference this problem is about.

## 6. Judgment

(a) No. It should integrate the DCM (or a quaternion), which is non-singular and cheap, and compute yaw, pitch and roll from it only when a display or a human needs them; then a singularity costs a jumpy display, not a corrupted attitude.

(b) A wrong parity sign ε in the extraction: 3-2-1 is an anti-cyclic order (ε = −1) and 1-2-3 a cyclic one (ε = +1), so code that hard-wires the 3-2-1 signs works for 3-2-1 and its anti-cyclic relatives only.

(c) Re-orthonormalisation makes the matrix *a* rotation, the one nearest to the integrated matrix. It cannot know the true attitude, so the integration error in attitude remains (the tests show it unchanged to five per cent).

(d) Its matrices map vectors from the rotated frame into the original frame, the transpose of this lesson's C_b/n. Compare its tables with Cᵀ (and check the angle labels and sequence order) before deciding that either is wrong; the scan also has OCR noise, so check a few entries by multiplying out.

(e) Gimbal lock is not physical for a strapdown system: the body can rotate any way at any attitude, and only the three-angle coordinates are singular. In a mechanical three-gimbal platform it *is* physical: when two gimbal axes align, the platform loses a rotational degree of freedom and can be tumbled, which is why such platforms carried a fourth gimbal or avoided the attitude.

(f) A 30° margin keeps |det B| ≥ 0.5, so the angle rates are at most twice the body rates and integrating 3-2-1 would work for the nominal profile. The DCM is still the better choice: off-nominal events (a tumble, a failed thruster, a recovery manoeuvre) are exactly when the attitude must be right, and the DCM costs only a few more operations.

## Check-yourself answers

1. Twelve. A sequence like 3-3-1 makes two successive turns about the same axis, which is one turn; it has only two independent angles and cannot reach every attitude.
2. The first rotation is R_{s₁}(a₁), and it is the rightmost factor.
3. 3-2-1: θ = ±90°. 3-1-3: a₂ = 0° or 180°.
4. |det B| = |cos a₂|: the triple product of the three rotation axes, which measures how far the first axis leans out of the plane of the other two; zero means the three axes are coplanar.
5. CCᵀ = E for all t ⇒ ĊCᵀ + CĊᵀ = 0 ⇒ ĊCᵀ = −(ĊCᵀ)ᵀ.
6. It restores CᵀC = E and det C = 1 (the matrix is a rotation again); it does not remove the attitude error.
7. The attitude error by 2⁴ = 16; the orthogonality drift by 2⁵ = 32.
8. No: near θ = ±90° the 3-2-1 angles of nearby attitudes differ by up to 180° in yaw and roll. The DCM is continuous; the angle coordinates are not.
