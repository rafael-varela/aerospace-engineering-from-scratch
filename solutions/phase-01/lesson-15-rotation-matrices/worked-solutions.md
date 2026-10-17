# P01.15  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for $s = 15$; numerical results are from `solution.py` and given to the figures shown.

## 1. One attitude

**Assumptions.** Passive 3-2-1 convention of the lesson; angles ψ = 45°, θ = 5°, φ = −7.5° for $s = 15$.

**Derivation.** R₃(45°) has top-left block (0.70711, 0.70711; −0.70711, 0.70711). R₂(5°) has (1,1) = (3,3) = 0.99619, (1,3) = −0.08716, (3,1) = 0.08716. R₁(−7.5°) has (2,2) = (3,3) = 0.99144, (2,3) = −0.13053, (3,2) = 0.13053. The product R₁R₂R₃ is

$$
C_{b/n} = \begin{pmatrix} 0.7044 & 0.7044 & -0.0872 \\ -0.7091 & 0.6930 & -0.1300 \\ -0.0312 & 0.1534 & 0.9877 \end{pmatrix}.
$$

(b) Row 1: 0.7044² + 0.7044² + 0.0872² = 0.4962 + 0.4962 + 0.0076 = 1.0000 ✓. Row 1 · row 3 = −0.0220 + 0.1081 − 0.0861 = 0.0000 ✓. Determinant along row 1: 0.7044(0.6930 × 0.9877 − (−0.1300)(0.1534)) − 0.7044((−0.7091)(0.9877) − (−0.1300)(−0.0312)) + (−0.0872)((−0.7091)(0.1534) − 0.6930(−0.0312)) = 0.7044 × 0.7044 − 0.7044 × (−0.7044) − 0.0872 × (−0.0872) = 0.4962 + 0.4962 + 0.0076 = 1.0000 ✓.

(c) θ = −arcsin(−0.08716) = 5.0°; ψ = atan2(0.7044, 0.7044) = 45.0°; φ = atan2(−0.1300, 0.9877) = −7.5°. `euler_321_from_dcm` returns (45.0, 5.0, −7.5)°.

(d) Gravity in body axes: g × column 3 = (−0.855, −1.275, 9.686) m/s². Axis-angle: angle 46.17° about n = (−0.196, 0.039, 0.980).

**Result.** The matrix above; det = 1.000000; recovered angles exact to 10⁻¹³ degrees in the code.

**Check.** |g_b| = √(0.855² + 1.275² + 9.686²) = 9.807 m/s² ✓ (length preserved). cos 46.17° = 0.6928 and (tr C − 1)/2 = (0.7044 + 0.6930 + 0.9877 − 1)/2 = 0.6926 ✓.

**Common mistakes.** Multiplying R₃R₂R₁ (gives C₁₃ = −sin θ cos φ cos ψ..., not −sin θ); using degrees inside `np.cos`; reading the gravity vector from the third *row* instead of the third column.

## 2. Group properties and gimbal lock

**Derivation.** (a) (C₁C₂)ᵀ(C₁C₂) = C₂ᵀ(C₁ᵀC₁)C₂ = C₂ᵀC₂ = I; det(C₁C₂) = det C₁ det C₂ = 1. For the transpose: (C₁ᵀ)ᵀC₁ᵀ = C₁C₁ᵀ = I because a square matrix with CᵀC = I also has CCᵀ = I (Cᵀ is a left inverse, hence the two-sided inverse); det C₁ᵀ = det C₁ = 1. So rotations are closed under products and inverses: they form a group.

(b) (Cu)·(Cv) = (Cu)ᵀ(Cv) = uᵀCᵀCv = uᵀv = u·v. With u = v: |Cv|² = |v|².

(c) At θ = 90°, sin θ = 1, cos θ = 0. Row 2 becomes (sin φ cos ψ − cos φ sin ψ, sin φ sin ψ + cos φ cos ψ, 0) = (sin(φ − ψ), cos(φ − ψ), 0). Row 3 becomes (cos φ cos ψ + sin φ sin ψ, cos φ sin ψ − sin φ cos ψ, 0) = (cos(φ − ψ), −sin(φ − ψ), 0). Row 1 is (0, 0, −1). The matrix depends on φ and ψ only through φ − ψ: all attitudes with the same φ − ψ are the same attitude. The frame's x axis points straight down, and yawing about NED-down and rolling about body-x (which is also NED-down) are the same physical motion. At θ = −90°: row 2 = (−sin(φ + ψ), cos(φ + ψ), 0), row 3 = (cos(φ + ψ), sin(φ + ψ), 0); only φ + ψ is determined.

(d) C₁₁ = cos θ cos ψ and C₁₂ = cos θ sin ψ with cos θ > 0 for |θ| < 90°: atan2(C₁₂, C₁₁) = atan2(sin ψ, cos ψ) = ψ, including quadrant, since atan2 is invariant under a positive common scale of its arguments. Likewise C₂₃ = sin φ cos θ, C₃₃ = cos φ cos θ give φ. θ = −arcsin(C₁₃) is unique in [−90°, 90°].

**Check.** Numerically, `dcm_from_euler_321(0.8, pi/2, -0.5)` and `dcm_from_euler_321(0.3, pi/2, -1.0)` differ by less than 10⁻¹⁶: both have φ − ψ = −1.3 rad.

**Common mistakes.** Proving CᵀC = I and assuming CCᵀ = I needs no argument (it does, via the square-matrix inverse); claiming gimbal lock is a singularity of the *rotation*; forgetting the quadrant when "solving" ψ = atan(C₁₂/C₁₁).

## 3. Star tracker

**Assumptions.** Rigid mounting; the tracker's reported attitude is exact.

**Derivation.** (a) `dcm_from_axis_angle([1, 1, 1], 120°)` returns

$$
C_{t/b} = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 1 & 0 & 0 \end{pmatrix},
$$

exactly (the entries are 0 and 1 to 10⁻¹⁶). Row 1 says tracker-x is body-y; row 2 says tracker-y is body-z; row 3 says tracker-z is body-x. A 120° turn about the diagonal permutes the axes cyclically.

(b) C_b/i = C_b/t C_t/i = C_t/bᵀ C_t/i. With C_t/i = `dcm_from_euler_321(40°, 20°, 10°)`:

$$
C_{t/i} = \begin{pmatrix} 0.71985 & 0.60402 & -0.34202 \\ -0.58753 & 0.79258 & 0.16318 \\ 0.36964 & 0.08348 & 0.92542 \end{pmatrix}, \qquad
C_{b/i} = \begin{pmatrix} 0.36964 & 0.08348 & 0.92542 \\ 0.71985 & 0.60402 & -0.34202 \\ -0.58753 & 0.79258 & 0.16318 \end{pmatrix}.
$$

The rows of C_t/i are permuted, as they must be: a row is a tracker axis in inertial coordinates, and body-x is tracker-z. The other order, C_t/i C_t/bᵀ, would apply the mounting rotation to inertial components, which have nothing to do with how the tracker is bolted on; its Euler angles come out as (−29.5°, −46.0°, −57.8°), a different attitude.

(c) Thrust along body +x in inertial axes: C_b/iᵀ (1, 0, 0) = first row of C_b/i = (0.370, 0.083, 0.925). The thruster pushes mostly along inertial z.

(d) `euler_321_from_dcm(C_b/i)` = (12.73°, −67.73°, −64.49°). These bear no simple relation to (40°, 20°, 10°) plus 120°: Euler angles of a composed rotation are not sums of angles, because the elementary rotations about different axes do not commute. Only rotations about one common axis add.

**Check.** C_b/i is a rotation: rows orthonormal to 10⁻¹⁶, det = 1, `is_rotation_matrix` passes. The thrust direction has unit length (0.370² + 0.083² + 0.925² = 1.000), and the body-x row of C_b/i equals the tracker-z row of C_t/i, which is what the mounting says it must.

**Common mistakes.** Using C_t/b instead of its transpose; taking the thrust direction from a column instead of a row (transposes again); expecting the body pitch to be 20°.

## 4. The integration experiment ($s = 15$: ω = 15°/s)

**Result** (axis (1, 2, 2)/3, T = 20 s, total rotation 300°):

| N | Δt [s] | orth. error raw | misalignment raw [rad] | orth. error renorm. | misalignment renorm. [rad] | prediction [rad] |
|---|---|---|---|---|---|---|
| 200 | 0.1 | 1.305 × 10⁻¹ | 1.196 × 10⁻³ | 3.1 × 10⁻⁷ | 1.196 × 10⁻³ | 1.196 × 10⁻³ |
| 400 | 0.05 | 6.305 × 10⁻² | 2.990 × 10⁻⁴ | 2.0 × 10⁻⁸ | 2.990 × 10⁻⁴ | 2.991 × 10⁻⁴ |
| 800 | 0.025 | 3.099 × 10⁻² | 7.476 × 10⁻⁵ | 1.2 × 10⁻⁹ | 7.476 × 10⁻⁵ | 7.476 × 10⁻⁵ |
| 1600 | 0.0125 | 1.536 × 10⁻² | 1.869 × 10⁻⁵ | 7.6 × 10⁻¹¹ | 1.869 × 10⁻⁵ | 1.869 × 10⁻⁵ |

(a) Fitted orders: misalignment 2.00; raw orthogonality error 1.05 → 1.01 (tending to 1); renormalized orthogonality error 4.0 (it is the δ⁴ residual of one orthonormalization step). (b) Largest relative deviation from the prediction: 0.05 % at Δt = 0.1 s. (c) With exact steps: orthogonality error 7 × 10⁻¹⁵ (N = 200) to 9 × 10⁻¹⁴ (N = 1600), pure rounding; misalignment 0.0. Rotations about one fixed axis commute and their angles add, so N exact steps of ωΔt compose to exactly ωT; the only error is floating-point. (d) 0.1° = 1.745 × 10⁻³ rad; δ_max = √(3 × 1.745 × 10⁻³ / (0.2618 × 20)) = 0.0316 rad; Δt = δ/ω = 0.121 s, f = 8.3 Hz. For a slow roll the first-order step is fine at any realistic rate.

**Check.** The misalignment at N = 200 for ω = 30°/s (lesson table) is 9.554 × 10⁻³, eight times the 1.196 × 10⁻³ here: the formula has ω³. ✓

**Common mistakes.** Reading the misalignment from `euler_321_from_dcm` or `axis_angle_from_dcm` of the total matrix, which cannot see rotations beyond 180°; renormalizing and then reporting a "smaller error" because the orthogonality error went down.

## 5. Rodrigues

**Derivation.** (a) **v** = **n**(**n**·**v**) + **v**⊥, **v**⊥ = **v** − **n**(**n**·**v**). The vectors **v**⊥ and **n** × **v** are perpendicular, both perpendicular to **n**, and have equal length |**v**⊥| (because |**n** × **v**| = |**v**| sin(angle) = |**v**⊥|). Turning the *frame* by *a* about **n** (passive) changes the components of **v** in the plane as a clockwise turn of the vector would: the perpendicular part reads cos *a* **v**⊥ − sin *a* (**n** × **v**). Hence v_new = **n**(**n**·**v**) + cos *a* (**v** − **n**(**n**·**v**)) − sin *a* **n** × **v** = [cos *a* I + (1 − cos *a*) **nn**ᵀ − sin *a* [**n**×]] **v**. Sign check with **n** = e₃: first row (cos *a*, sin *a*, 0) = R₃(a). ✓ (The active matrix has + sin *a* [**n**×].)

(b) tr I = 3, tr **nn**ᵀ = **n**·**n** = 1, tr[**n**×] = 0: tr C = 3 cos *a* + 1 − cos *a* = 1 + 2 cos *a*. The symmetric terms cancel in C − Cᵀ, and [**n**×]ᵀ = −[**n**×], so C − Cᵀ = −2 sin *a* [**n**×].

(c) arccos has infinite slope at 1: arccos(1 − ε) ≈ √(2ε). A trace error of 2 × 10⁻¹⁶ gives ε = 10⁻¹⁶ in cos *a* and an angle error of √(2 × 10⁻¹⁶) = 1.4 × 10⁻⁸ rad, eight orders of magnitude worse than the input. atan2(2 sin *a*, 2 cos *a*) has slope 1 everywhere, so the angle error stays at 10⁻¹⁶. The code's exact-step run shows a misalignment of 0.0 with atan2 and about 10⁻⁷ with arccos.

(d) cos *a* = −1, sin *a* = 0: C = −I + 2**nn**ᵀ, symmetric; replacing **n** by −**n** leaves **nn**ᵀ unchanged. The axis sign is undetermined because a half-turn either way is the same rotation.

**Check.** For *a* = 180° about e₃, 2e₃e₃ᵀ − I = diag(−1, −1, 1) = R₃(180°). ✓

**Common mistakes.** Mixing active and passive signs on the skew term (the result is Cᵀ, which passes the orthogonality test and fails the SciPy comparison); normalising **w** when |**w**| = 0.

## 6. The wrong model

(a) δ = 2π × 10 / 100 = 0.628 rad = 36° per sample. Orthogonality error of one step: δ² = 0.395 (40 % of unity). Actual turn per step: atan(0.628) = 32.1° instead of 36°, a 3.9° lag per sample, 39°/s of drift. This is not a numerical-accuracy problem to be fixed by tolerances; the first-order model is being used outside the range (δ ≪ 1) where it is a model at all.

(b) 0.1° = 1.745 × 10⁻³ rad over T = 10 s at ω = 62.8 rad/s: δ_max = √(3 × 1.745 × 10⁻³ / 628) = 2.89 × 10⁻³ rad, f = ω/δ = 21.8 kHz. No flight gyro is sampled that fast for attitude; the sensible alternative is the exact Rodrigues step (or the equivalent quaternion update) per sample, which is exact for a fixed axis at any rate and leaves only the coning error from a changing axis, handled in Phase 05.

(c) Every function that works on the matrix (`rot_*`, `compose`, `transform`, `dcm_from_axis_angle`, `axis_angle_from_dcm`, `orthonormalize`, `integrate_constant_rate`) is valid at θ = 90°; `euler_321_from_dcm` is defined there only by convention, and any code that differentiates or filters the angles fails. Propagate the DCM itself (as this lesson's integrator does) or a unit quaternion, and convert to angles only for display.

(d) The non-rotating Earth: 7.29 × 10⁻⁵ rad/s × 7200 s = 0.525 rad = 30.1°. An attitude integrated in a frame treated as inertial for two hours is wrong by 30° unless Earth rate is subtracted from the gyros. (Wing flex, by comparison, is a few tenths of a degree and does not grow with time.)

**Common mistakes.** Answering (a) with "use smaller tolerances"; answering (b) by quoting the kHz number as the design.

## Check-yourself answers

1. Active: the vector turns, the frame stays; passive: the frame turns, the vector's components are re-expressed. This lesson uses passive, because flight computations re-express fixed physical vectors (gravity, velocity) in a moving body frame. The two differ by a transpose.
2. R₃(ψ) = ((cos ψ, sin ψ, 0), (−sin ψ, cos ψ, 0), (0, 0, 1)). Rows are unit length (cos² + sin² = 1) and mutually perpendicular (row 1 · row 2 = −cos sin + sin cos = 0; rows 1, 2 · row 3 = 0), so R₃R₃ᵀ = I, and for a square matrix this gives R₃ᵀR₃ = I too.
3. det(CᵀC) = (det C)² = 1 forces det C = ±1; a continuous turning of a rigid body from the identity (det +1) cannot jump to −1. A determinant of −1 is a reflection (mirror image), which swaps handedness.
4. C_b/n = R₁(φ) R₂(θ) R₃(ψ); the yaw R₃(ψ) is applied first because it is the right-most factor acting on v_n.
5. C₁₃ = −sin θ. At θ = ±90° the entries C₁₁, C₁₂, C₂₃, C₃₃ all vanish; only φ ∓ ψ is determined.
6. CᵀC − I = δ²(I − **nn**ᵀ), so the orthogonality error is δ² (= |δθ|²). Misalignment after N steps about a fixed axis: Nδ³/3 = ωT(ωΔt)²/3.
7. ω = 0.5236 rad/s, Δt = 0.01 s, T = 20 s: 0.5236 × 20 × (0.005236)²/3 = 9.57 × 10⁻⁵ rad = 0.0055°.
8. Orthonormalization moves a matrix to the *nearest* rotation, and the nearest rotation to a product of lagging steps is the lagging rotation itself; the lag is a property of the step, not of the orthogonality.
