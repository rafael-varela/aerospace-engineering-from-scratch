# P02.14  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: your own fuselage seam)

The lesson's fuselage (r = 2.0 m, t = 1.6 mm, T = 0.80 MN m, σ_y = 325 MPa) with cabin pressure p = 0.050 + s/1000 MPa and a seam at a = 20 + s degrees to the longitudinal axis.

(a) Write the stress tensor in MPa.
(b) By hand, from the Mohr formulas (9), give the normal stress across the seam and the shear along it. Confirm with `seam_stresses`, with `transform_tensor2` and with `mohr_2d`.
(c) Give the three invariants, the in-plane principal stresses and the angle of σ₁ from the axis, by hand (Mohr's circle) and with `principal_stresses`.
(d) Give the true maximum shear, the von Mises stress and the margin of safety on yield.

Answer form: a 3 × 3 matrix, two stresses, three invariants and three principal numbers, three further numbers.

## 2. (Derivation: index gymnastics)

(a) Prove the ε–δ identity (2) by counting cases, and use it to show (a × b) · (c × d) = (a · c)(b · d) − (a · d)(b · c) and |a × b|² = |a|²|b|² − (a · b)².
(b) Derive the Rodrigues form (3) for a passive rotation of axes by a about the unit vector n: decompose a vector into its parts along n and perpendicular to n and write where the new basis vectors lie.
(c) Starting from Cauchy's formula t_i = σ_ij n_j and the vector law (4), derive the tensor law (5). State where the orthogonality Q_ik Q_jk = δ_ij is used.
(d) Show from (5) that σ_ii, σ_ij σ_ji and det σ are invariant, and write out (5) for a rotation about e₃ to obtain the three Mohr formulas (9).

Answer form: four derivations.

## 3. (Unfamiliar scenario: a strain-gauge rosette)

A 0°/45°/90° rosette bonded to a wing skin (2024-T3, E = 73 GPa, ν = 0.33) reads ε_a = 800 + 10s, ε_b = 300, ε_c = −200 microstrain along the 0°, 45° and 90° gauges. The gauges read normal strains; no gauge reads shear directly.

(a) Using the strain transformation law (the same as (9) with ε_x, ε_y, ε_xy in place of σ_x, σ_y, τ_xy), find ε_xy from the 45° reading.
(b) Give the in-plane principal strains and the angle of ε₁ from gauge a, by Mohr's circle and with `principal_stresses` applied to the strain tensor.
(c) The skin is in plane stress, not plane strain: σ₃₃ = 0, so ε₃₃ ≠ 0. Find ε₃₃ from the inverse law, build the full 3 × 3 strain tensor and compute the stress with `isotropic_stiffness` and `hooke_stress`. Check against the plane-stress formulas σ_x = E(ε_x + νε_y)/(1 − ν²), τ_xy = Eε_xy/(1 + ν).
(d) What stress tensor would you get if you set ε₃₃ = 0 by mistake? Which component is wrong, and by how much?

Answer form: one strain, two principal strains and an angle, a 3 × 3 stress tensor with a check, one wrong tensor with the error named.

## 4. (Numerical, uses your code: seam angle and torque)

For the fuselage of problem 1:

(a) Sweep the seam angle from 0° to 180° in 0.1° steps. Report the angles and values of maximum and minimum normal stress and of maximum shear, and relate each to the principal direction of problem 1(c).
(b) Reverse the torque. What does your seam see now? Which seam angle under the original torque sees the same stresses?
(c) Find, with `scipy.optimize.brentq` on `von_mises`, the torque T* at which the von Mises stress reaches σ_y/1.5 (a factor of safety 1.5 on yield), and check it by hand from σ_vm² = σ_x² − σ_xσ_y + σ_y² + 3τ_xy². What are the seam stresses at T*?
(d) Tabulate, for T = 0, 0.8, 1.6, 2.4, 3.2, 4.0 MN m, the shear, the seam stresses, the three principal stresses, θ_p, σ_vm and the margin. Where does the skin first go into compression, and what failure mode would you then expect that the yield margin cannot see?

Answer form: three angle–value pairs, two stresses and an angle, one torque with a hand check and two stresses, a six-row table with a sentence.

## 5. (Numerical, uses your code: a composite lamina as a fourth-order tensor)

A unidirectional carbon–epoxy lamina has E₁ = 140 GPa along the fibres, E₂ = 10 GPa across, G₁₂ = 5 GPa, ν₁₂ = 0.30, and is in plane stress. Its reduced stiffnesses are Q₁₁ = E₁/(1 − ν₁₂ν₂₁), Q₂₂ = E₂/(1 − ν₁₂ν₂₁), Q₁₂ = ν₁₂E₂/(1 − ν₁₂ν₂₁), Q₆₆ = G₁₂ with ν₂₁ = ν₁₂E₂/E₁.

(a) Build the 3 × 3 × 3 × 3 tensor C with C₁₁₁₁ = Q₁₁, C₂₂₂₂ = Q₂₂, C₁₁₂₂ = C₂₂₁₁ = Q₁₂, C₁₂₁₂ = C₁₂₂₁ = C₂₁₁₂ = C₂₁₂₁ = Q₆₆, all else zero. Check that it has the minor and major symmetries and that `hooke_stress` on a uniaxial strain ε₁₁ gives σ₁₁ = Q₁₁ε₁₁, σ₂₂ = Q₁₂ε₁₁.
(b) The fibres make angle θ = 30 + s degrees with the loading axis x. Rotating the *axes* by −θ puts the fibres at +θ in the new frame. Use `transform_tensor4` with `rotation_matrix([0, 0, 1], -θ)` and report C̄₁₁₁₁, C̄₂₂₂₂, C̄₁₁₂₂, C̄₁₂₁₂, C̄₁₁₁₂. Compare with the classical formulas Q̄₁₁ = Q₁₁c⁴ + 2(Q₁₂ + 2Q₆₆)c²s² + Q₂₂s⁴, Q̄₁₆ = (Q₁₁ − Q₁₂ − 2Q₆₆)c³s − (Q₂₂ − Q₁₂ − 2Q₆₆)cs³ (c = cos θ, s = sin θ).
(c) Form the 3 × 3 matrix Q̄ in engineering notation (rows σ_x, σ_y, τ_xy; columns ε_x, ε_y, γ_xy = 2ε_xy, so the third column is 2C̄_ij12), invert it, and give the effective E_x, E_y, G_xy and ν_xy at your angle, and the strains ε_x, ε_y, γ_xy for 100 MPa applied along x. Why is γ_xy non-zero?
(d) Plot E_x against θ from 0° to 90° in 15° steps. At what angle has E_x fallen to one-tenth of E₁?

Answer form: two checks, five stiffnesses with a comparison, four moduli and three strains with a sentence, a seven-point table with one angle.

## 6. (Judgment: is the seam acceptable?)

(a) For the fuselage of problem 1, give the margin on yield by von Mises and by Tresca. Which do you quote for a riveted seam, and why?
(b) Rank the lesson's assumptions (no stress concentration at the rivet holes, no secondary bending at the lap, no frame and stringer load sharing, thin wall) by how much each changes the seam's local normal stress, with a direction and a rough factor for the first two.
(c) The margin on yield is above 3. Explain in three sentences why the seam may still be the critical item on the aircraft, naming the quantity that actually sizes it.
(d) A colleague says "principal stresses depend on the axes you choose, so always compute them in the fuselage axes". Reply in three sentences, naming which quantities in this lesson are frame independent and which are not.

Answer form: two margins with a choice, a ranked list with two factors, two short paragraphs.
