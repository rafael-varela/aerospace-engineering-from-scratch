# P02.14  Tensor Notation & Index Gymnastics

> A stress is not nine numbers; it is one object that nine numbers describe in one frame and nine others describe in another, and the transformation law is what makes it one object.

| | |
|---|---|
| **Learning outcomes** | 1. Read and write index notation fluently: the summation convention, free and dummy indices, the Kronecker delta and Levi-Civita symbol, the ε–δ identity, and the dictionary between index expressions and matrix operations. 2. State and derive the transformation laws of vectors, second-order tensors and the fourth-order stiffness tensor under a rotation of axes, and show why the stress invariants, the principal stresses and the von Mises stress are frame independent. 3. Implement and test, from scratch, the two symbols as arrays, every contraction as an `einsum` string, the transformation laws, Cauchy's traction formula, invariants, principal stresses, the deviator, the 2-D Mohr formulas as a special case, and isotropic Hooke's law. 4. Decide the stresses opening and sliding a riveted seam in a pressurised fuselage skin under torque, the margin against yield, and which assumption of the thin-wall model fails first. |
| **Prerequisites** | `P01.15` (direction-cosine matrices, rotations of axes), `P01.16` (eigenvalues of symmetric matrices, orthogonal diagonalisation), `P01.13` (matrices and determinants), `P01.07` (dot and cross products), `P00.03`, `P00.07`. Skill check: you can write a 3 × 3 matrix product as a double sum, diagonalise a symmetric matrix, and explain why the rows of a direction-cosine matrix are orthonormal. |
| **Engineering question** | A transport-aircraft fuselage of radius 2.0 m with 1.6 mm skin carries a cabin pressure of 0.060 MPa and, during a rudder kick, a torque of 0.80 MN m: what are the normal and shear stresses across a riveted lap seam inclined 30° to the fuselage axis, what are the principal stresses and their directions, what is the margin against yield of the 2024-T3 sheet, and how do we know these answers do not depend on the axes we chose? |
| **Assumptions and limits** | Thin-walled membrane theory (t ≪ r): the stress is uniform through the thickness and the radial stress is zero (plane stress). Linear elastic, isotropic, homogeneous sheet; small strains; no bending at the lap joint, no stress concentration at the rivet holes, no frame or stringer load sharing, no buckling of the skin in shear. Cartesian axes only: the transformation laws are for rotations of orthonormal frames; general curvilinear tensor calculus (covariant and contravariant components, Christoffel symbols) is mentioned at the end as the next step and is original instruction, not part of the cited sources. Plausible numbers, not a specific aircraft. |
| **Analytical work** | The summation convention and its rules; δ_ij and ε_ijk; the ε–δ identity and the vector identities it generates; the dictionary to matrix notation; the transformation law of a vector from the geometry of direction cosines, of a second-order tensor from Cauchy's formula, of the stiffness tensor from Hooke's law; invariance of σ_ii, of the characteristic polynomial and of s_ij s_ij; principal stresses as the eigenproblem; the Mohr formulas as the i, j = 1, 2 block of Q_ip Q_jq σ_pq; isotropic C_ijkl from δ's; the fuselage stress state from equilibrium. |
| **Computational deliverable** | `starter/main.py`: `kronecker_delta`, `levi_civita`, `epsilon_delta_identity_residual`, `dot`, `cross`, `determinant`, `trace`, `double_dot`, `matrix_product`, `symmetric_part`, `antisymmetric_part`, `axial_vector`, `rotation_matrix`, `is_orthogonal`, `transform_vector`, `transform_tensor2`, `transform_tensor4`, `traction`, `normal_and_shear_on_plane`, `stress_invariants`, `principal_stresses`, `hydrostatic_and_deviatoric`, `von_mises`, `max_shear`, `mohr_2d`, `strain_from_displacement_gradient`, `isotropic_stiffness`, `hooke_stress`, `isotropic_strain`, `fuselage_skin_stress`, `seam_stresses`, `yield_margin`. Tests prove: ε_ijk is totally antisymmetric and satisfies the ε–δ identity exactly; the index-notation cross product, determinant, trace and double contraction agree with NumPy to 10⁻¹²; the Rodrigues matrix is a proper rotation, rotates axes passively, fixes its axis and returns to the identity at 2π; the transformation laws equal their matrix forms and compose; Cauchy's traction transforms as a vector; the invariants are unchanged under random rotations and are the coefficients of the characteristic polynomial; principal stresses match `eigh` and their planes carry no shear; the deviator is traceless and von Mises is rotation-invariant, pressure-blind, equals σ for uniaxial and √3 τ for pure shear; the Mohr formulas are the tensor law restricted to the plane; a rigid rotation gives zero strain; the isotropic stiffness has all its symmetries, is invariant under every rotation while a stiffened C₁₁₁₁ is not; Hooke's law round-trips with its inverse and is frame independent; the fuselage stress scales correctly with p, t and T and the seam stresses match the rotated tensor and lie between the principal values. |
| **Practical activity** | **computational experiment**: a sweep of seam angle, a sweep of torque, a random-rotation invariance table, and an anisotropic-stiffness rotation; optionally a stick of blackboard chalk twisted to failure to see the 45° principal plane. |
| **Equipment** | None required. A computer with the environment from `P00.01`. Optional: one stick of chalk. |
| **Safety** | None beyond normal computer use. If you twist chalk, do it over a bin; the pieces are blunt and the dust is a mild irritant. |
| **Evidence** | Completed `main.py`, passing tests, `stress_rotation.png`, the seam-angle and torque tables, the invariance table, and a short note stating the seam stresses, the margin, and the assumption you trust least. |
| **Assessment** | `assessments/phase-02/grading-rubric.md`; competencies 1–4. Pass: the derivations of the transformation law of σ_ij and of the Mohr formulas from it are complete; all public tests pass; the hand values of the seam stresses, the principal stresses and the margin agree with the code to three figures; the note names the frame-independent quantities used in the margin. |
| **References** | `@mit-ocw-3-11-mechanics-of-materials` Appendix "Matrix and Index Notation" (required-free; summation convention, δ_ij, matrix–index correspondence); Module "Tensor Transformations" (required-free; rotation of stress and strain, principal values, Mohr's circle); `@mit-ocw-18-06-linear-algebra` Lecture 21: Eigenvalues and eigenvectors and Lecture 25: Symmetric matrices and positive definiteness (required-free; the eigenproblem behind principal stresses); `@mit-ocw-8-07-electromagnetism-ii` "Electromagnetism II, Lecture Notes 2" (supplementary-free; rotations, the Levi-Civita symbol and index gymnastics in a physics setting). The fourth-order stiffness tensor, von Mises stress and the fuselage application are original instruction in this lesson. |
| **Extensions** | ★ Build the reduced stiffness of an orthotropic lamina (E₁ = 140 GPa, E₂ = 10 GPa, G₁₂ = 5 GPa, ν₁₂ = 0.3) as a fourth-order tensor and rotate it with `transform_tensor4`; recover the classical Q̄ formulas and plot the effective modulus against fibre angle (problem 5). ★ Write the strain-compatibility equations ε_ij,kl + ε_kl,ij − ε_ik,jl − ε_jl,ik = 0 in index form and count the independent ones. ★ Read about covariant and contravariant components: in a non-orthogonal basis the dot product becomes g_ij a^i b^j and the Kronecker delta becomes the metric; that is where this lesson's "index gymnastics" leads for curvilinear coordinates. |
| **Time** | unmeasured |

## 1. The engineering problem

A pressurised fuselage is a thin tube. The cabin pressure stretches it around (hoop) and along (longitudinal); a rudder kick twists it, adding a shear. The skin is riveted together along seams, and a seam does not run along a convenient axis: a lap seam on a lofted fuselage may lie at any angle to the longitudinal axis. The rivets care about the stress *across the seam* (which pulls the two sheets apart) and *along the seam* (which makes them slide), and neither of those is one of the hoop, longitudinal or shear stresses we compute from equilibrium. We need to turn the stress we know in one set of axes into the stress we want in another.

The mathematics that does this is the transformation law of a second-order tensor, and the notation that makes it painless is **index notation**: every component equation is written once, with subscripts, and a repeated subscript means a sum (`@mit-ocw-3-11-mechanics-of-materials` Appendix "Matrix and Index Notation"). The payoff is not just the seam stress. The same notation gives, in one line each, the principal stresses (the eigenproblem), the three invariants (the coefficients of the characteristic polynomial), the von Mises stress (a contraction of the deviator with itself), Hooke's law for any material (σ_ij = C_ijkl ε_kl), and the proof that the margin against yield does not depend on which axes we drew. In Phases 05 and 06 the stress and strain tensors, the stiffness tensor of a composite lamina, and the inertia tensor of a vehicle are all handled with exactly these rules.

Success criteria, written down first: the two symbols, built as arrays, satisfy the ε–δ identity exactly; every index contraction, written as an `einsum` string, agrees with NumPy's matrix routines to 10⁻¹²; the stress on the 30° seam computed three ways (Cauchy's formula, the full tensor transformation, the 2-D Mohr formulas) agrees to 10⁻¹⁰; the invariants and the von Mises stress are unchanged by random three-dimensional rotations to 10⁻⁹; and the hand calculation of §2.7 reproduces the code's seam stresses and margin to three figures.

## 2. Understand and derive

### 2.1 The summation convention

Write a vector as its components a_i, i = 1, 2, 3, in an orthonormal basis e_i. Einstein's convention: **a subscript that appears twice in a term is summed over 1, 2, 3**, and the summation sign is dropped. So a_i b_i means a₁b₁ + a₂b₂ + a₃b₃ = a · b, and A_ij x_j means Σ_j A_ij x_j, the i-th component of Ax. Rules: a repeated (**dummy**) index can be renamed freely (a_i b_i = a_k b_k); an index must not appear more than twice in a term; an index appearing once (**free**) must appear once in every term of the equation, and the number of free indices is the order of the object (0: scalar, 1: vector, 2: second-order tensor, 4: stiffness). The expression a_i b_i c_i is illegal; A_ij B_jk C_kl has free indices i, l and is the matrix product ABC.

Two symbols do most of the work (`@mit-ocw-8-07-electromagnetism-ii` Lecture Notes 2):

δ_ij = 1 if i = j, 0 otherwise (the identity matrix); ε_ijk = +1 for (1,2,3), (2,3,1), (3,1,2), −1 for (3,2,1), (2,1,3), (1,3,2), 0 if any two indices are equal.   (1)

δ_ij substitutes: δ_ij a_j = a_i, δ_ii = 3. ε_ijk is the cross product and the determinant: (a × b)_i = ε_ijk a_j b_k, det A = ε_ijk A_1i A_2j A_3k. Swapping any two indices of ε changes its sign; cycling them does not. The one identity to memorise is

ε_ijk ε_imn = δ_jm δ_kn − δ_jn δ_km.   (2)

Proof by counting: the left side is non-zero only if (j, k) and (m, n) are both the two indices other than i, in the same order (+1) or swapped (−1); the right side says exactly that. Everything about cross products follows from (2): a × (b × c) has components ε_ijk a_j ε_kmn b_m c_n = ε_kij ε_kmn a_j b_m c_n = (δ_im δ_jn − δ_in δ_jm) a_j b_m c_n = b_i (a_j c_j) − c_i (a_j b_j), the "BAC − CAB" rule, in one line. Problem 2 asks for two more.

The dictionary to matrices (`@mit-ocw-3-11-mechanics-of-materials` Appendix): A_ij B_jk ↔ AB; A_ji B_jk ↔ AᵀB; A_ii ↔ tr A; A_ij B_ij ↔ A : B = tr(AᵀB); a_i b_j ↔ the outer product abᵀ. Any second-order tensor splits uniquely into a symmetric part S_ij = (A_ij + A_ji)/2 and an antisymmetric part W_ij = (A_ij − A_ji)/2, and an antisymmetric tensor has only three independent components, which form the **axial vector** w_i = −½ ε_ijk W_jk with W_jk v_k = (w × v)_j. This is how the small-rotation part of a displacement gradient becomes a rotation vector (§2.6).

### 2.2 Rotating the axes: the direction-cosine matrix

Let a second orthonormal basis e′_i be related to the first by Q_ij = e′_i · e_j, the direction cosines (`P01.15`): row i of Q is the new basis vector e′_i written in old components. Since both bases are orthonormal, Q_ik Q_jk = δ_ij (QQᵀ = I) and, for a right-handed pair, det Q = +1. Rodrigues gives Q for a rotation of the axes by angle a about the unit vector n:

Q_ij = cos a δ_ij + (1 − cos a) n_i n_j + sin a ε_ijk n_k.   (3)

Check the three terms: for a = 0, Q = δ; the axis is fixed, Q_ij n_j = n_i, because ε_ijk n_j n_k = 0; the sign of the last term makes a rotation of the *axes* by +a about e₃ send the vector e₁ to components (0, −1, 0), as it must (if you turn your head left, the tree in front of you appears to the right). The derivation is in problem 2.

A vector has components a_i in the old basis and a′_i = a · e′_i in the new one. Expanding a = a_j e_j: a′_i = a_j e_j · e′_i, so

a′_i = Q_ij a_j   (vector transformation law).   (4)

A tensor is whatever transforms like a product of vectors. Here is why stress is one. Cauchy's formula says the traction (force per area) on a plane with unit normal n is t_i = σ_ij n_j (the derivation is the equilibrium of a small tetrahedron, `@mit-ocw-3-11-mechanics-of-materials` Module "Tensor Transformations"). Both t and n are vectors and obey (4) in any frame, and the formula must hold in the new frame too: t′_i = σ′_ij n′_j. Substitute t′_i = Q_ip t_p = Q_ip σ_pq n_q and n_q = Q_jq n′_j (the inverse of (4), using QᵀQ = I): t′_i = Q_ip Q_jq σ_pq n′_j for every n′, so

σ′_ij = Q_ip Q_jq σ_pq   (second-order tensor transformation law; σ′ = Q σ Qᵀ).   (5)

One Q per index. Strain, the inertia tensor, the moment of area tensor all transform the same way. The same argument one level up: Hooke's law σ_ij = C_ijkl ε_kl must hold in every frame, and σ and ε obey (5), so

C′_ijkl = Q_ia Q_jb Q_kc Q_ld C_abcd   (fourth-order law).   (6)

Four Q's. In code this is one `einsum` string, "ia,jb,kc,ld,abcd->ijkl", and reading the string *is* reading (6).

### 2.3 What does not change: invariants

Contract (5) on its own indices: σ′_ii = Q_ip Q_iq σ_pq = δ_pq σ_pq = σ_pp. The trace is the same in every frame: it is the first **invariant**, I₁. The same trick on σ′_ij σ′_ji gives σ_pq σ_qp, so I₂ = ½(σ_ii σ_jj − σ_ij σ_ji) is invariant, and det σ′ = det Q det σ det Qᵀ = det σ gives I₃. These three are the coefficients of the characteristic polynomial

det(σ − s δ) = −s³ + I₁ s² − I₂ s + I₃ = 0,   (7)

whose roots, the **principal stresses** σ₁ ≥ σ₂ ≥ σ₃, are therefore the same in every frame. They are the eigenvalues of the symmetric matrix σ (`@mit-ocw-18-06-linear-algebra` Lecture 25): there is an orthonormal set of **principal directions** n^(k) with σ_ij n_j^(k) = σ_k n_i^(k), and on each principal plane the traction is parallel to the normal, so the shear on it is zero. Put the principal directions as the rows of Q and (5) gives a diagonal σ′: that is orthogonal diagonalisation, read as a change of axes.

Two more frame-independent numbers come from the **deviator**. Split σ_ij = p δ_ij + s_ij with p = σ_kk/3 the mean (hydrostatic) stress; then s_ii = 0 and s transforms by (5) on its own because δ is invariant (Q_ip Q_jq δ_pq = δ_ij). The scalar s_ij s_ij is invariant, and the **von Mises stress**

σ_vm = √(3/2 s_ij s_ij) = √(½[(σ₁ − σ₂)² + (σ₂ − σ₃)² + (σ₃ − σ₁)²])   (8)

reduces to σ for uniaxial tension (so it can be compared directly with a tensile-test yield stress) and to √3 τ for pure shear; it ignores p entirely, which matches the observation that metals do not yield under hydrostatic pressure. The alternative is Tresca's maximum shear τ_max = (σ₁ − σ₃)/2, which acts on planes at 45° to the extreme principal directions; the two criteria differ by at most 2/√3 − 1 = 15 %.

### 2.4 Plane stress: Mohr's formulas are a special case

For a thin sheet, σ₃₃ = σ₁₃ = σ₂₃ = 0 and only σ₁₁ = σ_x, σ₂₂ = σ_y, σ₁₂ = τ_xy remain. Rotate the axes by θ about e₃: from (3), Q = [[cos θ, sin θ, 0], [−sin θ, cos θ, 0], [0, 0, 1]]. Writing out (5) for i, j ∈ {1, 2}:

σ_x′ = σ_x cos²θ + σ_y sin²θ + 2τ_xy sin θ cos θ = (σ_x + σ_y)/2 + (σ_x − σ_y)/2 cos 2θ + τ_xy sin 2θ,
σ_y′ = (σ_x + σ_y)/2 − (σ_x − σ_y)/2 cos 2θ − τ_xy sin 2θ,
τ_x′y′ = −(σ_x − σ_y)/2 sin 2θ + τ_xy cos 2θ.   (9)

These are the formulas of `@mit-ocw-3-11-mechanics-of-materials` Module "Tensor Transformations", and (σ_x′ − c)² + τ_x′y′² = R² with c = (σ_x + σ_y)/2, R = √[((σ_x − σ_y)/2)² + τ_xy²] is **Mohr's circle**: as θ runs through 180°, the point (σ_x′, τ_x′y′) runs once round a circle of centre c and radius R. The in-plane principal stresses are c ± R, at tan 2θ_p = 2τ_xy/(σ_x − σ_y); the in-plane maximum shear is R, at 45° from them. Note that the third principal stress is 0, so if both in-plane principal stresses are positive the *true* maximum shear is σ₁/2, larger than R; the plane that carries it is inclined through the thickness. Forgetting the third principal stress is the classic plane-stress error.

### 2.5 Stress on a seam

A seam whose line makes angle a with the longitudinal axis has, in the 1–2 plane, the unit normal n = (−sin a, cos a, 0). Cauchy's formula gives t_i = σ_ij n_j; the **normal stress** that opens the seam is σ_n = t_i n_i, and the **shear** that slides it is τ = |t − σ_n n|. By (5) these are also σ′₂₂ and |σ′₁₂| after rotating the axes by a, and also σ_y′ and τ_x′y′ of (9) with θ = a. Three routes, one answer, and the agreement is the first thing the code checks.

### 2.6 Strain, and Hooke's law as a fourth-order tensor

With displacement u_i(x), the gradient H_ij = ∂u_i/∂x_j splits into the small **strain** ε_ij = (H_ij + H_ji)/2 and the **rotation** ω_ij = (H_ij − H_ji)/2. A small rigid rotation u = w × x has H_ij = −ε_ijk w_k: purely antisymmetric, zero strain, and its axial vector is w. That is why the symmetric part alone measures deformation.

For a linear material σ_ij = C_ijkl ε_kl, with C the stiffness tensor: 81 components, reduced by σ_ij = σ_ji, ε_kl = ε_lk and the existence of a strain energy (C_ijkl = C_klij) to 21 independent ones in general. An **isotropic** material must have a C invariant under every Q in (6), and the only isotropic tensors available are products of δ's, so

C_ijkl = λ δ_ij δ_kl + μ (δ_ik δ_jl + δ_il δ_jk),   λ = Eν/((1 + ν)(1 − 2ν)),   μ = E/(2(1 + ν)),   (10)

two constants, which gives σ_ij = λ ε_kk δ_ij + 2μ ε_ij and, inverted, ε_ij = [(1 + ν) σ_ij − ν σ_kk δ_ij]/E. The test that (10) is invariant under random rotations, while the same C with C₁₁₁₁ increased by half is not, is the code's statement of what "isotropic" means. Problem 5 does the anisotropic case that matters in aerospace, a composite lamina.

### 2.7 The fuselage skin by hand

Equilibrium of the thin tube: cut it lengthwise and the pressure on the projected area 2rL is held by two skin strips of area tL, so the hoop stress is σ₂₂ = pr/t; cut it across and the pressure on πr² is held by the ring 2πrt, so the longitudinal stress is σ₁₁ = pr/(2t); a torque T on a closed thin tube is carried by a uniform shear flow q = T/(2πr²) (the moment of q about the axis over the enclosed area), so σ₁₂ = T/(2πr²t). With p = 0.060 MPa, r = 2.0 m, t = 1.6 mm, T = 0.80 MN m:

σ₁₁ = 0.060 × 2.0/(2 × 0.0016) = 37.5 MPa, σ₂₂ = 75.0 MPa, σ₁₂ = 0.80 × 10⁶/(2π × 4.0 × 0.0016) = 19.89 MPa.

Seam at a = 30°, by (9) with θ = 30°: c = 56.25 MPa, (σ_x − σ_y)/2 = −18.75 MPa, cos 60° = 0.5, sin 60° = 0.866: σ_n = σ_y′ = 56.25 + 18.75 × 0.5 − 19.89 × 0.866 = **48.40 MPa**; τ = |18.75 × 0.866 + 19.89 × 0.5| = **26.19 MPa**. Principal: R = √(18.75² + 19.89²) = 27.34 MPa, σ₁ = 56.25 + 27.34 = **83.59 MPa**, σ₂ = 28.91 MPa, σ₃ = 0; tan 2θ_p = 2 × 19.89/(37.5 − 75) = −1.061, 2θ_p = 133.3°, θ_p = **66.65°** from the longitudinal axis. Invariants: I₁ = 112.5 MPa, I₂ = 37.5 × 75 − 19.89² = 2416.7 MPa², I₃ = 0. True maximum shear (σ₁ − σ₃)/2 = 41.79 MPa, not R = 27.34. Von Mises from (8): √(½[(83.59 − 28.91)² + 28.91² + 83.59²]) = **73.53 MPa**; with σ_y = 325 MPa the margin of safety σ_y/σ_vm − 1 = **3.42**. Pressure alone gives σ_vm = √3/2 × 75 = 64.95 MPa and margin 4.00: the torque costs 15 % of the margin.

### 2.8 Where the model stops being valid

Membrane theory needs t/r ≪ 1 (here 0.0008) and no bending; a lap seam bends locally because the two sheets are offset by one thickness, and the rivet holes concentrate the seam's normal stress by a factor near 3 at the hole edge. Frames and stringers carry a share of the longitudinal and shear load. None of this changes the tensor algebra; it changes σ_ij. The algebra is also restricted to rotations of orthonormal Cartesian frames: in cylindrical or curvilinear coordinates the basis vectors change from point to point, components come in covariant and contravariant flavours, and derivatives acquire Christoffel terms. That is the subject of general tensor calculus, which the cited sources do not treat and this lesson does not need, since the skin is analysed in a local Cartesian frame.

## 3. Build the computational model

Open `starter/main.py`. Functions marked GIVEN are complete. Run `pytest curriculum/phase-02/lesson-14-tensor-notation` after each step; the tests are grouped in the same order.

### Step 1: the symbols and the dictionary

`kronecker_delta` is `np.eye(3)`. `levi_civita` fills a 3 × 3 × 3 array from the three even permutations and their swaps. `epsilon_delta_identity_residual` computes both sides of (2) with `einsum("ijk,imn->jkmn")` and returns the largest difference, which must be exactly 0. Then `dot`, `cross` ("ijk,j,k->i"), `determinant` ("ijk,i,j,k", with the three rows), `trace`, `double_dot`, `matrix_product`, `symmetric_part`, `antisymmetric_part`, `axial_vector` ("ijk,jk->i" with the −½). Write each `einsum` string by copying the index equation from §2.1; the string is the documentation.

### Step 2: rotations and the three transformation laws

`rotation_matrix` is (3) term by term; raise `ValueError` for a zero axis. `is_orthogonal` checks Q_ik Q_jk = δ_ij and det = +1. `transform_vector` ("ij,j->i"), `transform_tensor2` ("ip,jq,pq->ij"), `transform_tensor4` ("ia,jb,kc,ld,abcd->ijkl"). The test compares the second-order law with Q σ Qᵀ and checks composition.

### Step 3: traction, invariants, principal stresses, deviator, Mohr

`traction` normalises n, checks symmetry (raise `ValueError` otherwise) and returns σ_ij n_j. `normal_and_shear_on_plane` returns (t · n, |t − (t · n)n|). `stress_invariants` returns (σ_ii, ½(σ_ii σ_jj − σ_ij σ_ji), det). `principal_stresses` uses `numpy.linalg.eigh`, sorts descending and flips the third direction if the triad is left-handed. `hydrostatic_and_deviatoric`, `von_mises` by (8) in its s_ij s_ij form, `max_shear` = (σ₁ − σ₃)/2. `mohr_2d` is (9), which you should derive before you type it.

### Step 4: strain, Hooke's law, the fuselage

`strain_from_displacement_gradient` returns the two parts of H. `isotropic_stiffness` builds (10) from three `einsum` outer products of δ; raise `ValueError` unless E > 0 and −1 < ν < ½. `hooke_stress` is "ijkl,kl->ij"; `isotropic_strain` is the inverse law. `fuselage_skin_stress` is §2.7's equilibrium (raise `ValueError` unless r > 0 and 0 < t < 0.1 r); `seam_stresses` builds n = (−sin a, cos a, 0) and calls `normal_and_shear_on_plane`; `yield_margin` is σ_y/σ_vm − 1.

### Verify

`python main.py` prints the identity checks (ε–δ residual 0, cross product and determinant against NumPy 0, QQᵀ − I at 2 × 10⁻¹⁶), the stress tensor in MPa, the seam stresses three ways (48.40 / 26.19 MPa by Cauchy, by the rotated tensor, by (9)), invariants before and after rotation (112.500, 2416.714, 0.000 both times), principal stresses 83.59 / 28.91 / 0 MPa with σ₁ at 66.7° from the axis, maximum shear 41.79 MPa, mean 37.5 MPa, von Mises 73.53 MPa, margin 3.42, the strains for 2024-T3 (E = 73 GPa, ν = 0.33: ε₁₁ = 175, ε₂₂ = 858, ε₃₃ = −509, ε₁₂ = 362 microstrain) and the round trip C_ijkl ε_kl = σ to 10⁻¹⁴ MPa, the isotropic C invariant under a random rotation to 3 × 10⁻¹⁶ and the stiffened one changed by 42 %, and writes `stress_rotation.png`: seam normal and shear stress against seam angle, and Mohr's circle with the 1–2 axes, the seam axes and the principal points marked.

### Compare

Everything in §2.7 should match. If `seam_stresses` disagrees with `mohr_2d`, your normal vector or the sign of the ε term in (3) is wrong: a seam *along* the axis (a = 0) must be opened by the hoop stress, 75 MPa, and that single check fixes most sign errors.

## 4. Experiment at home, or investigate real data

The practical activity is a **computational experiment**; the chalk is optional.

**Experiment A: seam angle sweep.** Using `seam_stresses` at 1° steps, tabulate the angle of maximum normal stress, of minimum normal stress, and of maximum shear, and compare each with the principal direction 66.65° from §2.7. Expected: maximum normal stress 83.59 MPa at 156.65° (the seam *line* perpendicular to σ₁: a seam at 156.65° has its normal at 66.65°), minimum 28.91 MPa at 66.65°, and shear exactly zero at both, with maximum shear 27.34 MPa (= R) at 21.65° and 111.65°, 45° from the principal seams. The dashed line at 30° in the figure shows the lesson's seam sitting close to the worst shear angle, which is the kind of thing a seam layout review is for.

**Experiment B: torque sweep.** Run `fuselage_skin_stress` for T = 0, 0.4, 0.8, 1.2, 1.6, 2.0 MN m and tabulate σ₁₂, the 30° seam stresses, σ₁, θ_p, σ_vm and the margin. From the code: the seam normal stress falls from 65.6 to 22.6 MPa while the seam shear rises from 16.2 to 41.1 MPa; θ_p moves from 90° (pressure only: hoop is principal) towards 45° (pure torsion); σ_vm rises from 64.95 to 107.9 MPa and the margin falls from 4.00 to 2.01. Note the torque *relieves* the 30° seam's opening stress while loading its rivets in shear: a seam at the other hand (−30°) sees the opposite. By 3.0 MN m, σ₂ reaches zero and the skin would be in diagonal compression, where the real failure is shear buckling, which this model cannot see.

**Experiment C: invariance table.** Rotate σ by three random three-dimensional rotations (seed 14). The components scramble completely (one of the three gives σ′₁₁ = 83.0, σ′₂₂ = 5.1, σ′₃₃ = 24.4 MPa with shears up to 10.6 MPa), yet I₁ = 112.500 MPa, I₂ = 2416.714 MPa², I₃ = 0 and σ_vm = 73.526 MPa come back unchanged every time to 10⁻⁹. This is the experimental meaning of "invariant": a quantity you can compute in whatever frame is convenient and quote in any other.

**Experiment D: anisotropy.** Multiply C₁₁₁₁ of the isotropic stiffness by 1.5 and rotate it by a random Q: the change is 42 % of the largest component. Then rotate it by 90° about e₃ only: C₁₁₁₁ and C₂₂₂₂ swap, which is what a stiff axis does when you turn the material. Problem 5 does this properly for a carbon lamina and reproduces the textbook Q̄ formulas with one `einsum`.

**Optional chalk.** Twist a stick of blackboard chalk in two hands. It breaks on a helical surface at 45° to the axis: pure torsion is σ₁₂ = τ only, (9) gives θ_p = 45°, and a brittle material fails on the plane of maximum *normal* stress, σ₁ = τ. Bending it instead breaks it straight across (σ₁ along the axis). This is the whole of §2.4 in one snap.

## 5. Validate and improve

| quantity | hand (§2.7) | code | difference | cause |
|---|---|---|---|---|
| σ₁₁, σ₂₂, σ₁₂ | 37.5, 75.0, 19.89 MPa | 37.5, 75.0, 19.894 MPa | 0 | same equilibrium |
| seam normal, shear (30°) | 48.40, 26.19 MPa | 48.396, 26.185 MPa (three routes agree to 10⁻¹⁰) | rounding | — |
| σ₁, σ₂, θ_p | 83.59, 28.91 MPa, 66.65° | 83.588, 28.912 MPa, 66.65° | rounding | `eigh` vs Mohr's circle |
| I₁, I₂, I₃ | 112.5, 2416.7, 0 | 112.500, 2416.714, 0.000, unchanged under rotation | 0 | — |
| τ_max | 41.79 MPa | 41.794 MPa | 0 | (σ₁ − σ₃)/2 with σ₃ = 0, not R |
| σ_vm, margin | 73.53 MPa, 3.42 | 73.526 MPa, 3.42 | 0 | — |
| isotropic C under rotation | invariant | max change 2.8 × 10⁻¹⁶ | round-off | — |

Every line agrees, which proves the algebra and nothing about the aircraft. What the model leaves out, ranked by how much it would move the seam numbers: the stress concentration at the rivet holes (K_t ≈ 3 on the seam normal stress at the hole edge, so 48 MPa becomes about 145 MPa locally, which is why fatigue cracks start at rivet holes and why the margin against *yield* is not the design driver for a seam); secondary bending at the lap joint (a further factor of order 1.5 on the local stress on one face); load sharing by frames and stringers (reduces the skin shear and longitudinal stress by tens of percent); and the thin-wall approximation itself, which at t/r = 0.0008 is exact to that order. The invariance results are exact and transfer to any of these refinements unchanged; that is the point of writing the calculation in tensor form.

Improvements worth making: a Tresca margin alongside von Mises (here σ_y/(2τ_max) − 1 = 2.89, 15 % lower, the conservative choice for a shear-dominated seam); a fatigue allowable at the hole in place of σ_y; an orthotropic C for a composite skin (problem 5). Improvement worth *not* making: a finer sweep of seam angle. One-degree steps already locate the extremes to the figure shown, and the extremes are known analytically from θ_p.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (2), (5) and (9); `stress_rotation.png`; the tables of experiments A–C; and a one-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) giving the seam stresses, the principal stresses and direction, the margin by both criteria, and the assumption you trust least.

### Adversarial self-checks

- **What if the torque reverses?** σ₁₂ = −19.89 MPa. The 30° seam then sees 82.85 MPa normal and only 6.29 MPa shear (code), exactly what the −30° seam sees under the original torque: the opening stress nearly doubles. θ_p moves to 113.35° (the mirror of 66.65°); σ_vm, the invariants and the margin are unchanged because they are even in σ₁₂. A seam layout must be checked for both signs of a reversible load.
- **What if the skin thickness halves?** Every stress doubles (the tensor scales as 1/t), the principal directions do not move, σ_vm doubles to 147.1 MPa and the margin falls to 1.21. The thickness for a margin of exactly 1.0 at this load is 0.72 mm (code). What if p doubles? σ_vm rises to 134.4 MPa, less than double, because the shear does not scale with p.
- **Which assumption fails first in the real world?** The absence of stress concentration at the rivet holes: a factor near 3 on the seam normal stress at the hole edge, larger than every other correction in §5.
- **What measurement could falsify the model?** A 0°/45°/90° strain-gauge rosette on the skin away from the seam (problem 3): it gives ε₁₁, ε₂₂ and ε₁₂ directly, and the inverse Hooke's law predicts 175, 858 and 362 microstrain for the lesson's state. A measured ε₂₂/ε₁₁ far from 4.9 means the pressure stresses are not in the 2:1 membrane ratio (frames are taking load) and a non-zero ε₁₂ at zero rudder means a built-in twist.
- **Does the result survive a change of solver or discretisation?** The seam stresses by Cauchy's formula, by the full `einsum` transformation and by Mohr's closed form agree to 10⁻¹⁰; `eigh` and the Mohr-circle principal stresses agree to the same; three random three-dimensional rotations change no invariant beyond 10⁻⁹. There is no discretisation in this lesson to vary.
- **Which decision would change if the uncertainty doubled?** The decision is whether the 30° seam is acceptable at this torque. With a margin of 3.42 on yield, doubling any uncertainty in p, T or t leaves the yield decision untouched; the decision that moves is the fatigue life at the rivet holes, which depends on the seam normal stress range (48 MPa per pressurisation cycle here, 145 MPa at the hole edge), and that would need the fatigue curve for the riveted joint, not this lesson's yield check.

## Check yourself

1. Which of these are legal index expressions, and what is each: a_i b_i c_j, A_ij B_ji, A_ii B_jj, ε_ijk a_j a_k, δ_ij δ_ij?
2. Prove that ε_ijk a_j a_k = 0 and that δ_ii = 3. Use (2) to show (a × b) · (c × d) = (a · c)(b · d) − (a · d)(b · c).
3. Write the transformation law for a second-order tensor and say, in one sentence, why stress obeys it.
4. The trace of a tensor is invariant. Is the sum of the squares of its components invariant? Is the sum of its components?
5. A plane-stress state has σ_x = 80, σ_y = 20, τ_xy = 40 MPa. Give the centre and radius of Mohr's circle, both in-plane principal stresses and the true maximum shear.
6. Why does von Mises stress ignore the mean stress, and what experimental fact does that encode?
7. How many independent components has a general fourth-order stiffness tensor, and how many has an isotropic one? Which symmetries remove the rest?
8. A seam runs along the fuselage axis. Which stress opens it? Which stress would open a seam running around the circumference? Which seam angle in the lesson's load case carries no shear?
