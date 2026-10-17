# P01.16  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. tr A = 30 + 3s/2, det A = (20 + s)(10 + s/2) − 16. 2. λ = (tr ± √(tr² − 4 det))/2. 3. For each λ, **v** ∝ (−4, λ − (20 + s)); normalise; angle = atan2(v₂, v₁). 4. Dot the two unit vectors; add and multiply the eigenvalues. 5. `jacobi_eigen(A)` and `np.linalg.eigh(A)`; compare up to the sign of each eigenvector.

## 2.
(a) Take A**v** = λ**v**, multiply on the left by **v*** (conjugate transpose), show the scalar **v***A**v** equals its own conjugate. (b) Form **v**₂ᵀA**v**₁ two ways. (c) Expand the quadratic form term by term; use |**r**|² − (**r**·**n**)² = |**r**⊥|². (d) Write each of I_xx, I_yy, I_zz as its sum; the difference I_xx + I_yy − I_zz = 2 Σ m z² ≥ 0; equality iff all z = 0; invariance because the principal frame is just another frame and the inequality was proved for any frame.

## 3.
(a) Substitute **x** = M^{−1/2}**y** into M**ẍ** + K**x** = 0, multiply on the left by M^{−1/2}, assume **y** = **ŷ** e^{iωt}. Symmetric form: real eigenvalues guaranteed, orthonormal eigenvectors, Jacobi applies. (b) Build M^{−1/2} = diag(1/√m₁, 1/√m₂); S = M^{−1/2} K M^{−1/2}; eigenvalues are ω²; f = ω/2π. (c) **x** = M^{−1/2}**ŷ**; divide by x₁. (d) Evaluate K**x** − ω²M**x** and **x**₁ᵀM**x**₂.

## 4.
(a) `np.append(masses, 6.0)`, `np.vstack([positions, [0.4, -0.3, s/50]])`; `inertia_tensor`; `principal_axes`. (b) Old and new third rows of C_pb; angle = arccos|**n**_old · **n**_new|; tilt from z = arccos|n_z|. (c) `power_iteration(I_new)[2]` and the formula with the new ratio. (d) `rayleigh_quotient(I_new, [0, 0, 1])` against the largest moment; compare the gap with the old case.

## 5.
(a) **H** = I **e**_z = third column of I; angle = arccos(H_z/|**H**|). The wobble (nutation) follows because **H** is fixed in space while **ω** is not. (b) Minor-axis spin is stable only without dissipation; cite Explorer 1 (1958). (c) Two equal moments (`ValueError`): the axis is not isolated; the symmetric-top solution with nutation applies. (d) Set the tank mass to 4 kg, recompute, compare the ordering and the relative gap (I₃ − I₂)/I₃ in both cases.

## 6.
(a) Characteristic polynomial λ² + 2ζωλ + ω² = 0; discriminant 4ω²(ζ² − 1) < 0 for ζ < 1. (b) Call the function; call `power_iteration(A, max_iter=...)` inside `try/except RuntimeError`; note |λ₁| = |λ₂| = ω. (c) `expm(A*5) @ [1, 0]`; bound the norm using the eigen-decomposition (the constant comes from the non-orthogonality of the eigenvectors, bounded by the condition number of the eigenvector matrix, or simply report the ratio). (d) A ≠ Aᵀ; complex eigenvalues; Jacobi rotations preserve symmetry and real diagonals, which A does not have.
