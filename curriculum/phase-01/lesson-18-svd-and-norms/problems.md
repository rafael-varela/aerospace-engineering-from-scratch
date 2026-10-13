# P01.18  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every physical quantity (norms of a stiffness matrix are in N/m; condition numbers are dimensionless).

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: four norms and a condition number by hand)

Let A = [[s, 1], [0, 2]].

(a) Compute ‖A‖₁, ‖A‖∞ and ‖A‖_F from the entries.
(b) Form AᵀA, find its two eigenvalues from the trace and determinant, and hence σ₁, σ₂ and ‖A‖₂.
(c) Compute κ₂(A) = σ₁/σ₂, and compare with κ₁ = ‖A‖₁‖A⁻¹‖₁ and κ∞ (find A⁻¹ by hand).
(d) Verify (a)–(c) with your code and with `numpy.linalg.norm` / `cond`.

Answer form: four norms, two singular values, three condition numbers.

## 2. (Derivation: the inequalities behind the condition number)

(a) From the definition ‖A‖ = max ‖A**x**‖/‖**x**‖, prove ‖A**x**‖ ≤ ‖A‖‖**x**‖ and ‖AB‖ ≤ ‖A‖‖B‖.
(b) Prove ‖A‖₁ = max_j Σ_i |a_ij| (both the inequality and the vector that attains it).
(c) Using A = UΣVᵀ, prove ‖A‖₂ = σ₁ and ‖A‖_F² = Σσ_i², and prove ‖A‖₂ ≤ ‖A‖_F ≤ √n ‖A‖₂.
(d) Derive the bound ‖δ**u**‖/‖**u**‖ ≤ κ₂(K) ‖δ**f**‖/‖**f**‖ and show it is an equality when **f** ∝ **u**₁ and δ**f** ∝ **u**_n. What does the pairing **f** ∝ **u**_n, δ**f** ∝ **u**₁ give?
(e) Prove κ₂ ≥ 1 and κ₂(Q) = 1 for orthogonal Q, and κ₂(cA) = κ₂(A) for any scalar c ≠ 0.

Answer form: five short proofs.

## 3. (Unfamiliar scenario: a strain-gauge balance with a dead gauge)

A force balance has four strain gauges reading **g** = C**F**, with **F** = (F_x, F_y, F_z) in N and

$$
C = \begin{pmatrix} 1 & 0.2 & 0.1 \\ 0.1 & 1 & 0.2 \\ 0.2 & 0.1 & 1 \\ 0.5 & 0.5 & 0.5 + s/100 \end{pmatrix}.
$$

Simulate readings for **F** = (1000, −500, 300) N with Gaussian noise of 2 units on each gauge (`np.random.default_rng(1800 + s)`).

(a) Recover **F** from all four gauges with `pseudoinverse` and compare with `np.linalg.pinv`. Report κ₂(C), the residual and the standard errors of the three recovered loads (Lesson 17's formula with m = 4, n = 3).
(b) Gauge 4 dies. Recover **F** from gauges 1–3 (a square solve). How did κ₂ change, and why does losing a gauge make the recovery *worse* even though the system is now exactly determined?
(c) Gauges 3 and 4 die. Recover **F** from gauges 1–2 with the pseudoinverse. Show that the result solves the two equations exactly, that it is the shortest such solution, and that it is nevertheless far from the true load. Find the null-space direction and the component of the true load along it.
(d) The balance designer can add a fifth gauge with any row of unit norm. Which row would most reduce κ₂(C)? Give a reasoned choice and verify numerically.

Answer form: three load vectors with uncertainties where defined, two condition numbers, one null vector, one design choice.

## 4. (Numerical, uses your code: your own damaged element)

Damage element e = 1 + (s mod 10) of the 10-element bar to r = 10⁻¹, 10⁻², 10⁻³, 10⁻⁴ of its stiffness.

(a) Tabulate κ₂, σ₁ and σ_n against r; fit the exponent of κ₂ in r and explain the value.
(b) For r = 10⁻³ and a 1 kN tip load, apply ten random 1 % load perturbations and the worst-case pair; tabulate actual error against the bound.
(c) Compute the tip deflection for r = 10⁻³ and compare with F(9/k + 1/(rk)): which single element dominates, and by what factor?
(d) Precondition: form D = √diag(K) and compute κ₂(D⁻¹KD⁻¹) for each r. How much of the conditioning was stiffness contrast and how much was geometry?

Answer form: two tables, one exponent, two comparisons.

## 5. (Judgment: what a condition number is and is not)

(a) The intact bar has κ₂ = 175 and a 1 % load error produced a 1.36 % displacement error. A colleague says "the model is fine up to κ = 10⁴ then". Say what is right and what is wrong in that sentence.
(b) With 1000 elements, κ₂ ≈ 1.6 × 10⁶. How many digits can a double-precision solve lose, and does that matter for a displacement reported to three figures?
(c) The pseudoinverse needs a threshold for "zero" singular values. For the bar with r = 10⁻⁵, what does `rcond = 1e-4` do, and what physical statement does that amount to?
(d) Name a load case for a real aircraft structure where the worst-case alternating pattern is *not* unphysical, and say what that implies for the mesh.

Answer form: four short paragraphs, each with a number.

## 6. (Hard, uses your code: modes of a pressure field)

Build a 30 × 50 matrix of chordwise pressure-difference snapshots: for 30 angles of attack α from −2° to 10° (radians inside) and 50 chord stations x = (i + ½)/50, ΔC_p(α, x) = 4α √((1 − x)/x) + 2√(x(1 − x)), plus Gaussian noise of standard deviation 0.002 s/15 (`np.random.default_rng(600 + s)`).

(a) Compute the singular values of the (transposed, 50 × 30) matrix with `svd_jacobi`. How many are above the noise floor, and what is the floor's expected size?
(b) Give the relative Frobenius error of the rank-1, rank-2 and rank-3 approximations and the 2-norm error of each; confirm Eckart–Young.
(c) Interpret **u**₁, **u**₂ (the chordwise modes) and **v**₁, **v**₂ (their dependence on α) in terms of the two terms you built in.
(d) What fraction of the "energy" Σσ_i² is in the first two modes? Explain why a flow solver's output could be stored as two vectors and two coefficients per snapshot, and what would break that.

Answer form: a list of singular values, a 3×2 table of errors, two interpretations, one fraction.
