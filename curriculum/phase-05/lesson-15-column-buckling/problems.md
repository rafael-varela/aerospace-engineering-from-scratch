# P05.15  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31). Material properties are the lesson's typical values; take F_cy equal to the typical F_ty unless told otherwise.

## 1. (Warm-up: a drag-brace tube)

A landing-gear drag brace is a 2024-T3 tube, 25.4 mm outside diameter with a 1.24 mm wall, of length L = (0.4 + 0.02s) m (E = 73.1 GPa, F_cy = 345 MPa).

(a) Compute A, I and r, and the transition slenderness λ_c.
(b) For pinned ends and for ideal fixed-pinned ends, compute KL/r, say which branch of the column curve applies, and give P_cr from the column curve and from Euler.
(c) In which case is Euler wrong, and by how much?

Answer form: three section properties, one slenderness limit, two pairs of loads, one percentage.

## 2. (Derivation: the other end conditions)

(a) From E I w'''' + P w'' = 0 and the boundary conditions of part 2.2, show that a column fixed at x = 0 and free at x = L buckles at P = π²EI/(4L²), with the shape 1 − cos(πx/2L).
(b) For a column fixed at x = 0 and pinned at x = L, show that the buckling condition is tan kL = kL. Solve it numerically for the smallest positive root β and give K = π/β to four figures. Write the mode shape.
(c) Derive the Johnson transition slenderness and constant (6) from the conditions of equal value and equal slope with the Euler curve, and show that the critical stress at the transition is F_cy/2.

Answer form: three derivations, one number to four figures, one mode shape.

## 3. (Unfamiliar scenario: a launch-adapter truss member)

A satellite's launch-adapter truss uses Ti-6Al-4V tubes (E = 113.8 GPa, F_cy = 880 MPa, density 4430 kg/m³) with a 1.0 mm wall, 0.80 m long. The joints are bonded titanium end fittings, credited in design with K = 0.65. Each member must carry 40 kN at ultimate load (the 1.5 factor is already included).

(a) Find, with `critical_stress` and `brentq`, the smallest outside diameter that gives P_cr = 40 kN. Which branch of the column curve is the design on?
(b) Check that the diameter-to-wall ratio stays below 50, a common limit for avoiding local wall buckling of round tubes in this range.
(c) Repeat (a) for pinned ends (K = 1) and compare the mass per metre. How much is the end-fixity credit worth?

Answer form: two diameters to 0.1 mm, one branch, two masses per metre and a percentage.

## 4. (Code: the finite-difference eigenproblem checked in detail)

(a) For the lesson's engine-mount tube (L = 500 mm) and each end condition, compute `fd_buckling` at n = 50 and 100 and form the Richardson extrapolation (4P₁₀₀ − P₅₀)/3. How close is it to the closed form, and how much better than P₁₀₀?
(b) For the fixed-pinned case, derive the exact mode shape w(x) = sin kx − β cos kx − kx + β (with k = β/L, the base fixed at x = 0) and compare it with the FD eigenvector at n = 200 after normalising both to a maximum of 1. Where along the strut is the maximum deflection?
(c) Explain, by looking at the rows and columns of B, why `scipy.linalg.eig` returns infinite eigenvalues. Count them for each end condition at n = 20 with `np.isfinite`, and say how many finite eigenvalues remain.

Answer form: a table of four extrapolated errors, one maximum difference, one position, four counts and a short argument.

## 5. (An imperfect strut)

Real struts are not straight. Give the lesson's pinned engine-mount strut an initial crookedness w₀ sin(πx/L) with w₀ = L/1000 = 0.5 mm.

(a) Show from E I (w − w_init)'' = −P w that the mid-length deflection under load is w₀/(1 − P/P_E), with P_E the Euler load.
(b) The extreme-fibre stress at mid-length is σ = (P/A)(1 + w c / r²) with c the outer radius. Find the load at which it reaches F_cy = 435 MPa, and compare with the Johnson load of the straight strut.
(c) Repeat for w₀ = 0.25, 1.0 and 2.0 mm, and comment on the imperfection sensitivity of this strut.

Answer form: one derivation, four loads, a comment.

## 6. (Judgment: which number goes in the stress report?)

A stress report for the engine mount sizes a compression tube with P_cr = π²EI/L² because "the joints are welded, which is stiffer than pinned, so Euler with K = 1 is conservative." The tube's L/r is 60 and the material is 4130.

Write a half-page review comment. Address: whether Euler applies at this slenderness at all; whether "welded so stiffer than pinned" justifies the formula used; what imperfections and load eccentricity at a welded cluster do; and what you would ask the author to change, with the evidence that would settle it.

Answer form: a structured argument of about 250 words.
