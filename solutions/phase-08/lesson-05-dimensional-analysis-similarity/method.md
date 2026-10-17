# P08.05  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
(a) `isa_troposphere(0)`, `sutherland(T)`, `reynolds`, `mach`. (b) Same at 5000 m, with d and L. (c) Compare M with 0.3 and 1; Re with the transition range (about 10⁵–10⁶ for a flat plate).

## 2.
(a) Take logarithms: the dimension exponents of Π are D e; dimensionless iff D e = 0; the nullspace has dimension n − rank D (rank–nullity). (b) Two bases of a k-dimensional space are related by an invertible k × k matrix; in exponent space that matrix is products of powers. (c) `pi_groups_repeating` with ρ, V, c; F/(ρV²c²) = g(Re, M); drop M; require ∂F/∂ρ = 0 at fixed μ, V, c: g ∝ 1/Re. (d) The SVD vectors satisfy D B = 0 and are independent, so they span a subspace of the nullspace of the right dimension.

## 3.
(a) Copy `DIMENSIONS`, add `n`; `dimension_matrix`, `pi_count`, `pi_groups_repeating(D, names, ["rho", "n", "d"])`, `pi_groups_sympy`, `same_span`. (b) n = 2700/60 rev/s; J and πnd/a at 3000 m. (c) n_m = M_tip a₀/(π d_m), V_m = J n_m d_m; Re_tip = ρ(πnd)c/μ at both scales. (d) Compare which effect depends most on each number.

## 4.
(a) `full_scale_numbers`, `speed_for_reynolds` with sea-level ρ and μ, `mach`. (b) `pressurised_tunnel`, `cryogenic_tunnel` at 101 325 and 202 650 Pa; water V = Re ν/c. (c) q = ½ρV² in each case. (d) Weigh Mach error, loads, cost and Sutherland's validity at cryogenic temperatures.

## 5.
(a) Copy `DIMENSIONS`, add `h`; rank by `np.linalg.matrix_rank`. (b) `pi_groups_sympy`, `pi_groups_svd`, `pi_groups_repeating(D, names, ["rho", "V", "d", "k"])`, `is_dimensionless`, `same_span`; rewrite each as Re^a Nu^b Pr^c. (c) `np.linalg.svd(D)[1]`. (d) Physical meaning and integers versus orthonormality.

## 6.
For each test, identify the dominant physics and its number; check whether the tunnel or tank can match it; propose a correction (transition trips, extrapolation in Re, matched Froude with corrected friction, etc.).
