# P05.16  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
(a) U = ∫₀ᴸ M²/(2EI) dx with M = q(L − x)²/2 + Q(L − x). ∂U/∂Q at Q = 0: ∫ q(L − x)³/(2EI) dx = qL⁴/(8EI). (b) Same with C: ∫ q(L − x)²/(2EI) dx = qL³/(6EI). (c) m(x) = (L/2 − x) for x < L/2, zero beyond; ∫ M m/EI dx. (d) Define U(P) on a fine grid with `strain_energy_bending`; `castigliano_deflection(U, [0, 0], 0)` and `(..., 1)`.

## 2.
(a) U = ½ΣΣF_jk P_j P_k; ∂U/∂P_j = ½Σ_k(F_jk + F_kj)P_k = (FP)_j = δ_j once F is symmetric. (b) Order a, b: W = ½F_aa P_a² + ½F_bb P_b² + F_ab P_a P_b (P_a rides through the displacement caused by P_b). Order b, a: the cross term is F_ba P_b P_a. Equal energies ⇒ F_ab = F_ba. (c) M = q(L − x)²/2 − R(L − x); ∂U/∂R = −∫ M (L − x)/EI dx = 0 ⇒ R. For the mid-span deflection, use the real M of the propped beam and the unit-load moment of the released cantilever (m = L/2 − x for x < L/2). (d) Write M, N_spar, N_strut in terms of S; ∂U/∂S = 0; collect S.

## 3.
(a) M(φ) = PR cos φ; U = ∫₀^{π/2} M² R dφ/(2EI); δ_v = ∂U/∂P. (b) m(φ) = R(1 − sin φ) for a unit horizontal load; δ_h = ∫ M m R dφ/EI. (c) PR³/(3EI); compare: the curved path is longer (πR/2) and the lever arm stays large over much of it. (d) Ratio = (αPL/GA)/(PL³/3EI); A = bh, I = bh³/12; evaluate.

## 4.
(a) F = 6.0 kN (cos θ, sin θ), θ = −90 + 2s degrees; `truss_force_method`, `truss_stiffness_method`. (b) Unit load along F/|F|: n from `np.linalg.lstsq(C, e)`; δ = Σ N n L/EA; compare with u · e. (c) C over all eight degrees of freedom; reactions = (C N) at the six support degrees of freedom; sum with F; moments about the origin x R_y − y R_x. (d) Loop over the EA factor; observe the middle bar's force tend to a limit as its EA grows.

## 5.
(a) For each a, `braced_spar(5.0, a, h, EI, EA_spar, EA_strut, q)`; record max |M|, S, R, tip; pick the minimum. (b) Same with the uniform lift. (c) At a = 3.0 m compare S, R and the root moment. Explain with the location of the lift's centroid (4L/3π against L/2).

## 6.
Candidate differences: strut-fitting or bolt flexibility (lowers S); spar root modelled as a hinge or a flexible fitting (changes S, direction depends); spar modelled with shear deformation (raises S slightly); the strut's angle or attachment offset (changes sin θ); the lift applied differently (chordwise position, so torsion; a different spanwise distribution); geometric non-linearity; units or a wrong EA. Check first the load path totals (does the FE model's strut vertical component plus root shear equal the lift?) and the strut's axial stiffness EA/L in the model. Measurement: strain gauges on the strut during a ground proof test with sandbags.
