# P05.15  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe scripts that call it, given to the figures shown.

## 1. A drag-brace tube (s = 15: L = 0.70 m)

**Assumptions.** Straight, prismatic, axially loaded; ideal end conditions; typical 2024-T3 properties.

(a) d_i = 22.92 mm. **A = 94.12 mm², I = 6885 mm⁴, r = 8.553 mm.** **λ_c = √(2π² × 73.1 GPa / 345 MPa) = 64.7.**

(b) Pinned: KL/r = 700/8.553 = **81.8 > λ_c → Euler: P_cr = 10.14 kN** (107.7 MPa); the column curve and Euler coincide. Fixed-pinned (K = 0.6992): KL/r = **57.2 < λ_c → Johnson: σ = 210.0 MPa, P_cr = 19.76 kN**; Euler would give 20.74 kN.

(c) Only in the fixed-pinned case, where Euler overestimates by **4.9 %**. Close to λ_c the two curves are close (they are tangent there), so the error grows quickly as the slenderness drops further.

**Check.** The Euler stress at λ_c is F_cy/2 = 172.5 MPa; the pinned case's 107.7 MPa is below it, consistent with being on the Euler branch.

**Common mistakes.** Using the outside radius instead of r; forgetting to check the branch; applying K to the area.

## 2. Derivations

(a) w = A sin kx + B cos kx + Cx + D, k² = P/EI. Fixed base: w(0) = 0 → B + D = 0; w'(0) = 0 → Ak + C = 0. Free tip: w''(L) = −k²(A sin kL + B cos kL) = 0; and E I w''' + P w' = EI k³(−A cos kL + B sin kL) + P(Ak cos kL − Bk sin kL + C) = P C (the trigonometric terms cancel because P = EI k²), so **C = 0**, hence A = 0, and then B cos kL = 0. The smallest non-trivial root is kL = π/2: **P_cr = π²EI/(4L²)**, K = 2, and w = B cos kx − B ∝ **1 − cos(πx/2L)**.

(b) Fixed at 0: B + D = 0, Ak + C = 0. Pinned at L: w''(L) = 0 → A sin kL + B cos kL = 0; w(L) = 0 → A sin kL + B cos kL + CL + D = 0 → CL + D = 0 → −AkL − B = 0 → B = −AkL. Then A sin kL − AkL cos kL = 0, i.e. **tan kL = kL**. `brentq` on (4.4, 4.6): **β = 4.4934**, **K = π/β = 0.6992**. With A = 1: B = −β, C = −k, D = β:

**w(x) = sin kx − β cos kx − kx + β,  k = β/L.**

(c) Johnson σ_J = F_cy − Cλ², Euler σ_E = π²E/λ². Equal slope at λ_c: −2Cλ_c = −2π²E/λ_c³ → C = π²E/λ_c⁴. Equal value: F_cy − π²E/λ_c² = π²E/λ_c² → **λ_c² = 2π²E/F_cy**. Then C = π²E/λ_c⁴ = π²E F_cy²/(4π⁴E²) = **(F_cy/2π)²/E**, and **σ(λ_c) = π²E/λ_c² = F_cy/2**.

**Common mistakes.** Taking w' = 0 as the free-end condition (that is a fixed end); treating the free end's shear condition as E I w''' = 0, which forgets that the axial load tilts with the end; dividing by sin kL without noting when it is zero.

## 3. A launch-adapter truss member

**Assumptions.** Typical Ti-6Al-4V properties, F_cy = 880 MPa; K = 0.65 as credited; the 40 kN already includes the factor 1.5.

(a) **d_o = 30.0 mm** (30.04 mm), A = 91.2 mm², KL/r = 50.6, against λ_c = 50.5: the design sits **at the transition**, where Johnson and Euler agree (both give F_cy/2 = 440 MPa), so the branch does not matter at this point; a slightly larger tube would move onto Johnson.

(b) **d/t = 30** < 50: local wall buckling is not expected to govern.

(c) With K = 1: **d_o = 39.7 mm**, mass **0.539 kg/m** against **0.404 kg/m** for K = 0.65: the end-fixity credit saves **25 %** of the mass. That is why the credit must be justified by test or analysis of the joint, not assumed.

**Check.** At the transition the critical stress is F_cy/2 = 440 MPa and 440 MPa × 91.2 mm² = 40.1 kN. ✓

**Common mistakes.** Applying the 1.5 factor twice; ignoring the local-buckling check for a thin, large-diameter tube, which is how a minimum-mass optimiser "wins" with an unbuildable tube.

## 4. The finite-difference eigenproblem in detail

(a) Richardson errors (relative, against (4)): **pinned −1.1 × 10⁻⁸, fixed-free −1.9 × 10⁻⁹, fixed-pinned −6.4 × 10⁻⁸, fixed-fixed −1.7 × 10⁻⁷**, against P₁₀₀ errors of 8.2 × 10⁻⁵, 2.1 × 10⁻⁵, 2.0 × 10⁻⁴ and 3.3 × 10⁻⁴: three to four orders of magnitude better, which confirms that the error of the scheme is a clean h² term.

(b) The FD eigenvector at n = 200 and the exact shape agree to **7.6 × 10⁻⁶** after normalisation. The maximum deflection is at **x = 0.6L** from the fixed base (the base's zero slope pushes it past mid-length).

(c) B multiplies λ = P/EI, so a row of B is zero wherever an equation does not contain P: the four boundary rows of a pinned or fixed end, but only three for a column with a free end, whose shear condition contains P w'. B also has two zero columns, the outermost ghost nodes w₋₂ and w_{n+2}, which only the fourth difference in A reaches. B is therefore singular, det(A − λB) is a polynomial of degree lower than n + 5, and the missing roots are reported as infinite. At n = 20 (25 unknowns), `eigvals` gives **6 infinite eigenvalues for pinned-pinned, fixed-pinned and fixed-fixed (19 = n − 1 finite) and 5 for fixed-free (20 = n finite)**: B has rank 21 for the first three cases and 22 for fixed-free. The code keeps whatever is finite, real and positive and takes the smallest, which is robust to the exact count.

**Common mistakes.** Extrapolating with the wrong order (the factor 4 is for h²); comparing mode shapes without fixing the sign of the eigenvector.

## 5. An imperfect strut

**Assumptions.** Pinned ends; initial shape w₀ sin(πx/L); elastic until the extreme fibre reaches F_cy (first yield, a conservative definition of failure); c = 11.1 mm.

(a) With w = a sin(πx/L): E I (a − w₀)(−π²/L²) = −P a, so (a − w₀)P_E = P a and **a = w₀/(1 − P/P_E)**. The imperfection is amplified by 1/(1 − P/P_E): at the 15 kN ultimate load, by 2.21.

(b) σ(P) = (P/A)(1 + a(P)c/r²) = 435 MPa at **P = 19.40 kN**, against the Johnson load of the straight strut, **19.77 kN**: L/1000 crookedness costs 2 %.

(c) w₀ = 0.25 mm: **21.26 kN**; 1.0 mm: **17.06 kN**; 2.0 mm: **14.27 kN**, below the 15 kN ultimate load. The strut is moderately imperfection-sensitive: a 4 mm-per-metre bow (bent in handling, or welded crooked) takes the margin negative. This is why tube straightness is specified on drawings and why column curves from tests include an allowance.

**Check.** As w₀ → 0, the first-yield load tends to min(P_E, A F_cy) = 25.9 kN, above Johnson, because first yield of a straight column ignores the knee of the stress–strain curve that Johnson represents.

**Common mistakes.** Using P/A alone for the stress; letting P exceed P_E in the root search.

## 6. Judgment: which number goes in the report

A model answer, in outline. **Branch.** For 4130, λ_c = 96.4; at L/r = 60 the strut is in the Johnson range. Johnson gives 350.8 MPa and Euler 562.0 MPa: Euler overstates the strength by 60 %. No choice of K fixes that, because the error is in the formula, not the length. **"Welded so stiffer."** A stiffer end lowers K and raises both curves; it does not move a strut from the Johnson to the Euler branch at the same KL/r. Taking K = 1 is conservative only when combined with the right column curve. **What welding does.** It bows tubes as the welds cool, offsets the load line where several tubes meet, and creates a heat-affected zone of lower yield next to the joint; all three lower the real strength below the ideal-column value. **Recommendation.** Recompute with the column curve (Johnson here) and K = 1, add an allowance for straightness (or use test-based column allowables for welded 4130 tubes), and report the margin; ask for the tube straightness tolerance on the drawing, and for a compression test of a representative welded strut if the margin is small.

## Check-yourself answers

1. **P_cr = π²EI/L².** Buckling is a loss of stability of the straight shape, governed by bending stiffness and length; the material's strength enters only if the stress at P_cr exceeds the proportional limit, which is the Johnson range.
2. Fixed-free **K = 2**, load **¼** of pinned; fixed-fixed **K = 0.5**, load **4 ×** pinned.
3. KL/r, the effective length over the radius of gyration √(I/A): **600/7.5 = 80.**
4. At **λ_c = √(2π²E/F_cy)**, where the critical stress is **F_cy/2**.
5. **Johnson.** Euler would overestimate the strength, increasingly so at lower slenderness (unconservative).
6. The local slope dσ/dε of the stress–strain curve. A column that begins to bend while the load is still rising has every fibre still loading on the curve, so its resistance to bending is E_t I, not E I; in the knee E_t is well below E.
7. MS = 24/(1.5 × 12) − 1 = **+0.33.**
8. The boundary-condition rows of B are zero (they contain no P), so B is singular and the generalised problem has eigenvalues at infinity, one per null direction. Discard them (and any non-real or negative values) and take the smallest finite positive eigenvalue.
