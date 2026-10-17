# P08.05  Dimensional Analysis & Similarity

> Physics does not care about units, so every law can be written in dimensionless groups, and the groups are the nullspace of a small matrix. Match the groups and a model behaves like the real thing; fail to match one and you must know what it costs.

| | |
|---|---|
| **Learning outcomes** | 1. Build the dimension matrix of a physical relation and show that its dimensionless groups are the nullspace, with n − rank groups (the Buckingham Pi theorem). 2. Find the groups three ways, exact rational nullspace by SymPy, orthonormal nullspace by the SVD, and the method of repeating variables, and prove that the three span the same space. 3. Implement and test a Pi tool that verifies every group is dimensionless and refuses repeating variables that are dimensionally dependent, and check that Reynolds and Mach numbers do not change with units. 4. Scale a wind-tunnel test of the light aeroplane's wing: show that Reynolds and Mach numbers cannot both be matched in the same air, and compute the pressurised, cryogenic and water tunnels that would. 5. Decide which groups a given model test must match and how to correct for those it cannot. |
| **Prerequisites** | `P01.14` (rank, nullspace, rank–nullity), `P01.18` (the SVD), `P00.05` (units and dimensions), `P00.09` (SymPy), `P03.02` (root finding), `P08.01` (ISA, Sutherland's law), `P08.04` (dynamic pressure, compressibility). Skill check: you can find the nullspace of a 3 × 6 matrix by row reduction and say what the rank of a matrix is. |
| **Engineering question** | Can a 1:5 model of the light aeroplane's wing be tested in a wind tunnel at the Reynolds and Mach numbers of the full-size wing in cruise (chord 1.5 m, 60 m/s, 3000 m), and if not, what would it take? The answers: full scale Re = 4.83 × 10⁶, M = 0.183. In a sea-level atmospheric tunnel the Reynolds number needs 235 m/s, which is M = 0.69: Mach is not matched. Both are matched in a tunnel pressurised to 3.84 bar at 288 K (62.1 m/s), or, as an indicative estimate, a cryogenic tunnel near 110 K at 1 atm; a water tunnel matches Re alone at 16.2 m/s. |
| **Assumptions and limits** | The relation is dimensionally homogeneous and the list of variables is complete (a missing variable gives a wrong set of groups that the algebra cannot detect); four base dimensions (M, L, T, Θ). The tunnel calculations use the ISA, Sutherland's law (an extrapolation at cryogenic temperatures, so the cryogenic figures are indicative only) and air as an ideal gas; geometric similarity of the model is assumed exact, including surface finish. Partial similarity is the rule in practice, not the exception. |
| **Analytical work** | The dimension matrix and the Pi theorem from rank–nullity; non-uniqueness of the basis; the method of repeating variables as a linear solve; the force on a wing as C(Re, M); the creeping-flow limit F ∝ μVc and Stokes's 3πμVd; why Re and M cannot both be matched (Re/M = ρac/μ); the pressurised and cryogenic tunnel conditions. |
| **Computational deliverable** | `starter/main.py`: `dimension_matrix`, `pi_count`, `pi_groups_sympy`, `pi_groups_svd`, `same_span`, `is_dimensionless`, `pi_groups_repeating`, `reynolds`, `mach`, `full_scale_numbers`, `speed_for_reynolds`, `pressurised_tunnel`, `cryogenic_tunnel`, `drag_coefficient` (`DIMENSIONS`, `format_group`, `isa_troposphere`, `sutherland`, `speed_of_sound`, `stokes_drag` and `main` are given). Tests prove: the matrix columns; the Pi count n − rank for five relations; SymPy groups integer and dimensionless; the SVD basis orthonormal with |DB| < 10⁻¹²; SymPy and SVD span the same space for three relations and for ten seeded random matrices, some rank deficient; `same_span` detects different spaces; repeating variables give the force coefficient, 1/Re and 1/M and refuse dependent choices; Re unchanged in cgs and imperial units and M in knots; full-scale numbers from the shared atmosphere; a Re-matched 1:k model runs at k times the Mach number; pressurised and cryogenic solutions match both numbers when substituted back, the pressure scaling with k; Stokes drag collapses to C_D Re = 24 in three fluids. |
| **Practical activity** | **computational experiment**: groups for four engineering relations by three routes; agreement of SymPy and SVD on 200 seeded random dimension matrices; a model-scale sweep of the Reynolds-matching speed and its Mach number; the pressure and temperature that match both, against scale; and the water-tunnel alternative. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `similarity.png`, the derivation of the Pi count and of F = ½ρV²c²C(Re, M), the experiment tables, and a half-page test plan for a model of the wing with what it gives up. |
| **Assessment** | `assessments/phase-08/grading-rubric.md`; competencies 1–4. Pass: the Pi theorem is derived from rank–nullity and the three routes are shown to agree; all public tests pass; the conflict between Re and M is shown algebraically and numerically; the test plan names the group it matches, the one it does not, and the correction. |
| **References** | `@mit-ocw-2-25-advanced-fluid-mechanics` class notes "The Buckingham Pi Theorem in Dimensional Analysis" (Historical Note; The Theorem) and "The Physical Basis of Dimensional Analysis", https://ocw.mit.edu/courses/2-25-advanced-fluid-mechanics-fall-2013/resources/mit2_25f13_the_buckingham/ (required-free: the theorem and its proof). `@leishman-intro-aerospace-flight-vehicles` chapters "Dimensional Analysis" (The Buckingham Π Method; Π Method Applied to a Wing; When Things Go Wrong in the Buckingham Π Method) and "Dynamic Similarity" (Similarity Criteria; Partial Similarity) (required-free; link only). `@mit-ocw-16-100-aerodynamics` Lecture 31 "Similarity in Wind Tunnel Testing", https://ocw.mit.edu/courses/16-100-aerodynamics-fall-2005/resources/16100lectre31_cj/ (required-free: why M and Re cannot both be matched). `@libretexts-bar-meir-fluid-mechanics` Ch. 9 "Dimensional Analysis": 9.2 Buckingham–π–Theorem (9.2.1 Construction of the Dimensionless Parameters; 9.2.4 Similarity and Similitude); 9.4 Summary of Dimensionless Numbers (supplementary-free). `@libretexts-marchman-aerodynamics-aircraft-performance` 1.13 "Similarity Parameters" (supplementary-free). `@nasa-glenn-beginners-guide-aeronautics` "Similarity parameters" (supplementary-free). `@sympy-tutorial` (supplementary-free: `Matrix.nullspace`). |
| **Extensions** | ★ Add a fifth base dimension (electric current) and find the groups of a magnetohydrodynamic flow (Hartmann number). ★ Write `best_repeating(names)` that searches all dimensionally independent triples and returns the set whose groups have the smallest integer exponents. ★ Fit a dimensionless drag law to the seeded synthetic data of `P08.16` in dimensional and in dimensionless form and compare the scatter. |
| **Time** | unmeasured |

## 1. The engineering problem

Wind tunnels test models because full-size aeroplanes do not fit in them and because a model costs less to change. The model's numbers are only worth something if the flow around it is **similar** to the full-size flow: the same streamline pattern, the same pressure coefficients, the same separation, so that a force coefficient measured on the model applies to the aircraft. Which numbers have to match for that? Dimensional analysis answers it without solving any flow: it reduces the dozen dimensional variables of a problem to a handful of dimensionless groups, and similarity means equal groups.

The light aeroplane's wing in cruise has a chord of 1.5 m at 60 m/s in the 3000 m air of `P08.01`. A tunnel model at 1:5 scale has a 0.3 m chord. To have the same Reynolds number, the ratio of inertial to viscous forces that controls the boundary layer, the tunnel must blow five times faster in the same air; to have the same Mach number, the ratio of flow speed to sound speed that controls compressibility, it must blow at the same speed. It cannot do both. The questions are by how much, and what tunnel would.

Success criteria, written first: the Pi tool finds n − rank groups for every relation, verifies each dimensionless, and its SymPy, SVD and repeating-variable bases span the same space; Re and M computed from SI, cgs and imperial inputs agree; and a tunnel condition that claims to match both numbers matches them when substituted back.

## 2. Understand and derive

### 2.1 The dimension matrix

Every variable has dimensions that are products of powers of base dimensions: mass M, length L, time T and temperature Θ. Force is M L T⁻², viscosity M L⁻¹ T⁻¹. Write the exponents as a column; the columns of all n variables form the **dimension matrix** D (4 × n). For the force on a wing, F = f(ρ, V, c, μ, a):

| | F | ρ | V | c | μ | a |
|---|---|---|---|---|---|---|
| M | 1 | 1 | 0 | 0 | 1 | 0 |
| L | 1 | −3 | 1 | 1 | −1 | 1 |
| T | −2 | 0 | −1 | 0 | −1 | −1 |
| Θ | 0 | 0 | 0 | 0 | 0 | 0 |

### 2.2 The Pi theorem is a statement about a nullspace

A product Π = F^{e₁} ρ^{e₂} V^{e₃} c^{e₄} μ^{e₅} a^{e₆} has the dimensions M^{(De)₁} L^{(De)₂} T^{(De)₃} Θ^{(De)₄}: taking logarithms turns products of powers into a matrix–vector product. So

**Π is dimensionless ⇔ De = 0,**   (1)

and the exponent vectors of all dimensionless groups form the **nullspace** of D (`P01.14`). By rank–nullity it has dimension

**k = n − rank D.**   (2)

That is the **Buckingham Pi theorem** (2.25, "The Buckingham Pi Theorem in Dimensional Analysis", The Theorem; Bar-Meir 9.2). A dimensionally homogeneous relation among the n variables can be rewritten as a relation among k independent groups. For the wing, rank D = 3 (M, L and T appear; Θ does not), so k = 3. Note that the count uses the rank, not the number of base dimensions: they differ whenever a base dimension is absent or two are tied together.

**The basis is not unique.** Any invertible combination of a basis is another basis, and in variables a combination of exponent vectors is a product of powers of groups. (Re, Nu, Pr) and (Re, Nu, St) with St = Nu/(Re Pr) describe the same heat-transfer physics. Which basis to report is a choice of convenience.

### 2.3 Three routes to the nullspace

1. **Exact (SymPy).** `sympy.Matrix(D).nullspace()` row-reduces D in rational arithmetic and returns one basis vector per free column; scaled to the smallest integers, these are exponent vectors you can read. For the wing it returns F/(ρV²c²), Fρ/μ² and V/a.
2. **Floating point (SVD).** D = UΣVᵀ; the right singular vectors whose singular values vanish span the nullspace (`P01.18`). They are orthonormal and irrational, so they are correct but unreadable; they are the robust choice when D comes from data, and the singular values themselves say how sharply rank-deficient D is.
3. **Repeating variables.** Choose as many dimensionally independent variables as the rank (for the wing ρ, V, c, which contain M, L and T between them). For each other variable q, find exponents x so that q·ρ^{x₁}V^{x₂}c^{x₃} is dimensionless: the linear system D_r x = −D_q. This is the textbook method (Leishman, "Dimensional Analysis", The Buckingham Π Method), and it gives one group per non-repeating variable:

   **Π_F = F/(ρV²c²),  Π_μ = μ/(ρVc) = 1/Re,  Π_a = a/V = 1/M.**   (3)

   If D_r is singular (for example V and a as two repeating variables: they have the same dimensions) the method fails, and the code refuses.

The three bases are different lists of groups and the same space: SymPy's Fρ/μ² is Π_F × Re². The code checks it with `same_span`, which compares rank(B₁), rank(B₂) and rank([B₁ B₂]).

### 2.4 What the groups say

From (3) the force on a wing is

**F = ½ρV²c² C(Re, M),**   (4)

with ½ a convention and c² standing for a reference area. Dimensional analysis cannot find C; it says that C depends on two numbers only, so a model test at the same Re and M measures the full-size C. Two limits show its power:

- **Creeping flow** (Re → 0): inertia is irrelevant, so F cannot depend on ρ. In (4) that forces C ∝ 1/Re, hence **F = K μVc**. For a sphere the exact solution gives K = 3π with c = d (Stokes; `P08.16`), so C_D = F/(½ρV²πd²/4) = 24/Re. The tests check that Stokes drag collapses to C_D·Re = 24 in water, air and glycerine-like oil: one dimensionless law, three fluids.
- **Incompressible flow** (M → 0): a drops out and C = C(Re).

### 2.5 Why Re and M cannot both be matched

Divide the two groups: **Re/M = ρac/μ**. For the same fluid at the same state, ρa/μ is fixed, so Re/M scales with the length c. A 1:k model in the same air at the Reynolds-matching speed runs at **k times the full-scale Mach number** (the tests check it exactly). For the light aeroplane, Re = 4.83 × 10⁶ and M = 0.183 at full scale; a 1:5 model in a sea-level tunnel needs 235.2 m/s, M = 0.69, well into the range where compressibility changes the pressure distribution (`P08.04`). To match both, ρa/μ must rise by k (16.100 Lecture 31; Leishman, "Dynamic Similarity", Partial Similarity):

- **Pressurise**: at fixed T, a and μ are fixed and ρ ∝ p. Mach fixes V = Ma(T); Reynolds then fixes ρ = Re μ/(Vc) and p = ρRT. For the 1:5 model at 288.15 K: V = 62.14 m/s, ρ = 4.637 kg/m³, **p = 3.835 bar**. Substituted back: Re = 4.8308 × 10⁶, M = 0.182605.
- **Cool**: lowering T lowers μ (by Sutherland) and a, and raises ρ at fixed p; Re at fixed M rises steadily as T falls, so a root finder on T works. At 1 atm the 1:5 model needs about **110 K** (V = 38.4 m/s), an indicative figure because Sutherland's law is extrapolated there; cryogenic tunnels use nitrogen and its measured properties.
- **Change fluid**: water's kinematic viscosity is 1/15 of air's, so a water tunnel matches Re at 16.16 m/s; Mach is irrelevant in water but cavitation and free-surface effects are not.

### 2.6 Where the method stops

The Pi theorem is only as good as the variable list. Leave out the speed of sound and you will never find M; include an irrelevant variable and you get a spurious group that experiments must show to be unimportant. Dimensional analysis also cannot tell which groups matter most, or in what range; that takes physics or experiments (Leishman, "Dimensional Analysis", When Things Go Wrong in the Buckingham Π Method). And geometric similarity includes surface roughness and transition, which a model rarely matches: at lower Reynolds numbers a model's boundary layer stays laminar longer, which is why tunnel models carry trip strips.

## 3. Build the computational model

Open `starter/main.py`. `DIMENSIONS` (exponents of 22 common variables), `format_group`, `isa_troposphere`, `sutherland`, `speed_of_sound`, `stokes_drag` and `main` are given.

### Step 1: the matrix and the count

`dimension_matrix(names, dims)` stacks the exponent columns; `pi_count` is (2) with `np.linalg.matrix_rank`.

### Step 2: three nullspaces

`pi_groups_sympy` converts D to a `sympy.Matrix` of integers, calls `nullspace()`, and scales each vector to integers with the least common multiple of the denominators and the greatest common divisor of the result, sign fixed so the first nonzero entry is positive. `pi_groups_svd` returns the columns of V with singular values below a relative tolerance:

```python
_, s, Vt = np.linalg.svd(D)
rank = int(np.sum(s > tol * s[0]))
return Vt[rank:].T
```

`same_span` and `is_dimensionless` are rank and residual checks. `pi_groups_repeating` solves D_r x = −D_q with `np.linalg.lstsq`, refuses a rank-deficient D_r, and verifies each result with `is_dimensionless`.

### Step 3: similarity

`reynolds`, `mach` and `drag_coefficient` are one line each; `full_scale_numbers` uses the ISA and Sutherland; `speed_for_reynolds` solves Re = ρVc/μ for V; `pressurised_tunnel` and `cryogenic_tunnel` are 2.5, the second with `scipy.optimize.brentq` on log(Re(T)/Re) between 20 K and 400 K, raising when no root is bracketed.

Running `main()` prints

```
wing force F = f(rho, V, c, mu, a): n = 6, rank = 3, 3 groups
  SymPy nullspace (integers):
    F^1 rho^-1 V^-2 c^-2                     D e = [0 0 0 0]
    F^1 rho^1 mu^-2                          D e = [0 0 0 0]
    V^1 a^-1                                 D e = [0 0 0 0]
  SVD nullspace: 3 orthonormal vectors, |D B| max 3.1e-16, same span as SymPy: True
  repeating variables rho, V, c:
    Pi_F: F^1 rho^-1 V^-2 c^-2
    Pi_mu: rho^-1 V^-1 c^-1 mu^1
    Pi_a: V^-1 a^1
  same span as SymPy: True

pipe pressure drop dp = f(rho, V, d, mu, eps, L): 4 groups:
    Pi_dp: dp^1 rho^-1 V^-2
    Pi_mu: rho^-1 V^-1 d^-1 mu^1
    Pi_eps: d^-1 eps^1
    Pi_L: d^-1 L^1

light aeroplane wing: chord 1.5 m, TAS 60.0 m/s at 3000 m: Re = 4.8308e+06, M = 0.1826
model tests in a sea-level atmospheric tunnel (rho 1.225, mu 1.7894e-05, a 340.29):
   scale   chord (m)   V for Re match (m/s)   tunnel Mach   Re at full-scale Mach
   1:1       1.500             47.0               0.138        6.381e+06
   1:2       0.750             94.1               0.276        3.191e+06
   1:5       0.300            235.2               0.691        1.276e+06
   1:10      0.150            470.4               1.382        6.381e+05
1:5 model, both matched, pressurised tunnel at 288.15 K: V = 62.14 m/s, rho = 4.637 kg/m^3, p = 3.835 bar (3.79 atm)
   substituted back: Re = 4.830847e+06, M = 0.182605
1:5 model, both matched, cryogenic tunnel at 1 atm: T = 110.2 K, V = 38.42 m/s, rho = 3.204 kg/m^3 (Sutherland's law extrapolated far below ordinary temperatures: indicative only)
1:5 model in a water tunnel (20 C): V = 16.16 m/s for the Re match (Mach irrelevant)
```

and writes `similarity.png`: the tunnel Mach number at the Reynolds-matching speed against scale, with the full-scale Mach number and M = 1, and the pressure that matches both against scale.

### Verify

Run `python -m pytest curriculum/phase-08/lesson-05-dimensional-analysis-similarity`. The 24 tests:

1. **Matrix**: shape and the force and viscosity columns; the Pi count n − rank for the wing (3), the pipe (4), a propeller power relation (2), a Froude relation (1) and a Prandtl relation (1).
2. **SymPy**: integer exponents, dimensionless, the right number of groups.
3. **SVD**: orthonormal to 10⁻¹² with |DB| < 10⁻¹².
4. **Agreement**: SymPy and SVD span the same space for three relations and ten seeded random integer matrices (half of them forced rank-deficient); `same_span` rejects two different lines and accepts a rescaled one.
5. **Repeating variables**: exactly (3) for ρ, V, c, spanning the SVD space; a `ValueError` for V, a, c.
6. **Helpers**: `is_dimensionless` on Re and on ρV; `format_group`.
7. **Units**: Re in SI, cgs and imperial (slug, ft, lbf·s/ft²) agree to 10⁻⁶; M in knots equals M in m/s.
8. **Similarity**: full-scale ρ, μ, Re from the shared atmosphere; a 1:k model at the Re-matching speed runs at exactly k times the Mach number; pressurised and cryogenic tunnels match both numbers when substituted back, the pressure doubling when the model halves; an impossible cryogenic case is refused; Stokes drag gives C_D·Re = 24 in three fluids.
9. `main` writes the figure.

### Compare

SymPy's exact rational row reduction and NumPy's floating-point SVD are independent algorithms; they agree on the span for every relation in the lesson and for 200 of 200 seeded random 4 × n integer matrices (0.12 s in all). The method of repeating variables, a third route, gives the same span with the familiar groups. The SVD has one advantage the exact route lacks: its singular values (5.64, 3.96, 0.61, 0.36 for the heat-transfer matrix of Problem 5) say how far D is from losing rank, which matters when the exponents come from a fitted or measured relation rather than a table.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction before running it.

1. **Four relations, three routes.** The wing force, the pipe pressure drop, the propeller thrust (Problem 3: F, ρ, n, d, V, μ, a) and the duct heat-transfer coefficient (Problem 5: h, ρ, V, d, μ, k, c_p). Predict the counts. *Finding:* 3, 4, 4, 3; for the pipe the groups are Δp/(ρV²), 1/Re, ε/d and L/d; for heat transfer, SymPy gives Re, Nu and St and repeating ρ, V, d, k gives Nu, 1/Re and Re Pr, all spanning the same space.
2. **Random matrices.** 200 seeded random integer 4 × n matrices (n = 4–7, entries −3 to 3, seed 5). Predict whether SymPy and SVD ever disagree. *Finding:* they agree on the count and the span in 200 of 200. A disagreement would mean a tolerance problem in the SVD's rank decision, the only judgment either algorithm makes.
3. **Scale sweep in an atmospheric tunnel.** The table above, and the figure from 1:1 to 1:10. Predict the scale at which the Re-matched model goes supersonic. *Finding:* M = 0.138, 0.276, 0.691, 1.382 at 1:1, 1:2, 1:5, 1:10; supersonic beyond about 1:7.3 (M ∝ k, 0.138k = 1). At the full-scale Mach number instead, the model's Re is 6.38 × 10⁶/k: 1.28 × 10⁶ at 1:5, a quarter of the flight value.
4. **What matches both.** Pressurised tunnel at 288.15 K and cryogenic tunnel at 1 atm for 1:1, 1:2, 1:5, 1:10. Predict how the pressure scales. *Finding:* p = 0.767, 1.534, 3.835, 7.671 bar (exactly proportional to k; a full-size wing would need less than atmospheric pressure, because sea-level air is denser than 3000 m air); cryogenic T = 359.7, 209.9, 110.2, 70.6 K (the last two indicative only). The model's dynamic pressure in the pressurised tunnel is 1790, 3581, 8952, 17 905 Pa, 1.1 to 11 times the flight value of 1636 Pa: loads, not just air supply, limit pressurised testing.
5. **Water.** The Re-matching speed in 20 °C water at the four scales. *Finding:* 3.23, 6.47, 16.16, 32.33 m/s. Water tunnels are excellent for flow visualisation at the right Reynolds number, but cavitation, free-surface effects and the absence of compressibility limit what they can say about an aeroplane.

Raw data format: CSV tables `relation, n, rank, groups_sympy, groups_repeating, same_span` (1); `trial, n, rank, agree` (2); `scale, chord_m, V_Re_ms, M_tunnel, Re_at_M` (3); `scale, p_bar, q_Pa, T_cryo_K` (4); `scale, V_water_ms` (5); plus `similarity.png`.

## 5. Validate and improve

| quantity | analysis (§2) | code | difference | cause |
|---|---|---|---|---|
| Pi count, wing | n − rank = 3 | 3 by all three routes | 0 | exact |
| dimensionless check | De = 0 | max \|DB\| = 3.1 × 10⁻¹⁶ (SVD), 0 (SymPy) | round-off | — |
| span agreement | same nullspace | 200/200 random matrices | 0 | exact |
| Mach of a Re-matched 1:k model | k × M_full | k × M_full to 10⁻¹² | round-off | — |
| pressurised 1:5 | Re = 4.8308 × 10⁶, M = 0.1826 | substituted back: same | < 10⁻¹² | — |
| Stokes collapse | C_D Re = 24 | 24 in three fluids to 10⁻¹² | round-off | — |
| Re in three unit systems | identical | agree to 10⁻⁶ (imperial factors rounded) | 10⁻⁶ | conversion constants |

The engineering answer: **at full scale the wing works at Re = 4.83 × 10⁶ and M = 0.183. A 1:5 model in an atmospheric sea-level tunnel can match one but not both: the Reynolds number at 235 m/s (M = 0.69), or the Mach number at 62 m/s (Re = 1.28 × 10⁶, a quarter of flight). Matching both takes a pressurised tunnel at 3.84 bar at 288 K, at 5.5 times the flight dynamic pressure on the model, or a cryogenic tunnel near 110 K at 1 atm (indicative).** For a low-speed wing the usual choice is to match Mach, which is easy, and correct for Reynolds number with transition trips and data at several Reynolds numbers, because the compressibility error of the alternative is larger and harder to correct.

What to improve: a check that the variable list is complete, which no algebra provides; in practice, a test at two Reynolds numbers shows whether the Re-dependence is weak enough to extrapolate. What not to improve: the nullspace algorithm; three independent routes agree exactly.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivation of (2), of (3) by repeating variables and of the Stokes limit; `similarity.png`; the experiment tables; and a half-page test plan (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) for a model of the wing at your scale, saying which group you match, which you do not, and how you correct for it.

### Adversarial self-checks

- **What if the units change?** Nothing dimensionless may change: the tests compute Re in SI, cgs and imperial units and M in knots, and they agree. A group that changes with units is not dimensionless, and the error is in the dimension table.
- **What if a variable is missing?** The algebra still returns groups, and they are wrong. Leaving a out of the wing relation gives two groups and no Mach number; nothing in the code can notice. Physics, not linear algebra, chooses the list.
- **Which assumption fails first in the real world?** Geometric similarity of the surface and of transition: a model at a quarter of the flight Reynolds number keeps a laminar boundary layer longer, so its drag and maximum lift differ even when Mach is matched. Trip strips force transition at a chosen place.
- **What measurement could falsify the scaling?** Tests of the same model at two Reynolds numbers (two tunnel pressures) at one Mach number: if C changes with Re in the tested range, extrapolation to flight Re carries that trend's uncertainty; if C did not depend only on Re and M, two different (ρ, V, c, μ) with the same groups would give different C, which a pressurised tunnel can test directly.
- **Does the result survive a change of method?** SymPy, SVD and repeating variables give the same space; the tunnel conditions are substituted back.
- **Which decision changes if the uncertainty doubles?** If Sutherland's law is uncertain by 10 % at cryogenic temperatures instead of 5 %, the cryogenic tunnel temperature for a Re match moves by several kelvin; a test campaign planned on the edge of a tunnel's temperature limit would have to be redesigned around the pressurised option.

## Check yourself

1. State the Buckingham Pi theorem: how many independent dimensionless groups does a relation have?
2. What are the columns of the dimension matrix, and where in it does a dimensionless group live?
3. How do the SymPy and SVD bases of the groups differ, and why are both correct?
4. Define the Reynolds and Mach numbers and say what each compares.
5. Why can a small model in the same air not match both Re and M?
6. Name three ways a wind tunnel can raise the Reynolds number of a model at fixed Mach number.
7. Why do Reynolds and Mach numbers come out the same in any system of units?
8. What is partial similarity, and give an example of a test that relies on it.
