# P05.15  Column Buckling

> A strut in compression does not wait for the material to fail. At a load set by its stiffness and its length, straight stops being the only way to stand, and it bows.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the Euler critical load of a pinned column from the equilibrium of its bent shape, and the effective-length factors K = 2, 0.699 and 0.5 of the fixed-free, fixed-pinned and fixed-fixed columns from the fourth-order equation and its boundary conditions. 2. Derive the Johnson parabola from the condition that it meets the Euler curve with the same value and slope, and use the resulting column curve to size a strut with an aerospace margin of safety. 3. Compute the tangent-modulus (Engesser) critical stress from a Ramberg–Osgood compression curve and explain what the Johnson parabola approximates. 4. Implement the buckling load as the smallest eigenvalue of a finite-difference generalised eigenproblem with ghost nodes and show that it converges at second order to every closed form. |
| **Prerequisites** | `P05.05` (second moment of area, radius of gyration), `P05.06` (the stress–strain curve, compressive yield), `P05.12` (E I w'' = M and its boundary conditions), `P05.14` (margins of safety, the factor 1.5), `P02.02` (second-order linear ODEs with constant coefficients), `P01.16` (eigenvalues and eigenvectors), `P03.13` and `P03.08` (assembling and solving discretised equations), `P03.10` (root finding). Skill check: you can solve y'' + k²y = 0 with y(0) = y(L) = 0 and say for which k a non-zero solution exists. |
| **Engineering question** | A 4130 steel engine-mount strut, a round tube 22.2 mm in outside diameter with a 0.889 mm wall and 500 mm between pinned ends, carries 10 kN of compression at limit load. What load buckles it, what is its margin of safety at ultimate load, and how wrong would the Euler formula be? The answers: slenderness L/r = 66.3, below the Johnson–Euler transition at 96.4; Johnson gives 19.8 kN and MS = +0.32; Euler would give 27.4 kN and MS = +0.83, overestimating the strength by 39 %. |
| **Assumptions and limits** | A straight, prismatic, homogeneous column loaded exactly along its centroidal axis; small lateral deflections (the buckling load is the bifurcation point, not the post-buckled shape); ideal end conditions (real joints sit between pinned and fixed, so design uses K = 0.65 rather than 0.5 and 0.8 rather than 0.699); the Johnson parabola is an empirical fit, and the tangent-modulus stress depends on a Ramberg–Osgood exponent that is assumed here (n = 15), not measured. Initial crookedness, load eccentricity and local buckling of a thin wall all lower the real strength and are not modelled; the lesson says how big each can be. |
| **Analytical work** | The pinned Euler load from E I w'' + P w = 0; the general equation E I w'''' + P w'' = 0 and its boundary conditions for pinned, fixed and free ends; K for each case, including tan β = β for fixed-pinned; slenderness and radius of gyration; the Johnson parabola from tangency; the tangent-modulus equation σ = π²E_t(σ)/(KL/r)²; the strut margin MS = P_cr/(1.5 P_limit) − 1. |
| **Computational deliverable** | `starter/main.py`: `effective_length_factor`, `euler_load`, `euler_stress`, `transition_slenderness`, `johnson_stress`, `critical_stress`, `strut_margin`, `ramberg_osgood_strain`, `tangent_modulus`, `tangent_modulus_stress`, `fd_buckling`, `euler_mode_residual` (`tube_properties` and `main` are given). Tests prove: the Euler load scales as EI/(KL)² and the end conditions stand in the ratios ¼ : 1 : 2.046 : 4; tan β = β for the fixed-pinned K; the finite-difference eigenvalue converges at order 1.8–2.3 to the closed form for all four end conditions; the pinned mode is a half sine that satisfies E I w'' + P w = 0 at P_cr and not at 1.2 P_cr; the fixed-free mode is 1 − cos(πx/2L); Johnson starts at F_cy, meets Euler at the transition with equal value and slope, and never exceeds Euler or F_cy; the tangent modulus matches the derivative of the Ramberg–Osgood curve; the tangent-modulus stress solves its own equation, lies below Euler and meets it for slender columns. |
| **Practical activity** | **computational experiment**: a convergence study of the finite-difference buckling eigenvalue for four end conditions; the column curve of the engine-mount tube against length (Euler, Johnson, tangent modulus); the transition slenderness of four aerospace alloys; the sensitivity of the tangent-modulus curve to the Ramberg–Osgood exponent and of the strut's strength to the effective-length factor. |
| **Equipment** | None. A computer with the environment from `P00.01`. (Optional demonstration: press down gently on a plastic ruler standing on a table and watch it snap sideways at a definite load; safety glasses, and use a flexible ruler, not a brittle one.) |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `column_buckling.png`, the derivations of the Euler load and of the Johnson tangency, the four experiment tables, and a half-page strut sizing note. |
| **Assessment** | `assessments/phase-05/grading-rubric.md`; competencies 1–4. Pass: the Euler load and K for one non-pinned case are derived; all public tests pass; the observed convergence order is reported; the note picks the right branch of the column curve and states the margin with its basis. |
| **References** | `@mit-ocw-16-001-unified-materials-structures` lecture notes "Buckling of Beams" (lectures 33–35): Introduction to bifurcation of equilibrium and structural instability; A single degree-of-freedom example; Columns; Column buckling (Euler buckling); Example: Simply-supported column; Other boundary conditions; Effects of initial imperfections; Failure of columns (required-free). `@libretexts-wierzbicki-structural-mechanics` Ch. 8 "Stability of Elastic Structures": 8.3 Stability of Elastic Column Using the Energy Method; 8.4 Effect of Structural Imperfections; 8.6 Plastic Buckling of Columns (supplementary-free): the energy view, imperfections, and the tangent modulus. `@mit-ocw-3-11-mechanics-of-materials` Roylance, module "Stresses in Beams": section "Buckling" (required-free). `@nasa-tm-x-73305-astronautic-structures-manual` Section C1 "COLUMNS" and A3.3.0 Interaction for Beam-Columns (supplementary-free): how aerospace column allowables are organised. |
| **Extensions** | ★ Add a variable EI(x) to `fd_buckling` and find the buckling load of a tapered strut; compare with a uniform strut of the same mass. ★ Add an initial crookedness w₀ sin(πx/L) to a pinned column and plot the lateral deflection against P/P_cr (the Southwell plot); show how to extract P_cr from a test that never reaches it. ★ Solve the same eigenproblem with the energy (Rayleigh–Ritz) method using one and two sine terms for the fixed-pinned column and compare with K = 0.699. |
| **Time** | unmeasured |

## 1. The engineering problem

An engine mount is a small space frame of steel tubes that carries a piston engine's weight, thrust and torque into the firewall. Under a hard landing or a high-g pull-up, some of its tubes are in compression. A tube in tension fails when its stress reaches the ultimate strength; a tube in compression can fail at a much lower stress by **buckling**: at a definite load, the straight shape stops being stable, a small sideways disturbance grows instead of dying away, and the tube bows out of line and stops carrying load. The same failure governs landing-gear drag braces, wing lift struts, fuselage longerons and stringers between frames, and the truss members of a spacecraft structure. It is why compression members are short, fat and hollow, and why a strut's strength depends on how its ends are held.

The classical answer is Euler's: a pinned column of length L and bending stiffness EI buckles at P_cr = π²EI/L². It is exact for a perfect elastic column and it is wrong, dangerously, for a short one, because a short column reaches the material's compressive yield stress before the Euler load. Aerospace practice uses a **column curve**: Euler for slender columns, and for intermediate ones a curve that bends over to the compressive yield stress F_cy. The **Johnson parabola** is the standard one for metallic struts, and the physics behind it is the **tangent modulus**: once the stress is in the non-linear part of the stress–strain curve, the column's effective stiffness is the local slope E_t, not E.

Success criteria, written first: the finite-difference eigenvalue reproduces the Euler load for four end conditions at second order; Johnson meets Euler with the same slope at the transition slenderness; the strut's buckling load, its margin, and the size of Euler's error are reported with the branch of the curve that applies. A strut sized by Euler in the Johnson range is an under-designed strut: in the engine-mount example, 39 % weaker than its designer believes.

## 2. Understand and derive

### 2.1 The bent shape in equilibrium

Take a straight pinned column of length L carrying an axial compression P, and suppose it has bowed a small amount w(x). Cut it at x and look at the lower piece. The axial force at the end is P, but at the cut the column has moved sideways by w, so P has a lever arm w about the cut: the internal moment is **M = −P w** (it bends the column back toward straight when w > 0). With the moment–curvature law of `P05.12`, E I w'' = M:

**E I w'' + P w = 0,  w(0) = w(L) = 0.**   (1)

With k² = P/(EI) the general solution is w = A sin kx + B cos kx. w(0) = 0 gives B = 0; w(L) = 0 gives A sin kL = 0. Either A = 0 (the straight column, always an equilibrium) or **sin kL = 0**, which has non-zero solutions only for kL = π, 2π, … The smallest load at which a bent equilibrium exists is

**P_cr = π² E I / L²,  w = A sin(πx/L).**   (2)

This is a **bifurcation**: below P_cr the straight shape is the only equilibrium; at P_cr a second, bent one appears with an undetermined amplitude A (the linear theory cannot fix A; that is the post-buckling problem). It is also an **eigenvalue problem**: P/EI is an eigenvalue of −d²/dx² with those boundary conditions, which is exactly how part 3 will compute it.

### 2.2 Other ends: the fourth-order equation

For other supports the moment at a cut also contains unknown end reactions. Differentiating twice removes them:

**E I w'''' + P w'' = 0,**   (3)

with general solution w = A sin kx + B cos kx + Cx + D and two boundary conditions at each end:

| end | conditions | meaning |
|---|---|---|
| pinned | w = 0, w'' = 0 | no deflection, no moment |
| fixed | w = 0, w' = 0 | no deflection, no rotation |
| free | w'' = 0, E I w''' + P w' = 0 | no moment, no transverse shear (the load stays vertical while the end rotates) |

Each case gives four homogeneous equations in A, B, C, D, and a non-zero solution needs their determinant to vanish. The **fixed-free** column (a flagpole, or a strut cantilevered off a frame) gives cos kL = 0, kL = π/2; the **fixed-fixed** column gives kL = 2π (the antisymmetric mode needs a larger load); the **fixed-pinned** column gives **tan kL = kL**, whose smallest positive root is β = 4.4934. In every case the result can be written as Euler's formula with an **effective length** KL, the length of the equivalent pinned column:

**P_cr = π² E I / (K L)²,**   (4)

| ends | K (ideal) | P_cr / P_pinned | K used in design |
|---|---|---|---|
| fixed-free | 2 | 0.25 | 2.1 |
| pinned-pinned | 1 | 1 | 1 |
| fixed-pinned | π/4.4934 = 0.6992 | 2.046 | 0.8 |
| fixed-fixed | 0.5 | 4 | 0.65 |

No real joint is perfectly rigid, so design practice uses the larger K values of the last column; a welded tube cluster on an engine mount is often taken as pinned outright, which is conservative.

### 2.3 Slenderness and the stress at buckling

Divide (4) by the area A and write I = A r², with **r = √(I/A)** the radius of gyration:

**σ_cr = π² E / (K L / r)².**   (5)

The **slenderness ratio** KL/r is the one number that decides a column's behaviour. For a thin-walled round tube r ≈ d/(2√2): for the engine-mount tube (d_o = 22.2 mm, t = 0.889 mm) A = 59.52 mm², I = 3385 mm⁴, r = 7.54 mm, and at L = 500 mm, **L/r = 66.3**. Euler's stress there is π² × 205 GPa / 66.3² = **460 MPa**: above 4130's typical compressive yield of 435 MPa. The Euler formula is predicting that the tube buckles *elastically* at a stress the material cannot carry elastically. Something else must happen first.

### 2.4 The Johnson parabola

For intermediate columns the standard aerospace curve is a parabola in the slenderness that starts at F_cy for a zero-length column and joins the Euler curve smoothly:

σ_J = F_cy − C (KL/r)².

Choose C and the joining slenderness λ_c so that the two curves have the **same value and the same slope** at λ_c (a smooth handover, no kink in the design curve):

F_cy − C λ_c² = π²E/λ_c²,  −2Cλ_c = −2π²E/λ_c³.

The second gives C = π²E/λ_c⁴; putting it in the first, F_cy = 2π²E/λ_c², so

**λ_c = √(2π² E / F_cy),  σ_J = F_cy − (F_cy/(2π))² (KL/r)² / E,**   (6)

and at the transition both curves give **F_cy/2**. For 4130 (typical E = 205 GPa, F_cy = 435 MPa) λ_c = 96.4; for 2024-T3 (73.1 GPa, 345 MPa) 64.7; for 7075-T6 (71.7 GPa, 503 MPa) 53.0. The strong aluminium alloys go "Johnson" at lower slenderness because their yield is high relative to their stiffness. Below λ_c use Johnson, above it Euler: that is the column curve, `critical_stress`.

For the engine-mount strut at L/r = 66.3: **σ_J = 332 MPa, P_cr = 19.77 kN**, against Euler's 27.39 kN. Euler overestimates the strength by 39 %.

### 2.5 Why the curve bends over: the tangent modulus

A metal's compressive stress–strain curve is not a straight line up to F_cy and then flat; it bends gradually. A **Ramberg–Osgood** fit in the yield-anchored form is

**ε = σ/E + 0.002 (σ/F_cy)ⁿ,**   (7)

which gives exactly 0.2 % plastic strain at F_cy (the offset yield of `P05.06`); the exponent n sets how sharp the knee is (large n: sharp, like mild steel; n ≈ 10–20: gradual, like most aerospace alloys). Its local slope, the **tangent modulus**, is

E_t = dσ/dε = 1 / (1/E + 0.002 n σ^(n−1) / F_cy^n).

Engesser's argument: a column that starts to bend while the load is still increasing has every fibre still loading, so its stiffness against bending is E_t, and

**σ_t = π² E_t(σ_t) / (KL/r)²**   (8)

is an equation for σ_t (solve it with `brentq`). Shanley showed that this tangent-modulus load is the load at which a real, perfect column starts to bend; it is the right lower bound for design. With n = 15 assumed for 4130 (illustrative; a measured compression curve would replace it), the engine-mount strut gives σ_t = 335 MPa, within 1 % of Johnson's 332 MPa. At that stress E_t = 155 GPa, three quarters of E: the strut buckles in the knee of the curve. The Johnson parabola is, in effect, a one-parameter stand-in for (8) that needs only E and F_cy.

### 2.6 The strut's margin

With a limit compression of 10 kN and the ultimate factor 1.5 (`P05.14`), the margin is

**MS = P_cr / (1.5 P_limit) − 1 = 19.77/15 − 1 = +0.32.**

By Euler it would be +0.83. In design terms the tube passes, with room; and the strut could be lengthened to 667 mm before the margin reached zero.

### 2.7 Where the models stop

- **Imperfections.** A column with initial crookedness or an eccentric load starts bending at once and has no sharp bifurcation; its deflection grows like 1/(1 − P/P_cr) and the material yields before P_cr. Aerospace column curves are calibrated on tests of real struts, which include typical imperfections (MIT 16.001 "Effects of initial imperfections"; Wierzbicki 8.4).
- **Local buckling.** A thin-walled tube can wrinkle locally (the wall crumples) before the whole tube bows. For round tubes, d/t below about 50–60 keeps local buckling out of the way in this load range; the engine-mount tube has d/t = 25.
- **Ends.** The K values are idealisations; the real one comes from the stiffness of the joint and is often found by test or by finite-element analysis of the assembly.
- **Torsional and flexural-torsional buckling** of open sections (channels, angles) can govern stringers; a round tube is immune.
- **The Ramberg–Osgood exponent** in this lesson is assumed. The tangent-modulus curve is sensitive to it; Part 4 measures how sensitive.

## 3. Build the computational model

Open `starter/main.py`. `tube_properties` and `main` are given.

### Step 1: Euler and the effective length

`effective_length_factor(ends)` returns K from the table (the fixed-pinned value as π/4.493409457909064). `euler_load(E, I, L, K)` is (4); `euler_stress(E, slenderness)` is (5).

### Step 2: Johnson and the column curve

`transition_slenderness` and `johnson_stress` are (6); `critical_stress` picks Johnson below λ_c and Euler above it (use `np.where`, so it works on arrays); `strut_margin` is the margin of 2.6.

### Step 3: the tangent modulus

`ramberg_osgood_strain` is (7); `tangent_modulus` its inverse slope (write it as (σ/F_cy)^(n−1) to avoid overflow for large n); `tangent_modulus_stress(slenderness, E, F_cy, n)` solves (8) with `brentq` between a tiny stress and min(Euler stress, 3F_cy).

### Step 4: the finite-difference eigenproblem

`fd_buckling(L, EI, ends, n)` discretises (3) on n intervals of h = L/n, with nodes x₀ … x_n and **two ghost nodes beyond each end**, n + 5 unknowns. At every real node the central differences

w'''' ≈ (w_{i−2} − 4w_{i−1} + 6w_i − 4w_{i+1} + w_{i+2})/h⁴,  w'' ≈ (w_{i−1} − 2w_i + w_{i+1})/h²

turn (3) into a row of **A w = λ B w** with λ = P/EI, A holding the fourth difference and B minus the second. Four more rows impose the boundary conditions with central differences (w' ≈ (w_{i+1} − w_{i−1})/2h, w''' ≈ (w_{i+2} − 2w_{i+1} + 2w_{i−1} − w_{i−2})/2h³). Their B rows are zero, except the free end's shear condition, whose P w' term belongs in B:

```python
A[row, k-2:k+3] = np.array([-1.0, 2.0, 0.0, -2.0, 1.0]) / (2.0 * h**3)   # w'''
B[row, k-1:k+2] = -np.array([-1.0, 0.0, 1.0]) / (2.0 * h)               # -(P/EI) w'
```

`scipy.linalg.eig(A, B)` then returns some finite eigenvalues and some infinite ones, because B is singular: its boundary rows contain no P, and it never touches the two outermost ghost nodes. Keep the finite, real, positive eigenvalues and take the smallest, times EI. The eigenvector at the real nodes is the buckled shape. `euler_mode_residual` checks a pinned mode against (1).

Running `main()` prints

```
4130 strut 22.2 x 0.889 mm, L = 500 mm: A = 59.52 mm^2, I = 3385 mm^4, r = 7.54 mm, L/r = 66.3; transition slenderness 96.4
pinned ends: Euler 27.39 kN (460 MPa), Johnson 19.77 kN (332 MPa); Euler overestimates by 39 %
limit compression 10.0 kN, ultimate 15.0 kN: MS (Johnson) = +0.318, MS (Euler, wrong here) = +0.826
tangent modulus (Ramberg-Osgood n = 15): 19.97 kN (335 MPa)

ends            K       closed form [kN]   FD n=50     FD n=100    FD n=200    observed order
pinned-pinned  1.0000        27.3931      27.3841     27.3909     27.3926   2.00
fixed-free     2.0000         6.8483       6.8477      6.8481      6.8483   2.03
fixed-pinned   0.6992        56.0395      55.9943     56.0282     56.0367   2.00
fixed-fixed    0.5000       109.5726     109.4284    109.5365    109.5632   1.94
```

followed by the experiment tables, and writes `column_buckling.png` (the 4130 column curve with Euler, Johnson and tangent modulus and the strut marked; the four first modes from the eigenproblem).

### Verify

Run `python -m pytest curriculum/phase-05/lesson-15-column-buckling`. The tests:

1. **Euler scaling**: P ∝ E, ∝ I, ∝ 1/L², and K = 2 quarters it.
2. **End conditions**: K values, tan β = β for fixed-pinned, the ratio 2.046.
3. **Stress form**: σ_cr A = P_cr; the tube's radius of gyration.
4. **FD convergence** (independent check), for each of the four end conditions: observed order between 1.8 and 2.3 from n = 40 to 160, error below 5 × 10⁻⁴ at n = 160.
5. **Pinned mode**: a half sine to 10⁻⁶, zero at both ends, and in equilibrium (E I w'' + P w = 0 to 10⁻³) at P_cr but not at 1.2 P_cr.
6. **Fixed-free mode**: zero slope at the base, maximum at the tip, equal to 1 − cos(πx/2L).
7. **Linearity**: the FD load is proportional to EI.
8. **Johnson**: equal value and slope with Euler at λ_c, Euler stress F_cy/2 there, F_cy at zero slenderness.
9. **Column curve**: decreasing, below F_cy and below Euler, continuous across λ_c.
10. **Margin** definition.
11. **Ramberg–Osgood**: 0.2 % offset at F_cy; E_t equals the numerical derivative; E_t → E at low stress and non-increasing.
12. **Tangent modulus stress**: solves (8), stays below Euler, equals Euler for slender columns.
13. **main** writes the figure.

### Compare

There is no library call for "buckling load" in NumPy or SciPy, and the lesson does not use a structural code; the established tool here is the eigenvalue solver `scipy.linalg.eig`, applied to an operator you assembled. Its answer converges to the closed forms at second order; Richardson extrapolation of the pinned values at n = 50 and 100, (4P₁₀₀ − P₅₀)/3, gives 27.3931 kN, the exact value to five figures. A finite-element code would use beam elements with a geometric stiffness matrix and solve the same kind of generalised eigenproblem (K − λK_G)φ = 0; its first eigenvalue for one of these columns would agree with yours to the discretisation error.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction first.

1. **Convergence for four end conditions.** `fd_buckling` at n = 50, 100, 200 for each case; relative error and observed order. Predict the order (the stencils are second order) and predict which case converges slowest relative to its own load (the stiffest, whose mode has the most curvature per length).
2. **The column curve against length.** For the engine-mount tube, L = 150 to 1200 mm: Euler, Johnson (only where it applies), tangent modulus (n = 15), and the column curve. Predict at which length Euler and the column curve meet.
3. **Transition slenderness of four alloys.** λ_c for 2024-T3, 7075-T6 and 4130 from (6). Predict their order from F_cy/E.
4. **Sensitivity to the Ramberg–Osgood exponent.** σ_t at slenderness 30, 46.4, 66.3 and 90 for n = 10, 15, 25, 50, against Johnson. Predict the limit as n → ∞ (an elastic–perfectly plastic material).
5. **The effective-length factor.** The strut's column-curve load for K = 0.65, 0.8, 1.0 and 2.0. Predict the factor between the strongest and weakest.

Raw data format: CSV tables `ends, K, n, P_fd, P_exact, rel_error` (1); `L_mm, slenderness, euler_MPa, johnson_MPa, tangent_MPa, curve_MPa` (2); `material, E_GPa, Fcy_MPa, lambda_c` (3); `n, slenderness, sigma_t_MPa, johnson_MPa` (4); `K, slenderness, sigma_MPa, P_kN` (5); plus `column_buckling.png`.

What you should find. (1) The table of part 3: orders 1.94–2.03, the fixed-fixed case pre-asymptotic at n = 50 because its mode packs a full wave into the length. (2):

```
L [mm]   L/r     Euler [MPa]  Johnson [MPa]  tangent-modulus [MPa]  column curve [MPa]
   150    19.9         5114            426                    427                 426
   250    33.2         1841            409                    395                 409
   350    46.4          939            385                    371                 385
   500    66.3          460            332                    335                 332
   700    92.8          235            234                    234                 234
   900   119.3          142          (n/a)                    142                 142
  1200   159.1           80          (n/a)                     80                  80
```

Euler and the column curve meet at the transition, L = λ_c r = 727 mm; below it Euler is useless (at 250 mm it promises 4.5 times the real strength). (3) λ_c = 53.0 (7075-T6), 64.7 (2024-T3), 96.4 (4130). (4):

| n | σ_t at 30 | 46.4 | 66.3 | 90 [MPa] |
|---|---|---|---|---|
| 10 | 401.6 | 357.8 | 312.1 | 239.4 |
| 15 | 401.5 | 371.2 | 335.5 | 248.4 |
| 25 | 406.1 | 387.0 | 361.3 | 249.8 |
| 50 | 414.5 | 404.4 | 388.5 | 249.8 |
| Johnson | 414.0 | 384.7 | 332.2 | 245.6 |

As n grows the tangent-modulus curve tends to the elastic–perfectly plastic "min(F_cy, Euler)" with a corner; a gradual knee (small n) pulls it below Johnson in the middle of the range. At L/r = 46 Johnson is 3.6 % *unconservative* against n = 15: a reminder that an empirical curve should be checked against the real material's compression curve. (5) K = 0.65, 0.8, 1.0, 2.0 give 23.3, 22.0, 19.8 and 6.85 kN: a factor 3.4 between a well-fixed and a cantilevered strut, and only 18 % between the design "fixed" value and pinned, because at these slenderness ratios Johnson flattens the curve.

## 5. Validate and improve

| quantity | theory (§2) | code | difference | cause |
|---|---|---|---|---|
| pinned Euler load | π²EI/L² = 27.3931 kN | FD n = 200: 27.3926 kN | 2.1 × 10⁻⁵ | h² truncation |
| Richardson (n = 50, 100) | 27.3931 kN | 27.3931 kN | < 10⁻⁷ | — |
| fixed-free / fixed-pinned / fixed-fixed | ¼, 2.046, 4 × pinned | FD n = 200 within 5 × 10⁻⁶, 5 × 10⁻⁵, 9 × 10⁻⁵ | as listed | h² truncation, orders 1.94–2.03 |
| pinned mode shape | sin(πx/L) | FD eigenvector | < 10⁻⁶ | — |
| Johnson–Euler tangency | equal value and slope at λ_c = 96.4 | same to 10⁻⁶ | — | — |
| stress at transition | F_cy/2 = 217.5 MPa | 217.5 MPa | 0 | — |
| strut, Johnson against tangent modulus | — | 332 against 335 MPa | 1 % | assumed n = 15 |

The finite-difference eigenproblem is a genuinely independent route: it knows nothing of sines, cosines or tan β = β, only the differential equation and the boundary conditions, and it reproduces all four closed forms. The column curve has no such independent check inside this lesson, because it is an empirical design rule; its validation is the strut tests behind it, which is why Section C1 of the Astronautic Structures Manual tabulates column allowables per material and product form. The engineering answer is **P_cr = 19.8 kN by Johnson, MS = +0.32 at 10 kN limit; Euler's 27.4 kN does not apply at L/r = 66**.

What to improve: K for the real welded joint (a finite-element model of the cluster, or a test); the Ramberg–Osgood exponent from a measured compression curve (`P05.06`); an imperfection allowance. What not to improve: the FD grid beyond n = 100, whose error (at most 3 × 10⁻⁴) is already about a thousand times smaller than the end-fixity uncertainty.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (2), of K for one other end condition, and of (6); `column_buckling.png`; the five experiment tables; and a half-page strut sizing note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7): loads, geometry, section properties, slenderness, which branch of the column curve applies and why, P_cr, the margin, and the largest uncertainty.

### Adversarial self-checks

- **What if the load doubles?** The buckling load does not depend on the applied load; the margin falls from +0.32 to 19.77/30 − 1 = −0.34. A strut either buckles or does not; it does not "half buckle", so there is no warning before a negative-margin strut fails.
- **What if the length goes to zero?** Euler goes to infinity and is meaningless; Johnson goes to F_cy, the squash load of a stub. What if it goes to infinity? Both go to zero as 1/L², and the strut cannot even carry its own weight horizontally; the column curve is Euler there.
- **Which assumption fails first in the real world?** The ideal end condition and the perfect straightness: a 1 mm crookedness over 500 mm and a few millimetres of joint eccentricity are normal, and both lower the load. That is why test-based column curves sit below the theoretical ones.
- **What measurement could falsify the model?** A compression test of the strut with a sideways dial gauge at mid-length: the Southwell plot (Extension 2) gives P_cr from the pre-buckling data; a value far above Johnson would say the joint is stiffer than pinned, far below would say imperfection or local buckling.
- **Does the result survive a change of discretisation or method?** The FD load changes by at most 2.4 × 10⁻⁴ relative from n = 100 to 200 (fixed-fixed); the column curve and the tangent-modulus curve agree within 1 % at the design point and within 4 % across the Johnson range for n = 15.
- **Which decision changes if the end fixity is uncertain?** If the joint were credited as fixed (K = 0.65), the strut would appear 18 % stronger and someone might lengthen or lighten it; if it then behaved as pinned, the margin would be overestimated by that much. Take pinned unless a test says otherwise.

## Check yourself

1. Write the Euler load of a pinned column and explain in one sentence why it does not depend on the material's strength.
2. What are the effective-length factors of a fixed-free and a fixed-fixed column, and how do their buckling loads compare with a pinned one of the same length?
3. Define the slenderness ratio and give its value for a tube with r = 7.5 mm and an effective length of 600 mm.
4. At what slenderness does the Johnson parabola meet the Euler curve, and what is the critical stress there?
5. A strut's slenderness is below the transition value. Which formula do you use, and what happens if you use the other?
6. What is the tangent modulus, and why does it, rather than E, govern a column that buckles in the knee of the stress–strain curve?
7. A strut has P_cr = 24 kN and a limit load of 12 kN. What is its margin of safety?
8. Why does the finite-difference buckling problem produce infinite eigenvalues, and what do you do with them?
