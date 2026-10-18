# P05.12  Beam Deflection

> The bending moment sets the curvature; integrate the curvature twice and the boundary conditions tell you where the beam went.

| | |
|---|---|
| **Learning outcomes** | 1. Derive EI w″ = M and (EI w″)″ = q from the results of `P05.10` and `P05.09`, and state the boundary conditions of clamped, pinned and free ends. 2. Integrate the curvature by hand for the standard cantilever and simply supported cases, and combine them by superposition, including to solve a statically indeterminate propped cantilever and the elliptically loaded wing. 3. Implement deflection by numerical double integration of M/EI for any load and any spanwise stiffness, and an independent finite-difference solution of the fourth-order beam equation for any of four support conditions. 4. Show that the finite-difference solution converges at second order for smooth loads, explain why it converges more slowly for an elliptic lift, and check superposition and Maxwell's reciprocal theorem numerically. 5. Decide whether a wing's tip deflection at limit load is acceptable and what changes it most. |
| **Prerequisites** | `P05.10` (curvature κ = M/(EI), the spar section and its EI), `P05.09` (bending-moment diagrams, singularity functions, sign conventions), `P02.02` (second-order linear ODEs and boundary conditions), `P03.11` (finite differences, truncation error and observed order), `P03.08` (direct linear solvers), `P00.14`. Skill check: you can integrate a polynomial twice and fix the two constants from two conditions. |
| **Engineering question** | At its 3.8 g limit load, how far does the tip of the light aeroplane's 5 m wing rise if the spar of `P05.10` (EI = 6.96 × 10⁵ N·m²) carries the elliptic lift alone, how much do the wing's weight and a 40 kg tank relieve it, and how much does tapering the spar add? The answers: 243 mm (4.9 % of the semi-span) with a tip slope of 3.6°; 187 mm with the weight and tank; 284 mm if EI falls linearly to 30 % at the tip. |
| **Assumptions and limits** | Euler–Bernoulli beam (plane sections, shear deformation neglected); linear elastic material; small deflections and slopes (w′² ≪ 1, so the curvature is w″); loads keep their direction as the beam deflects; supports are ideal (a clamp allows no rotation at all); EI may vary along the span but the section shape stays similar. The wing's skins and stringers, which add to EI, are left out, so the deflection is an upper bound for the spar-alone model. |
| **Analytical work** | The moment–curvature relation and its integration; boundary conditions; closed forms for cantilevers and simply supported beams under point loads, uniform loads and end moments; superposition; the propped cantilever's prop force 3qL/8; the elliptic-lift tip deflection L s³(1/8 − 4/(45π))/EI by superposing point loads; Maxwell's reciprocal theorem; the truncation error of the finite-difference scheme. |
| **Computational deliverable** | `starter/main.py`: `cantilever_point`, `cantilever_udl`, `cantilever_end_moment`, `simply_supported_point`, `simply_supported_udl`, `propped_cantilever_reaction`, `propped_cantilever_udl`, `bending_moment`, `deflection_by_integration`, `fd_beam`, `elliptic_tip` (`_cumtrapz`, `elliptic_lift`, `tapered_EI`, `observed_order` and `main` are given). Tests prove: each closed form satisfies EI w⁗ = q and its boundary conditions; textbook extremes (PL³/3EI, PL³/48EI, 5qL⁴/384EI); continuity of deflection and slope under a point load; Maxwell reciprocity; the propped cantilever by superposition agrees with finite differences; numerical integration reproduces the closed forms; finite differences converge at order 2.0 ± 0.2 for all four support conditions, handle variable EI, point loads and superposition, and converge more slowly (order 1.2–1.8) under the elliptic lift. |
| **Practical activity** | **computational experiment**: a convergence study of the finite-difference beam against the closed forms, the effect of the lift distribution's shape and of spar taper on tip deflection, the relief from a wing-mounted mass at different stations, and the sensitivity of tip deflection to the root fitting's stiffness. |
| **Equipment** | None. A computer with the environment from `P00.01`. (Lesson `P05.18` measures a cantilever's deflection at home and compares it with these formulas.) |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `wing_deflection.png`, the derivations of the cantilever closed forms, the propped cantilever and the elliptic tip deflection, the five experiment tables, and a half-page note on the wing's tip deflection with what changes it most. |
| **Assessment** | `assessments/phase-05/grading-rubric.md`; competencies 1–4. Pass: EI w″ = M is derived from `P05.10` and integrated with correct boundary conditions; all public tests pass; the observed convergence orders are reported for both a smooth load and the elliptic lift, with the reason for the difference; the note states the tip deflection and the assumption that matters most. |
| **References** | `@mit-ocw-3-11-mechanics-of-materials` Roylance module "Beam Displacements" (https://ocw.mit.edu/courses/3-11-mechanics-of-materials-fall-1999/resources/mit3_11f99_bdisp/): Introduction; Multiple integration; Superposition (required-free): EI v″ = M integrated with boundary conditions and the table of standard cases. `@mit-ocw-16-001-unified-materials-structures` lecture notes "Transverse Loading of Slender Structural Elements; Beam Theory" (required-free). `@libretexts-wierzbicki-structural-mechanics` Ch. 4 "Solution Method for Beam Deflections" (https://eng.libretexts.org/Bookshelves/Mechanical_Engineering/Structural_Mechanics_(Wierzbicki)/04%3A_Solution_Method_for_Beam_Deflections) (supplementary-free): the fourth-order equation and its boundary conditions at graduate depth. |
| **Extensions** | ★ Add a rotational spring at the root to `fd_beam` (a mixed boundary row EI w″ = k w′) and fit the spring stiffness that reproduces a measured tip deflection. ★ Solve the clamped–clamped and propped cases with point loads by superposition and singularity functions, and compare with `fd_beam` at the nodes. ★ Replace the small-slope curvature w″ by w″/(1 + w′²)^{3/2} and solve the large-deflection cantilever with `scipy.integrate.solve_bvp`; find the load at which the small-slope tip deflection is 5 % in error. |
| **Time** | unmeasured |

## 1. The engineering problem

A wing bends upward in flight. On an airliner the tips rise by metres at the top of a gust; on a light aeroplane, by a few centimetres in cruise and tens of centimetres at limit load. How much is acceptable is not a strength question: a wing can be strong enough and still too flexible. Deflection changes the angle of the outboard sections and so the lift distribution; it couples with twist into aeroelastic problems (divergence, aileron reversal, flutter); it moves control runs and fuel lines; and on the ground it decides whether a wing tip clears a hangar door. Stiffness requirements size many structures more than strength does, and every stiffness requirement starts with a deflection calculation.

The tools are the curvature relation from `P05.10`, κ = M/(EI), and the bending-moment diagram from `P05.09`. Curvature is the second derivative of the deflection, so the deflection is the curvature integrated twice, with the two constants of integration fixed by how the beam is supported. This lesson does that integration three ways: by hand for the standard cases, numerically for any load and stiffness, and by finite differences on the fourth-order beam equation, which needs no moment diagram at all and so checks the other two. It then applies them to the wing of `P05.10`.

Success criteria, written first: every closed form satisfies the beam equation and its boundary conditions; numerical integration matches the closed forms to 10⁻⁵; finite differences converge at second order for smooth loads; superposition of closed forms matches the finite-difference solution with all loads applied together; the propped cantilever's prop force is 3qL/8 by superposition and by finite differences.

## 2. Understand and derive

### 2.1 From curvature to the beam equation

Let w(x) be the deflection of the beam axis, positive up, x along the beam from the root. For small slopes (w′ ≪ 1) the curvature of the axis is w″, and the beam of `P05.10` with M positive sagging (concave up) has κ = M/(EI). So

**EI w″ = M.**   (1)

Equilibrium of a slice (`P05.09`) gives dV/dx = q and dM/dx = V with q the distributed load positive up and V the upward force to the left of a cut; together, M″ = q. Substituting (1),

**(EI w″)″ = q,  and for constant EI,  EI w⁗ = q.**   (2)

Equation (1) is a second-order ODE that needs M(x) first, which needs the reactions: fine for a statically determinate beam. Equation (2) is fourth-order and needs only the load, with four boundary conditions, two at each end; it works whether the beam is determinate or not. Units: N·m² × m⁻² = N·m, a moment.

The **boundary conditions** come from what each support prevents:

| end | conditions | why |
|---|---|---|
| clamped (built in) | w = 0, w′ = 0 | no deflection, no rotation |
| pinned (simple support) | w = 0, M = EI w″ = 0 | no deflection, free to rotate |
| free | M = EI w″ = 0, V = (EI w″)′ = 0 | nothing applies a moment or a force |

![Three beams and their boundary conditions: a cantilever clamped at the left (w = w′ = 0) and free at the tip (M = V = 0) under an upward load q, deflecting upward by w; a pinned–pinned beam with w = 0 and M = 0 at both ends; and a propped beam clamped at the left (w = w′ = 0) and pinned at the right (w = 0, M = 0).](figures/beam-boundary-conditions.svg)
*Figure 1. Each support prescribes two of w, w′, M and V: a clamp fixes deflection and slope, a pin fixes deflection and leaves the moment zero, and a free end carries neither moment nor shear.*

### 2.2 Integrating by hand: the cantilever

A cantilever clamped at x = 0, length L, uniform upward load q. From the free-body diagram of the piece outboard of x, M(x) = q(L − x)²/2 (positive: lift sags the cantilever the way it sags the wing in `P05.10`). Integrate (1) twice:

EI w′ = −q(L − x)³/6 + C₁,  EI w = q(L − x)⁴/24 + C₁x + C₂.

The clamp gives w′(0) = 0 ⇒ C₁ = qL³/6, and w(0) = 0 ⇒ C₂ = −qL⁴/24. Expanding,

**w = q x²(6L² − 4Lx + x²)/(24EI),  w_tip = qL⁴/(8EI),  θ_tip = qL³/(6EI).**   (3)

The same steps with a point load P at x = a give w = Px²(3a − x)/(6EI) up to a and a straight line beyond it (no moment acts beyond the load), with the tip value **PL³/(3EI)** when a = L. A tip moment M₀ gives constant curvature and w = M₀x²/(2EI). For a simply supported beam the constants come from w(0) = w(L) = 0: the midspan deflection under a uniform load is **5qL⁴/(384EI)** and under a central point load **PL³/(48EI)** (Roylance, "Beam Displacements", Multiple integration). The lesson's code has all of these; the tests check that each satisfies EI w⁗ = q and its boundary conditions.

**A worked example with a singularity function.** A simply supported beam of length L carries an upward point load P at x = a (b = L − a). The left reaction is −Pb/L (downward, to balance the moment about the right support), so with the singularity function ⟨x − a⟩ = max(x − a, 0) from `P05.09`,

M(x) = −Pb x/L + P⟨x − a⟩.

Integrating twice, EI w = −Pb x³/(6L) + P⟨x − a⟩³/6 + C₁x + C₂. w(0) = 0 gives C₂ = 0; w(L) = 0 gives C₁ = Pb(L² − b²)/(6L). For x ≤ a the bracket is zero and

**w = P b x (L² − b² − x²)/(6 L EI),**

which is `simply_supported_point`; for a = b = L/2 it gives PL³/(48EI) at midspan. One expression covers both sides of the load: that is what singularity functions buy.

Every deflection has the form (load) × (length)ⁿ/(EI): doubling the span of a cantilever under a uniform load multiplies its tip deflection by 16. Span is the strongest lever on stiffness, which is why long, slender wings need deep spars.

### 2.3 Superposition

Equation (2) is linear: the deflection under several loads is the sum of the deflections under each one alone, provided the deflections stay small and the supports do not change. This turns a table of standard cases into a tool for any combination (Roylance, "Beam Displacements", Superposition).

**An indeterminate beam.** A cantilever under uniform load q with a prop (a pin) at the free end has three unknown reactions and two useful equilibrium equations. Remove the prop and replace it by an unknown force R. The free cantilever's tip rises by qL⁴/(8EI); the force R alone moves it by −RL³/(3EI). The prop allows no tip deflection, so the two must cancel:

**R = 3qL/8,  w = q x²(3L² − 5Lx + 2x²)/(48EI),  M(0) = qL²/2 − RL = qL²/8.**   (4)

The compatibility condition (zero deflection at the prop), not equilibrium, gave the extra equation: the second of the phase's three ingredients doing its job. For L = 1.2 m, q = 2000 N/m and EI = 2 × 10⁴ N·m²: R = 900 N, a clamp moment of 360 N·m, and a largest deflection of 1.123 mm.

**The elliptically loaded wing.** A point load P at η on a cantilever of length s lifts the tip by Pη²(3s − η)/(6EI) and rotates it by Pη²/(2EI). Treat the elliptic lift w(η)dη as a continuum of point loads and add their effects:

w_tip = ∫₀ˢ w(η) η²(3s − η)/(6EI) dη,  θ_tip = ∫₀ˢ w(η) η²/(2EI) dη.

With w = (4L/πs)√(1 − (η/s)²) and the two integrals ∫₀¹ t²√(1 − t²) dt = π/16 and ∫₀¹ t³√(1 − t²) dt = 2/15:

**θ_tip = L s²/(8EI),  w_tip = (L s³/EI)(1/8 − 4/(45π)) = 0.09671 L s³/EI.**   (5)

For comparison, the same lift spread uniformly gives L s³/(8EI) = 0.125 L s³/EI, and concentrated at the elliptic centre of lift 4s/(3π) it gives 0.0773 L s³/EI. Superposition with a kernel like this, the deflection at one point due to a unit load at another, is the idea behind influence coefficients and flexibility matrices in aeroelasticity.

**Maxwell's reciprocal theorem.** The deflection at b due to a unit load at a equals the deflection at a due to a unit load at b. The cantilever formula shows it directly (Pa²(3b − a)/(6EI) for a < b both ways), and it holds for any linear elastic structure; the tests check it.

### 2.4 The wing at limit load

For the spar of `P05.10` (7075-T6, E = 71.7 GPa typical, I = 9.704 × 10⁶ mm⁴ at the root), EI = 6.958 × 10⁵ N·m². With L = 13 974.5 N on a 5 m semi-span, (5) gives **w_tip = 0.2428 m and θ_tip = 3.596°**: the tip rises 4.86 % of the semi-span. The uniform-lift estimate would have been 0.3138 m, 29 % more; the point-load estimate 0.1941 m, 20 % less. The shape of the load matters as much as its size.

Superposition then adds the inertia relief. At 3.8 g a 45 kg wing spread uniformly pushes down with −3.8 × 45 × 9.80665/5 = −335.4 N/m, lowering the tip by 37.7 mm; a 40 kg tank at 2 m pushes down with 1491 N, lowering it by 18.6 mm. Together: **0.2428 − 0.0377 − 0.0186 = 0.1866 m**, and the root moment falls from 29 655 to 22 481 N·m. A wing that carries its fuel and engines outboard is lighter for the same strength: that is span loading, and it is why airliners hang their engines on the wings.

### 2.5 Where the model stops

- **Shear deformation.** The web's shear strain adds roughly ∫V dx/(G A_web) = M₀/(G A_web) ≈ 2.3 mm to the tip of the spar, 0.9 % of the bending deflection (Problem 6). It matters for short, deep beams and sandwich panels.
- **Large deflection.** With a 3.6° tip slope, w′² is 4 × 10⁻³, and replacing the true curvature by w″ is an error of order 10⁻³. Airliner wings at 10° and more need the nonlinear form (Extension 3).
- **The supports.** A real root fitting rotates a little under load; Experiment 5 shows a modest root rotation can add 10–20 % to the tip deflection. In the home lab of `P05.18` the clamp's compliance is the first suspect when theory and measurement disagree.
- **The section.** Skins and stringers add to EI, sometimes doubling it; the deflection here is for the spar alone.
- **Loads that follow the structure.** Lift changes as the wing twists and bends (aeroelasticity), so the load is not fixed in reality; this lesson treats it as fixed.

## 3. Build the computational model

Open `starter/main.py`. `_cumtrapz`, `elliptic_lift`, `tapered_EI`, `observed_order` and `main` are given.

### Step 1: the closed forms

`cantilever_point`, `cantilever_udl`, `cantilever_end_moment`, `simply_supported_point`, `simply_supported_udl` are 2.2. Use `np.where` for the two pieces of the point-load formulas so they accept arrays. `propped_cantilever_reaction` and `propped_cantilever_udl` are (4) built by superposition from the first two functions, not typed in.

### Step 2: the bending moment

`bending_moment(x, q, point_loads, support)` integrates M″ = q twice from the left with `_cumtrapz`, adds P(x − a)₊ for each point load (a singularity function, `P05.09`), and then adds a linear function C₁(x − x₀) + C₀ chosen to satisfy the two moment conditions of the support: M = M′ = 0 at the free end of a cantilever, M = 0 at both ends of a simply supported beam.

### Step 3: deflection by integration

`deflection_by_integration(x, M, EI, support)` integrates κ = M/EI twice with `_cumtrapz`, starting from θ = w = 0 at the clamp. For a simply supported beam it starts with θ = w = 0 at the left end and then adds the rigid rotation c x that brings w back to zero at the right end. EI may be an array.

### Step 4: the independent check, finite differences

`fd_beam(L, EI, q, n, bc, point_loads)` does not use M(x) from statics. It solves (2) written as two second-order equations,

m″ = q,  w″ = m/EI,

on n + 1 equally spaced nodes, with the central difference (f_{i−1} − 2f_i + f_{i+1})/h² at interior nodes. That gives 2(n − 1) equations for 2(n + 1) unknowns (w_i and m_i); the four boundary conditions close the system. Conditions on w or m are one row each. A slope or shear condition uses the end interval plus the known second derivative there, from Taylor's theorem: f′(0) ≈ (f₁ − f₀)/h − (h/2) f″(0), with f″ = m/EI for w and f″ = q for m (mirrored at x = L). That is second-order accurate and touches only the last interval, so a point load one node from the end does not spoil it. (A three-point one-sided difference such as (−3f₀ + 4f₁ − f₂)/(2h) is also second order for smooth f, but it straddles the kink a point load puts in m: the first version of this solver used it and gave a cantilever tip deflection 50 % too large for a load one node from the tip. A test with a load near the end found it.) A point load P at an interior node enters as P/h in that node's m″ row, which makes m exact at the nodes; at a free end it enters the shear condition m′ = −P. Then `np.linalg.solve` (a banded solver would be faster; `P03.08`).

```python
for i in range(1, N - 1):
    A[row, Mc(i-1)], A[row, Mc(i)], A[row, Mc(i+1)] = 1/h**2, -2/h**2, 1/h**2;  b[row] = q[i];  row += 1
    A[row, W(i-1)],  A[row, W(i)],  A[row, W(i+1)]  = 1/h**2, -2/h**2, 1/h**2;  A[row, Mc(i)] = -1/EI[i];  row += 1
```

### Step 5: the elliptic closed form

`elliptic_tip(L_semi, s, EI)` is (5).

Running `main()` prints

```
spar EI = 6.958e+05 N m^2 (7075-T6, I = 9.704e+06 mm^4), semi-span 5.0 m, lift 13974.5 N
uniform lift: tip 0.3138 m; all lift at the elliptic centre of lift 4s/(3 pi): tip 0.1941 m
elliptic lift: root moment 29654.7 N m, tip slope 3.596 deg, tip deflection 0.2428 m by integration (4.86 % of semi-span)
closed form (superposed point loads): tip slope 3.596 deg, tip 0.24279 m; integration error -1.2e-06 m
finite differences, uniform EI: |FD tip - closed form| (m)
     n    uniform lift    elliptic lift
    25    5.021e-04       2.202e-03
    50    1.255e-04       8.114e-04
   100    3.138e-05       2.948e-04
   200    7.846e-06       1.062e-04
   400    1.961e-06       3.803e-05
  observed orders, uniform:  [2. 2. 2. 2.]
  observed orders, elliptic: [1.44 1.46 1.47 1.48] (load ~ sqrt(s - y) at the tip: order tends to 1.5)

superposition at n = 3.8: lift 0.2428 + wing weight -0.0377 (45 kg uniform) + tank -0.0186 (40 kg at 2.0 m) = 0.1866 m; FD with all loads 0.1865 m
root moment relieved from 29655 to 22481 N m by superposition (FD 22472 N m)

tapered spar (EI tip = 0.3 EI root): tip 0.2840 m by integration, 0.2840 m by FD

propped cantilever L = 1.2 m, q = 2000 N/m: prop force 3qL/8 = 900.0 N, clamp moment qL^2/8 = 360.0 N m (FD 360.0), max deflection 1.123 mm (FD 1.123 mm)
```

and writes `wing_deflection.png`: the deflected spar under elliptic lift by integration and by finite differences, the uniform-lift closed form and the tapered spar; and the log–log convergence plot with slopes 2 and 1.5.

### Verify

Run `python -m pytest curriculum/phase-05/lesson-12-beam-deflection`. The tests:

1. **Closed forms**: EI w⁗ = q by a five-point fourth difference; w = w′ = 0 at the clamp and EI w″ = EI w‴ = 0 at the free end; qL⁴/8EI, 5qL⁴/384EI, PL³/3EI, PL³/48EI; continuous deflection and slope under a point load; a straight line beyond a cantilever's load; the end-moment parabola.
2. **Reciprocity**: w(b) due to a load at a equals w(a) due to the same load at b, for both supports.
3. **Superposition and indeterminacy**: R = 3qL/8; the propped deflection vanishes at both ends and equals (4); finite differences reproduce it with a clamp moment qL²/8.
4. **Moment and integration**: `bending_moment` against statics for three cases; numerical integration against the cantilever and simply supported closed forms; the elliptic closed form (5) against integration to 10⁻⁴ and bracketed by the uniform and point-load cases.
5. **Finite differences**: observed order 1.8–2.2 against closed forms for clamped–free, pinned–pinned, clamped–clamped and clamped–pinned beams; a point load at a node reproducing the simply supported closed form with the exact nodal moment; a cantilever point load one node from the tip and at the tip against the closed form, with the root moment exact; superposition of two loads to round-off; zero load gives zero; doubling EI halves w and leaves m unchanged; a tapered spar under a smooth load matching integration at second order; the elliptic lift converging with order between 1.2 and 1.8.
6. **main** writes the figure.

### Compare

Three independent routes agree. The closed form (5) and the numerical integration differ by 1.2 × 10⁻⁶ m (trapezoid rule on 4001 points). Finite differences converge to the closed forms at order 2.00 under a uniform load, as the central difference's truncation error, (h²/12) times the fourth derivative of the solution, predicts. Under the elliptic lift the order is 1.44 at coarse grids, rising to 1.48 at n = 400 and 1.49 at n = 3200: the lift falls to zero at the tip like √(s − y), whose derivatives are infinite there, so the truncation-error estimate fails in the last cells and the order settles at 1.5. The error is still small (3.8 × 10⁻⁵ m at n = 400, 0.02 % of the deflection); it is the *rate* that a smooth-load test would have misjudged. That is the reason to test a numerical method on more than one load.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction before running it. The calls are one-liners on your module.

1. **Convergence.** `fd_beam` against the closed forms for uniform and elliptic lift, n = 25–400 (and up to 3200 for the elliptic case). Predict both orders. *Finding:* the table above; uniform 2.00; elliptic 1.44 → 1.49 (n up to 3200), tending to 1.5.
2. **Shape of the lift.** The same 13 974.5 N spread uniformly, elliptically, and triangularly (largest at the root, zero at the tip). Predict the ranking of root moment and tip deflection. *Finding:* uniform 34 936 N·m and 313.8 mm (4.80° tip slope); elliptic 29 655 N·m and 242.8 mm (3.60°); triangular 23 291 N·m and 167.4 mm (2.40°; the closed form w₀s⁴/(30EI) gives 167.37 mm). The tip deflection spreads by a factor of 1.9 for one total lift: the lift distribution has to be known before the deflection is worth computing to three figures.
3. **Taper.** EI falling linearly from the root value to (1 − taper) × root at the tip, taper 0–0.9, elliptic lift, integration against finite differences (n = 400). Predict whether taper 0.9 doubles the deflection. *Finding:* 242.8, 257.7, 269.7, 284.0, 302.2 mm for taper 0, 0.3, 0.5, 0.7, 0.9 (ratios 1.00, 1.06, 1.11, 1.17, 1.24), finite differences within 0.05 mm. Even a tip stiffness one tenth of the root's adds only 24 %, because the outboard wing has small moments to bend it: material near the tip contributes little stiffness and less strength, which is why spars taper.
4. **Where to put a mass.** A 40 kg tank at 3.8 g placed at y = 0–5 m: tip deflection and root moment relief. Predict how the tip relief scales with position. *Finding:* 0, −5.0, −18.6, −38.6, −62.8, −89.3 mm and 0, −1491, −2981, −4472, −5962, −7453 N·m at y = 0, 1, 2, 3, 4, 5 m. The moment relief is linear in y and the tip relief grows as y²(3s − y): outboard masses relieve far more (and they also lower the wing's natural frequencies, which is the flutter engineer's concern).
5. **Root-fitting stiffness.** A root that rotates by M₀/k_θ under the root moment adds M₀s/k_θ to the tip. For k_θ = 10⁶ to 10⁸ N·m/rad, compute the root rotation and the extra tip deflection. Predict the stiffness at which the root adds 10 %. *Finding:* 1.70°, 148.3 mm (61 %) at 10⁶; 0.57°, 49.4 mm (20 %) at 3 × 10⁶; 0.17°, 14.8 mm (6.1 %) at 10⁷; 0.6 % at 10⁸; the root adds 10 % at k_θ = 6.1 × 10⁶ N·m/rad, a root rotation of only 0.28°. Small rotations at the root are large deflections at the tip: measure the root when you measure the tip.

Raw data format: CSV tables `n, err_uniform_m, err_elliptic_m` (1); `shape, M_root_Nm, w_tip_mm, theta_tip_deg` (2); `taper, w_int_mm, w_fd_mm` (3); `y_tank_m, dw_tip_mm, dM_root_Nm` (4); `k_theta, root_rot_deg, extra_tip_mm` (5); plus `wing_deflection.png`.

## 5. Validate and improve

| quantity | theory (§2) | code | difference | cause |
|---|---|---|---|---|
| elliptic tip deflection | (5): 0.24279 m | integration 0.24279 m; FD (n = 400) 0.24275 m | 1.2 × 10⁻⁶ m; 3.8 × 10⁻⁵ m | trapezoid; FD truncation, order → 1.5 |
| uniform-lift tip | qL⁴/8EI = 0.313824 m | FD (n = 200) 0.313832 m | 7.8 × 10⁻⁶ m | FD truncation, order 2 |
| superposition, three loads | sum 0.1866 m, root moment 22 481 N·m | FD all together 0.1865 m, 22 472 N·m | 10⁻⁴ m, 4 × 10⁻⁴ | FD at n = 200 (elliptic tip) |
| tapered spar | integration 0.28404 m | FD 0.28399 m | 5 × 10⁻⁵ m | FD at n = 400 |
| propped cantilever | R = 3qL/8 = 900 N, M(0) = qL²/8 = 360 N·m | FD 360.0 N·m, max w 1.123 mm both | < 10⁻³ | — |

Every result is reproduced by an independent route with the discrepancy the method's error predicts. The engineering answer: **at limit load the spar alone lets the tip rise 243 mm (4.9 % of the semi-span) with a 3.6° tip slope under the elliptic lift; the wing's weight and a 40 kg tank at 2 m bring it to 187 mm; tapering the stiffness to 30 % at the tip adds 17 %.** Whether that is acceptable depends on what it does to the lift distribution and the controls, but the uncertainty is dominated by the lift shape (Experiment 2: a factor of 1.9 between plausible shapes) and the root stiffness (Experiment 5), not by the mathematics.

What to improve: the root fitting stiffness and the skins' contribution to EI, which together could change the answer by 20–50 %; then the lift distribution itself. What not to improve: the finite-difference grid beyond n = 200, or the small-slope approximation at 3.6°.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (3), (4) and (5); `wing_deflection.png`; the five experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) stating the tip deflection at limit load, the effect of relief masses and taper, and the assumption that most changes it.

### Adversarial self-checks

- **What if the load doubles?** The deflection doubles (linearity), to 486 mm and a 7.2° tip slope; w′² reaches 0.016, so the small-slope approximation is still good to a couple of per cent, but the stress (550 MPa, `P05.10`) is past yield, so the linear deflection is no longer valid either. Stiffness and strength fail together here.
- **What if EI goes to infinity?** w → 0 everywhere; the moments of a determinate beam do not change (the tests check that doubling EI halves w and leaves m alone), but those of an indeterminate one would, if EI changed in only part of it.
- **Which assumption fails first in the real world?** The rigid root: a few tenths of a degree of root rotation add 10 % or more. Then the spar-only section (skins stiffen the wing), then the fixed lift distribution (aeroelastic washout reduces outboard lift as the wing bends and twists).
- **What measurement could falsify the model?** A ground load test with the wing loaded in a known distribution and the deflection measured at several stations: the shape w(y) should follow (3)-like curves with the measured root rotation subtracted. A deflection shape with a kink means a joint or a crack; one proportional to y means the root, not the spar, is flexing.
- **Does the result survive a change of method or grid?** Closed form, integration and finite differences agree to 5 × 10⁻⁵ m; halving the grid changes the FD answer by the predicted factor.
- **Which decision changes if the uncertainty doubles?** If EI is uncertain by ±20 % instead of ±10 %, the predicted tip deflection spans 202–303 mm instead of 221–270 mm; a clearance or control-run requirement set within that band would need a stiffness test.

## Check yourself

1. Write the relation between bending moment and deflection, and say what each symbol means and which sign convention it uses.
2. State the two boundary conditions at a clamped end, a pinned end and a free end.
3. A cantilever's tip deflection under a uniform load is qL⁴/(8EI). By what factor does it change if the length doubles? If the depth of a rectangular section doubles?
4. Why does superposition hold for beam deflections, and when does it stop holding?
5. How does superposition give an extra equation for a propped cantilever, and which of equilibrium, compatibility and the constitutive law supplies it?
6. State Maxwell's reciprocal theorem in one sentence.
7. Why did the finite-difference solution converge at order 2 for a uniform load but more slowly for an elliptic one?
8. A wing-mounted tank is moved outboard. What happens to the root bending moment and to the tip deflection at the same load factor?
