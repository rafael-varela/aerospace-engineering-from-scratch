# P05.09  Shear & Bending Moment Diagrams

> Cut the beam anywhere and ask what holds the left piece still. The answer, at every x, is the shear and bending moment diagram; everything a beam does later in this phase is read off it.

| | |
|---|---|
| **Learning outcomes** | 1. Find the reactions of statically determinate beams (simply supported, overhanging, cantilevered) and recognise when the supports make a beam indeterminate or a mechanism. 2. Derive the relations dV/dx = w and dM/dx = V and the jumps at point forces and couples from equilibrium of a slice, and write V(x) and M(x) for point forces, couples and linearly varying distributed loads with singularity (Macaulay) functions. 3. Implement and test a V–M diagram generator, verified against an independent method of sections with quadrature, the classic closed forms, and the extreme moments found from the zeros of V. 4. Compute the shear and bending moment along a wing spar under an elliptical lift distribution with inertia relief from structure, fuel and an engine, show how fast a panelled load converges, and use the generator to rank design choices that reduce the root bending moment. |
| **Prerequisites** | `P05.01` and `P05.02` (free-body diagrams, support reactions, determinacy), `P05.05` (resultant and centroid of a distributed load), `P01.03` (integrals as areas), `P03.03` (piecewise-linear interpolation), `P03.04` (quadrature), `P01.13` (2 × 2 linear systems). Skill check: you can find the reactions of a simply supported beam with one off-centre point load by taking moments about one support. |
| **Engineering question** | A 1900 kg light twin (numbers illustrative of the class) pulls up at a limit load factor of 3.8. Each wing (semispan 6.7 m, measured from the aircraft centre-line) carries half the lift, distributed elliptically, and the inertia of its own structure (80 kg), fuel (75 kg between 0.5 and 3.0 m) and engine (150 kg at 1.9 m), each n times its weight, pushes down. What are the shear force and bending moment along the spar, how much does inertia relief reduce the root bending moment, and how would moving the engine or the fuel change it? The answers: lift alone would give a root shear of 35.4 kN and a root bending moment of 100.7 kN·m; inertia relief reduces them to 24.0 kN and 75.2 kN·m (to 74.7 %), with the engine contributing the largest single relief, 10.6 kN·m; moving the engine from 1.2 to 3.3 m lowers the root moment by 11.7 kN·m, and putting the fuel outboard lowers it by another 7.0 kN·m. |
| **Assumptions and limits** | Straight, slender beams loaded transversely (no axial force, no torsion); loads and supports in one plane; statically determinate supports; small deflections, so equilibrium is written on the undeformed beam; the lift as an ideal elliptical distribution over the semispan with the root at the centre-line (the fuselage carry-through structure is treated as part of the spar); inertia loads as uniform or point loads at the stated positions; a symmetric pull-up (no roll, no gust, no aileron); concentrated loads as true points (a real fitting spreads them over its footprint). |
| **Analytical work** | Reactions from ΣF = 0 and ΣM = 0, and the count of unknowns; the internal V and M by the method of sections; the differential relations and jump conditions from a slice; singularity functions and their integration; the resultants and centroids of uniform, triangular and trapezoidal loads; closed forms wL²/8, Pab/L, PL, wL²/2 and w₀L²/(9√3); the root shear πw₀s/4 and moment w₀s²/3 of an elliptical lift distribution and its centre of lift 4s/3π; the convergence order of a panelled load with a square-root tip. |
| **Computational deliverable** | `starter/main.py`: `macaulay`, `load_resultant`, `solve_reactions`, `shear_moment`, `sections_check`, `piecewise_linear_load`, `extreme_moments` (the load dataclasses, `breakpoints`, the `WingCase` and floor-beam data and `main` are given). Tests prove: reactions satisfy ΣF = 0 and ΣM = 0 about any point on random beams; indeterminate and unstable supports are refused; V and M vanish beyond the end of the beam; the singularity-function diagram equals the method of sections with quadrature at random cuts to 10⁻⁶; V jumps by P at a force and M by −M₀ at a couple; dM/dx = V and dV/dx = w by finite differences; wL²/8, Pab/L and w₀L²/(9√3) come out at the right places; cantilever root values; symmetry and superposition; a vanishing load gives a vanishing diagram; the panelled elliptical lift converges to its closed-form root moment at order about 1.5; inertia relief, the engine jump and an outboard engine behave as derived. |
| **Practical activity** | **computational experiment**: a convergence study of the panelled elliptical lift against its closed form, a method comparison between singularity functions and sections, and a design sweep of engine position and fuel placement against root bending moment. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `vm_diagrams.png`, the derivations of the slice relations, the singularity-function terms for a linearly varying load and the elliptical root loads, the experiment tables, and a half-page note on the wing's root loads with the inertia-relief breakdown and a recommendation on engine and fuel placement. |
| **Assessment** | `assessments/phase-05/grading-rubric.md`; competencies 1–4. Pass: the relations dV/dx = w, dM/dx = V and the jump conditions are derived from a free-body diagram of a slice; all public tests pass; the sign convention is stated and used consistently; the wing note gives the root loads with and without relief and the convergence evidence. |
| **References** | `@libretexts-baker-haynes-engineering-statics` Ch. 8 "Internal Loadings": 8.1 Internal Loads; 8.2 Sign Conventions; 8.3 Internal Loads at a Point; 8.4 Shear and Bending Moment Diagrams; 8.5 Section Cut Method; 8.6 Relation Between Loading, Shear and Moment; 8.7 Graphical Method; 8.8 Integration Method; and 7.8 Distributed Loads (required-free, https://eng.libretexts.org/Bookshelves/Mechanical_Engineering/Engineering_Statics%3A_Open_and_Interactive_(Baker_and_Haynes)/08%3A_Internal_Loadings): the method of sections, the sign conventions and the graphical and integration methods. Note that Baker & Haynes take distributed loads positive downward, so their relation reads dV/dx = −w; this lesson takes w positive up and writes dV/dx = w. `@mit-ocw-3-11-mechanics-of-materials` Roylance module "Shear and Bending Moment Diagrams" (required-free, https://ocw.mit.edu/courses/3-11-mechanics-of-materials-fall-1999/resources/mit3_11f99_statics/): a second, compact treatment of the same diagrams, useful to compare conventions. |
| **Extensions** | ★ Add a second spar station at the wing strut of a strut-braced high wing and show that the beam becomes indeterminate; solve it in `P05.12` once deflections are available. ★ Replace the elliptical lift by the average of an elliptical and a planform-shaped (trapezoidal) distribution, a classic first estimate for tapered wings, and compare the root moment. ★ Add a rolling pull-up: an extra antisymmetric lift distribution, and the V–M diagrams of the two wings. |
| **Time** | unmeasured |

## 1. The engineering problem

A wing is a cantilever: bolted to the fuselage at the root, free at the tip, loaded all along its length by lift. Before anyone can size the spar caps (`P05.10`), the webs (`P05.11`) or check the tip deflection (`P05.12`), they need two curves along the span: the **shear force** V(x), which the spar webs carry, and the **bending moment** M(x), which the spar caps carry as tension and compression. Their values at the root, where they are largest, size the heaviest and most critical part of the airframe, the wing-to-fuselage attachment.

Two things make the wing case more than a textbook exercise. The lift is not uniform: it is distributed roughly elliptically, highest at the root and falling to zero at the tip with a square-root shape, which is not a polynomial and must be approximated or integrated numerically. And the wing is not weightless: in a 3.8 g pull-up its own structure, its fuel and its engine are all pulled down with 3.8 times their weight, and those downward inertia loads *relieve* the bending caused by lift. Where the mass sits along the span matters: mass far out on the wing relieves more. That is why airliners carry their fuel in the wings, why engines hang under them, and why a maximum zero-fuel weight is a structural limit.

This lesson builds the diagram generator. Part 2 derives the relations between load, shear and moment from a free-body diagram of a slice, then writes the diagram for any combination of point forces, couples and linearly varying loads with singularity functions. Part 3 implements it with a method-of-sections check that shares none of its algebra. Part 4 applies it to the wing and to a fuselage floor beam. Success criteria, written first: the reactions balance the loads to round-off; V and M close to zero beyond the end of every beam; the singularity-function diagram equals the method of sections at random cuts on random beams; the classic closed forms appear at the right places; and the panelled elliptical lift converges to its exact root moment at a measured rate.

The organising idea of the phase is *equilibrium, compatibility and constitutive law*. This lesson is pure equilibrium: the shear and moment diagrams of a determinate beam do not depend on the material or the cross-section. That is also their limit: an indeterminate beam needs compatibility, which waits for deflections in `P05.12`.

## 2. Understand and derive

### 2.1 Sign convention

Put x along the beam from its left end. Take **forces positive up**, **distributed loads w(x) positive up** (in N/m), and **applied couples positive counterclockwise**. Cut the beam at x and keep the part to the left. Define

- the **shear force V(x)** as the sum of the upward forces acting on the left part (so the right part pushes down on it with V), and
- the **bending moment M(x)** as positive when it bends the beam concave up ("sagging", top in compression).

With these choices M(x) is the sum of the moments, about the cut, of the forces on the left part taken clockwise positive, minus any counterclockwise couples on the left part. Sign conventions differ between books (Baker & Haynes 8.2 takes downward loads positive); any convention works if it is stated and kept.

![A beam cut at distance x from its left end into a left part and a right part. The left part carries an upward distributed load w(x) and the reaction R_A at its end; on its cut face the shear V acts downward and the moment M counterclockwise, and on the right part's face V acts upward and M clockwise. Below, a positive moment bends a beam concave up (sagging) and a negative one concave down (hogging).](figures/shear-moment-sign-convention.svg)
*Figure 1. The sign convention of this lesson. V and M act in opposite senses on the two faces of the same cut (Newton's third law); a positive M is the one that makes the beam sag.*

### 2.2 Reactions and determinacy

A beam loaded transversely in a plane has two independent equilibrium equations, ΣF = 0 and ΣM = 0 about any point (the horizontal equation is trivially satisfied when there are no axial loads). A pin and a roller give two unknown vertical forces; a fixed end gives a vertical force and a couple. Either way there are two unknowns, and the beam is **statically determinate**: the 2 × 2 system has a unique solution. Three or more unknowns (a third support, or a fixed end with an extra roller) make it **indeterminate**; one unknown, or two forces at the same point (which cannot resist a couple), make it a **mechanism**. The generator refuses both (`P05.02`).

A distributed load enters the equilibrium equations through its resultant and the position of its centroid (Baker & Haynes 7.8). For a load varying linearly from w_a at x = a to w_b at x = b (length ℓ = b − a), the resultant is (w_a + w_b)ℓ/2 and it acts at a + ℓ(w_a + 2w_b)/(3(w_a + w_b)) from the left end; for a uniform load that is the midpoint, for a triangle rising to w_b it is two-thirds of the way along.

### 2.3 The slice: dV/dx = w, dM/dx = V, and the jumps

Take a slice of the beam between x and x + dx carrying a distributed load w. On its left face act V and M (as defined for the left part, so on the slice's left face the left part pushes up with V); on its right face act −(V + dV) and the moment M + dM in the opposite sense. Vertical equilibrium: V + w dx − (V + dV) = 0, so

**dV/dx = w.**   (1)

Moments about the right face: M + V dx + w dx²/2 − (M + dM) = 0, and dropping the second-order term,

**dM/dx = V.**   (2)

Integrate (1) and (2) across a point force P at x = a (a slice of vanishing width): V jumps by **+P**, M is continuous. Across a counterclockwise couple M₀: V is continuous and M jumps by **−M₀**. Three consequences carry the whole of hand-sketched diagrams (Baker & Haynes 8.6–8.7): the slope of V is the load, the slope of M is the shear, and **M has its extreme values where V = 0** or changes sign, or at a couple, or at an end.

### 2.4 Singularity functions

The **Macaulay bracket** is

**⟨x − a⟩ⁿ = (x − a)ⁿ for x ≥ a, and 0 for x < a  (n ≥ 0),**

so ⟨x − a⟩⁰ is a unit step at a and ⟨x − a⟩¹ a ramp that starts there. It integrates like an ordinary power: ∫⟨x − a⟩ⁿ dx = ⟨x − a⟩ⁿ⁺¹/(n + 1). (The cases n = −1 and n = −2, the Dirac delta and the doublet, describe a point force and a point couple in the load function itself; this lesson starts one integration later and never needs to evaluate them.) Integrating (1) and (2) from the left end, where V = M = 0 just before x = 0, each load contributes one closed-form term:

| load | V(x) | M(x) |
|---|---|---|
| point force P at a | P⟨x − a⟩⁰ | P⟨x − a⟩¹ |
| couple M₀ (ccw) at a | 0 | −M₀⟨x − a⟩⁰ |
| linear load w_a → w_b on [a, b], k = (w_b − w_a)/(b − a) | w_a⟨x − a⟩¹ + (k/2)⟨x − a⟩² − w_b⟨x − b⟩¹ − (k/2)⟨x − b⟩² | (w_a/2)⟨x − a⟩² + (k/6)⟨x − a⟩³ − (w_b/2)⟨x − b⟩² − (k/6)⟨x − b⟩³ |

The distributed load is written as a load that starts at a and continues forever, minus the same load (with its value and slope at b) starting at b, so that the two cancel beyond b. **Reactions are loads too**: add them to the list, and the whole diagram is one sum. Beyond the right end, a correct set of reactions makes V and M vanish identically: the whole beam is in equilibrium. That closure is the first check on any diagram.

### 2.5 Worked closed forms

Simply supported span L, uniform load w (down, so w < 0 in our convention): reactions −wL/2 each; V = −wL/2 + wx, zero at L/2; **M_max = −wL²/8 at x = L/2** (500 N·m for 1000 N/m over 2 m). Point load P (down) at a, with b = L − a: **M_max = Pab/L under the load** (1125 N·m for 3000 N at 0.5 m on 2 m). A triangular load rising from 0 to w₀ over a simple span: V = 0 at L/√3, **M_max = w₀L²/(9√3)**. Cantilever of length L fixed at x = 0 with a tip load P (down): the root reaction pushes up with P, so V(0) = P, and **M(0) = −PL** (hogging): 3000 N on 2 m gives V = 3000 N and M = −6000 N·m. With a uniform load w as well, M(0) = −PL − wL²/2.

### 2.6 The elliptical wing load

Take the lift per unit span as **w(x) = w₀√(1 − (x/s)²)** on 0 ≤ x ≤ s, with x measured from the root (the aircraft centre-line) and s the semispan. Its integral is a quarter of an ellipse's area, **πw₀s/4**, which must equal half the lift, nmg/2:

w₀ = (nmg/2)/(πs/4) = 6728 N/m for the light twin (m = 1900 kg, n = 3.8, s = 6.7 m).

The root shear is the whole lift on the half-wing, **35 402 N**. The root moment is ∫₀ˢ w x dx; substituting u = 1 − (x/s)², it is **w₀s²/3 = 100 668 N·m**. Their ratio puts the centre of lift at **4s/(3π) = 2.8436 m** from the root, the centroid of a quarter-ellipse.

The elliptical shape is not a polynomial, so the singularity functions of 2.4 cannot take it directly. `piecewise_linear_load` replaces it by n linear panels through its values at n + 1 equally spaced points. Because w has a square-root behaviour at the tip (its slope is infinite there), the panels converge more slowly than the second order you would expect of linear interpolation of a smooth function: the error in the root moment falls as n^(−1.5), as part 4 measures. The method of sections, integrating w(x) itself with `scipy.integrate.quad`, needs no panels at all, which makes it the independent check.

**Inertia relief.** At load factor n, every mass on the wing feels n times its weight downward. Uniform wing structure of 80 kg over 6.7 m gives −2981 N; 75 kg of fuel spread over 0.5–3.0 m gives −2795 N; a 150 kg engine at 1.9 m gives a point load of −5590 N. Their moments about the root are −9987, −4891 and −10 621 N·m: the engine, being heavy and well outboard, relieves the most. Net root shear: 35 402 − 11 366 = 24 036 N; net root moment 100 668 − 25 499 = 75 169 N·m (computed exactly by the method of sections; the 400-panel singularity-function solution gives 75 158 N·m).

With the lesson's sign convention, every cut on the wing has the root reaction (which pushes down, against the lift) to its left, so V(x) is negative along the span while M(x) is positive (sagging: lift bends the wing up, compressing the upper spar cap). The magnitudes are what the spar carries.

### 2.7 A fuselage floor beam

A floor cross-beam spans 3.6 m between two fuselage frames, pinned at one end and on a roller at the other (an idealisation). Four seat-track loads of 2.0 kN act downward at 0.45, 1.25, 2.35 and 3.15 m, the floor panels add 0.5 kN/m, and an equipment fitting at 2.0 m applies a counterclockwise couple of 0.6 kN·m (illustrative values). Equilibrium gives reactions of 5066.7 N and 4733.3 N. The moment reaches **4533.3 N·m just to the left of the couple** and drops to 3933.3 N·m just to the right of it: the maximum is at the couple, not where V = 0, which is why the extreme-moment search must check both sides of every couple.

### 2.8 Where the model stops

- The diagram of a determinate beam depends only on equilibrium; it is exact for the stated loads and supports. The errors are in the loads: the real lift distribution of a tapered, twisted wing differs from an ellipse (the aerodynamic lessons of `P10` compute it), and gust and manoeuvre loads are set by the flight envelope.
- The wing root is idealised as a fixed end at the centre-line; in reality the spar continues through a carry-through box and the fuselage reacts the loads through frames.
- Point loads are idealisations: an engine pylon spreads its load over two or three fittings, which rounds the shear jump.
- Small deflections: a very flexible wing (a sailplane bends several per cent of its span) changes the geometry on which equilibrium is written.

## 3. Build the computational model

Open `starter/main.py`. The load dataclasses (`PointForce`, `PointMoment`, `DistributedLoad`, `DistributedFunction`), `breakpoints`, the `WingCase` (with `w0`, `lift`, `inertia_loads`, `loads` and `exact_lift_root`), the floor-beam data and `main` are given.

### Step 1: singularity functions and resultants

`macaulay(x, a, n, left=False)` implements the bracket; with `left=True` the step is taken as 0 at x = a, which gives one-sided limits from the left. `load_resultant(load)` returns the force and its counterclockwise moment about x = 0: P and aP for a point force, 0 and M₀ for a couple, the trapezoid integrals for a linear load, and two `quad` integrals for a `DistributedFunction`.

### Step 2: reactions

`solve_reactions(supports, loads)` counts the unknowns (one per pin or roller, two per fixed end), refuses anything but two, assembles

```python
# columns: the unknowns; rows: sum F = 0 and sum M about x = 0 = 0
A = [[1 if force else 0 ...], [x if force else 1 ...]]
r = np.linalg.solve(A, -np.array([F_applied, M_applied]))
```

and returns the reactions as load objects.

### Step 3: the diagram

`shear_moment(x, loads, left=False)` sums the table of 2.4 over the loads (reactions included). `sections_check(x, loads)` is the independent route: for one cut it adds the upward forces to the left and their moments about the cut, integrating every distributed load (linear or functional) over [a, min(x, b)] with `quad`, with no Macaulay brackets anywhere. `piecewise_linear_load` builds the panels. `extreme_moments(loads, L)` evaluates M on both sides of every breakpoint, then in each interval between breakpoints recovers V (a polynomial of degree at most 2 there) exactly from three samples, finds its real roots in the interval with `np.roots`, and evaluates M at them; the largest and smallest of all candidates are the extremes.

Running `main()` prints the closed-form checks of 2.5, the floor beam of 2.7 with its sections check (difference 0.0), and the wing:

```
  panels   root V (N)   root M (N m)   rel. error in M (lift only)
      10    -34984.2       98155.1     2.50e-02
      20    -35254.1       99745.9     9.16e-03
      40    -35349.7      100334.1     3.32e-03
      80    -35383.5      100548.1     1.19e-03
     160    -35395.5      100625.2     4.26e-04
     320    -35399.7      100652.8     1.52e-04
  observed orders [1.45 1.47 1.48 1.48 1.49]
```

followed by the relief breakdown, the singularity-against-sections table along the span and the design sweep of part 4, and writes `vm_diagrams.png` (V and M of the floor beam; V and M of the wing with and without relief).

### Verify

Run `python -m pytest curriculum/phase-05/lesson-09-shear-bending-moment-diagrams`. The tests:

1. **Brackets**: steps, ramps and parabolas in the right places; right-continuity; zero delta away from a.
2. **Resultants**: a trapezoidal load's force and centroid; the same load as a function through `quad`.
3. **Reactions**: ΣF = 0 and ΣM = 0 about the origin and about a random point, for random beams on four support layouts.
4. **Determinacy**: three supports, fixed plus roller, one roller and two forces at one point all refused.
5. **Closure**: V = M = 0 beyond the right end; nothing before x = 0.
6. **Sections**: singularity functions against the method of sections at 12 random cuts on random beams, to 10⁻⁶.
7. **Jumps**: ΔV = P at a force with M continuous; ΔM = −M₀ at a couple with V continuous.
8. **Derivatives**: central differences of M equal V, and of V equal w, inside a smooth interval.
9. **Simple spans**: wL²/8 at L/2, Pab/L at a, w₀L²/(9√3) at L/√3.
10. **Cantilever**: root shear and root moment, and the minimum at the root.
11. **Symmetry and superposition**: symmetric loads give symmetric M and antisymmetric V; diagrams of two load sets add.
12. **Vanishing load**: a zero load gives a zero diagram.
13. **Wing**: half the lift and the 4s/3π centre; panels converging at order above 1.3 to the closed-form root moment, below 3 × 10⁻⁴ at 320 panels; relief reducing the root shear by exactly the inertia loads; the engine jump; an outboard engine relieving more; sections with the exact ellipse agreeing to 10⁻³.
14. **Floor beam**: reactions balance; `extreme_moments` matches the maximum of a 20 001-point evaluation.
15. **main** writes the figure.

### Compare

The method of sections and the singularity functions are both equilibrium, but they share no code and no algebra: one integrates the load with adaptive quadrature at a single cut, the other evaluates closed-form brackets everywhere. Their agreement to 10⁻⁶ on random beams is the evidence that the bracket table of 2.4 is right. For the wing, `quad` on the exact ellipse is also the reference against which the panels converge. Every statics text's graphical method (Baker & Haynes 8.7) is a third route you can use by hand on the floor beam.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction first.

1. **Panel convergence.** The elliptical lift alone with 10 to 320 panels; root shear and moment against the closed forms. Predict the order (2 for a smooth load; what does the square-root tip do?).
2. **Singularity functions against sections.** The wing with relief, at x = 0, 1.0, 1.9, 3.0 and 5.0 m: V and M from 400 panels against sections with the exact ellipse. Predict the size of the difference from experiment 1.
3. **Design sweep.** Root bending moment with the engine at 1.2, 1.9, 2.6 and 3.3 m, and with the fuel at 0.5–3.0 m, at 3.0–5.5 m, or absent (the same aircraft mass, the payload in the fuselage instead). Predict the sensitivity to each.

Raw data format: CSV tables `panels, root_V_N, root_M_Nm, rel_error` (1); `x_m, V_sing, V_sections, M_sing, M_sections` (2); `engine_x_m, fuel_case, root_M_Nm` (3); plus `vm_diagrams.png`.

What you should find. (1) Root moment errors of 2.50 × 10⁻² at 10 panels down to 1.52 × 10⁻⁴ at 320, observed orders 1.45 to 1.49: **order 1.5**, the square-root tip, not 2. (2) Differences of about 1.7 N in V and 2.7–11 N·m in M (for example 75 158.3 against 75 169.3 N·m at the root), the size of the 400-panel error, and decreasing outboard. (3) Root moments:

| engine at (m) | fuel 0.5–3.0 m | fuel 3.0–5.5 m | no fuel |
|---|---|---|---|
| 1.2 | 79 071 | 72 084 | 83 962 |
| 1.9 | 75 158 | 68 171 | 80 049 |
| 2.6 | 71 245 | 64 258 | 76 137 |
| 3.3 | 67 333 | 60 345 | 72 224 |

Each 0.7 m of engine movement outboard removes 3913 N·m (the engine's inertia force, 5590 N, times 0.7 m, exactly, because the diagram is linear in the load positions); moving the fuel outboard by 2.5 m removes 6987 N·m; removing the fuel from the wing while keeping the aircraft mass adds 4891 N·m. The zero-fuel case is the worst: that is the structural meaning of a maximum zero-fuel weight.

## 5. Validate and improve

| quantity | closed form or independent route | generator | difference |
|---|---|---|---|
| simply supported UDL | wL²/8 = 500.0 N·m at 1.000 m | 500.0 N·m at 1.000 m | 0 |
| simply supported point load | Pab/L = 1125.0 N·m at 0.500 m | 1125.0 N·m at 0.500 m | 0 |
| cantilever tip load | V = P, M = −PL = −6000 N·m | 3000.0 N, −6000.0 N·m | 0 |
| floor beam, 4 cuts | method of sections | singularity functions | 0.0 |
| elliptical lift root moment | w₀s²/3 = 100 668 N·m | 100 652.8 N·m (320 panels) | 1.52 × 10⁻⁴, order 1.5 |
| wing with relief, root moment | sections with exact ellipse: 75 169.3 N·m | 75 158.3 N·m (400 panels) | 1.5 × 10⁻⁴ |
| shear jump at the engine | n m g = 5590 N | −5590 N | 0 |

Every check is met to round-off or to the panel error that experiment 1 measures. **The engineering answer**: at the 3.8 g limit load the light twin's wing carries a root shear of 24.0 kN and a root bending moment of 75.2 kN·m, 74.7 % of the 100.7 kN·m that lift alone would produce; the engine provides the largest relief (10.6 kN·m), the wing structure 10.0 kN·m and the fuel 4.9 kN·m. Moving the engine outboard and keeping fuel in outboard tanks both reduce the root moment, at the price of roll inertia and a larger load at the engine pylon; the critical case for the spar is the zero-fuel one. These are limit loads; the spar is designed to 1.5 times them at ultimate (`P05.14`).

What to improve: a realistic lift distribution for the actual planform (`P10`), the wing's own mass distributed in proportion to the chord rather than uniformly, and the pylon's two-point load introduction. What not to improve: the number of panels beyond a few hundred, which already puts the discretisation error two orders of magnitude below the uncertainty of the lift distribution itself.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; derivations of (1), (2), the jump conditions, the linear-load row of the table in 2.4, and the root loads of the elliptical distribution; `vm_diagrams.png`; the three experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) on the wing's root loads with the relief breakdown and a recommendation on engine and fuel placement.

### Adversarial self-checks

- **What if the load doubles?** Every term in V and M is linear in the loads, so the whole diagram doubles (superposition, tested). The load factor scales lift and inertia relief together, so the root moment at n = 7.6 is twice that at 3.8 if the mass distribution is unchanged.
- **What if the load goes to zero?** The diagram vanishes (tested). With lift but no relief the root moment is 100 668 N·m; with relief but no lift (the aircraft resting on its wheels at n = 1, gear in the fuselage) the wing hangs down and M changes sign.
- **Which assumption fails first in the real world?** The lift distribution: a tapered wing at high angle of attack loads its tip more than an ellipse, moving the centre of lift outboard and raising the root moment. Then the point-load idealisation of the engine.
- **What measurement could falsify the model?** Strain gauges on the spar caps near the root in a pull-up (or a ground test with whiffle-tree loading) give the bending moment through the flexure formula of `P05.10`; a measured root moment that disagrees with the integral of the applied loads beyond the gauge calibration error would mean the loads, not the statics, are wrong.
- **Does the result survive a change of method?** Singularity functions with 400 panels and sections with the exact ellipse differ by 11 N·m out of 75 kN·m at the root; the floor beam's two methods agree to 0.0.
- **Which decision changes if the uncertainty doubles?** If the lift distribution's centre of lift were uncertain by ±5 % instead of ±2.5 %, the lift-alone root moment would be uncertain by ±5 kN·m, comparable to the effect of moving the fuel outboard: the fuel-placement decision would need a better lift estimate first.

## Check yourself

1. State dV/dx and dM/dx in this lesson's sign convention. How do they change if distributed loads are taken positive downward?
2. A simply supported beam of span 4 m carries a uniform downward load of 2 kN/m. What are the reactions and the maximum bending moment, and where is it?
3. What happens to V and M across a point force P, and across a counterclockwise couple M₀?
4. Write ⟨x − 2⟩¹ at x = 1 and at x = 3. What is its integral?
5. Why must V and M be zero just beyond the right end of a beam whose reactions are correct?
6. A cantilever of length L carries a uniform downward load w. What are the root shear and root moment?
7. Why does a heavy engine mounted on the wing reduce the root bending moment in flight, and why does it increase it on the ground?
8. The elliptical lift panels converge at order 1.5 rather than 2. What feature of the load causes it?
