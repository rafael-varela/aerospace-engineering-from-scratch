# P05.01  Forces, Moments & Free-Body Diagrams

> Cut the part free, draw every force that touches it, and the structure tells you what it carries.

| | |
|---|---|
| **Learning outcomes** | 1. Represent forces as vectors, compute the moment of a force about a point (M = r × F) and about an axis (û · (r × F)), and show that the axis moment does not depend on the point chosen on the axis. 2. Replace any system of forces and couples by a force-couple at a chosen point, move it with M_B = M_A + (A − B) × R, and find its central axis (the line along which it is a single force, or a force plus a parallel couple). 3. Reduce a distributed load to its resultant and line of action, by hand for uniform, triangular and elliptic loads and numerically for any sampled load. 4. Draw the free-body diagram of a landing-gear leg with a trunnion pin and a drag brace, recognise the brace as a two-force member, and solve for the three unknowns both as a 3 × 3 linear system and by one moment equation about the trunnion. 5. Implement and test these tools in Python and use them to explore how the brace load depends on the brace geometry. |
| **Prerequisites** | `P04.01` (Newton's laws, free-body diagrams of particles), `P04.04` (torque as r × F), `P01.07` (dot and cross products), `P01.13` (linear systems). Skill check: you can compute a × b for two 3-vectors by hand and solve a 3 × 3 system with `numpy.linalg.solve`. |
| **Engineering question** | On a level touchdown, the ground pushes up on one main-gear tyre of a 1100 kg light aircraft with 14.6 kN and, while the wheel spins up, drags it aft with 3.6 kN. How much do the trunnion pin and the drag brace carry, is the brace in tension or compression, and how sensitive is its load to where the brace is attached? |
| **Assumptions and limits** | Rigid bodies; static equilibrium (the touchdown is treated as a quasi-static load case, not a dynamic impact); pins are frictionless and transmit no moment; the drag brace is pinned at both ends and carries no load between them; the leg's own weight and the wheel's inertia are neglected; planar (side-view) loads only for the gear. The load case is illustrative, not a certification case. |
| **Analytical work** | Moment of a force about a point and an axis; Varignon's theorem; couples; transport of a force-couple system; the central axis and its pitch; resultant and centroid of uniform, triangular and elliptic line loads; the free-body diagram of the gear leg; the brace force from ΣM_A = 0; the trunnion reactions from ΣF = 0. |
| **Computational deliverable** | `starter/main.py`: `moment_about_point`, `moment_2d`, `moment_about_axis`, `resultant`, `transport_moment`, `central_axis`, `distributed_load_resultant`, `elliptic_lift`, `equilibrium_residual_2d`, `solve_fbd_2d`, `gear_reactions`, `brace_force_by_moments`, `brace_lever_arm`, `gear_fbd_forces`, `landing_load_case`, `brace_sweep` (`unit`, `to3`, `force_vector`, `ground_force`, `GearGeometry` and `main` are given). Tests prove: moments are perpendicular to lever and force; Varignon's theorem; the axis moment is independent of the point on the axis; transport and couple invariance; the central axis carries a moment parallel to the resultant; distributed-load closed forms and convergence; zero residual about any point; refusal of concurrent or parallel unknowns; the 3 × 3 solve equals the single moment equation; superposition and the zero-drag limit. |
| **Practical activity** | **computational experiment**: a sensitivity study of the drag-brace load against its attachment height, the drag-to-vertical load ratio and ±5 mm of geometric tolerance. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `landing_gear_fbd.png`, a hand-drawn free-body diagram of the leg with every force labelled, the hand calculation of the brace force, the three experiment tables, and a half-page note answering the engineering question. |
| **Assessment** | `assessments/phase-05/grading-rubric.md`; competencies 1–4. Pass: the free-body diagram is complete and drawn before any equation; the brace force is found by hand from one moment equation and matches the code; all public tests pass; the sensitivity results are explained, not just listed. |
| **References** | `@libretexts-baker-haynes-engineering-statics` Ch. 2 "Forces and Other Vectors" (2.5 Unit Vectors, 2.7 Dot Products, 2.8 Cross Products), Ch. 4 "Moments and Static Equivalence" (4.1 Direction of a Moment, 4.2 Magnitude of a Moment, 4.3 Varignon's Theorem, 4.4 Moments in Three Dimensions, 4.5 Couples, 4.6 Equivalent Transformations, 4.7 Statically Equivalent Systems) and Ch. 5 "Rigid Body Equilibrium" 5.2 Free Body Diagrams (required-free: the vector toolkit and the free-body habit); `@openstax-university-physics-vol1` Ch. 12, 12.1 Conditions for Static Equilibrium (required-free: the Phase 04 free-body diagram, now with torque); `@tudelft-alderliesten-aerospace-structures-materials` Ch. 6 "Aircraft & spacecraft loads": 6.2 Externally loaded airframe, 6.3 Load path, 6.4 Loads and load paths in an airframe (required-free, link only: how loads enter and travel through an airframe). The landing-gear load cases and the solver are original to this lesson. |
| **Extensions** | ★ Add a side load and a side brace and solve the gear in 3-D (lesson `P05.02` gives the tools). ★ Make the leg raked (contact point ahead of the trunnion) and find the drag ratio at which the brace carries nothing. ★ Replace the trapezoidal rule in `distributed_load_resultant` by a substitution y = s sin u that removes the square-root tip of the elliptic load and restores fast convergence. |
| **Time** | unmeasured |

## 1. The engineering problem

Every structural calculation in this phase starts the same way: pick a part, cut it free from everything it touches, and replace each cut by the force or moment that the removed piece exerted. That drawing, the **free-body diagram**, is where most structural mistakes are made, and where most of them can be caught. A missing force makes the equations unsolvable or, worse, solvable with a wrong answer.

The part in this lesson is the main landing-gear leg of a light aircraft. It hangs from the airframe on a **trunnion**, a pin at its top about which it can swing; on a retractable gear this is the pin it folds about. A **drag brace**, a strut pinned to the leg partway down and to the airframe ahead of it, stops the leg from swinging aft. On touchdown the runway pushes up on the tyre and, while the stationary wheel is spun up to the landing speed, friction drags the tyre aft. The design question is simple to state: what forces do the trunnion pin and the drag brace carry? A brace sized for the wrong force, or in the wrong sense (a strut designed as a tie will buckle when it is pushed), is a classic failure.

The lesson builds the tools to answer it: forces and moments as vectors, moments about axes, equivalent force systems, distributed loads, and a small free-body-diagram solver. For a 1100 kg aircraft in an illustrative level-landing case (14 563 N up and 3 641 N aft on one main tyre), the brace turns out to be a **tie in tension at 10 501 N**, nearly three times the drag load, and the trunnion pin carries **21 015 N**, more than the vertical ground load itself. Success criteria, written first: the brace force from the 3 × 3 linear system equal to the one found from a single moment equation to round-off; zero force and moment residual about any point; every number in this paragraph reproduced by the code.

## 2. Understand and derive

### 2.1 Forces as vectors

A force has a magnitude, a direction and a point of application. In components, **F** = (F_x, F_y, F_z) in newtons; a force of magnitude F along a direction **d** is F **d**/|**d**|, the direction made into a unit vector (Baker & Haynes 2.5). For a **rigid body**, one whose deformation we ignore, the point of application can slide along the force's **line of action** without changing its effect; it cannot be moved off that line. Two forces with the same magnitude and direction on different lines turn the body differently, and that difference is the moment.

Axes in this lesson: for the gear, x forward, y up and z to starboard. Check that it is right-handed: x̂ × ŷ = ẑ, and forward × up points to the right when you sit in the aircraft.

### 2.2 The moment of a force about a point

The **moment** of a force **F** acting at **r** about a point P is

**M_P = (r − P) × F**,   (1)

a vector in N·m. Its direction is the axis about which the force tries to turn the body (right-hand rule), and its magnitude is |r − P| |F| sin θ = |F| d, with d the perpendicular distance from P to the line of action (Baker & Haynes 4.1, 4.2). Three consequences follow from the cross product. **M_P is perpendicular to both the lever r − P and F.** **A force whose line passes through P has no moment about P** (sin θ = 0). And sliding the force along its line changes r − P by a multiple of F, which adds nothing to the cross product: the moment depends on the line, not on the point of application.

In a plane (z = 0) only the z component survives, and it is the familiar scalar

M_z = (x − P_x) F_y − (y − P_y) F_x,   (2)

positive anticlockwise seen from +z. The code keeps both forms: `moment_about_point` for 3-D and `moment_2d` for planar problems.

**Varignon's theorem** (4.3): the moment of a sum of forces acting at one point equals the sum of their moments, r × (F₁ + F₂) = r × F₁ + r × F₂. It is the distributive law of the cross product, and it is why you may resolve a force into components and take the moment of each component with its own simple lever arm.

A worked example in 3-D. The tyre contact point of the lesson's gear is r = (0, −1.2, 0.3) m from the trunnion (1.2 m below and 0.3 m outboard), and the ground force there is F = (−1000, 3000, 500) N. Then

M = r × F = (r_y F_z − r_z F_y, r_z F_x − r_x F_z, r_x F_y − r_y F_x) = (−1500, −300, −1200) N·m.

### 2.3 The moment about an axis

A pin lets a body turn about one line only, so what matters at a pin is the part of the moment along that line. The **moment about an axis** through P along the unit vector û is

**M_u = û · ((r − P) × F)**.   (3)

It does not depend on which point P of the axis you use: replace P by P + λû and the lever changes by −λû, whose cross product with F is perpendicular to û and drops out of the dot product. Two special cases are worth remembering: a force **parallel** to the axis, or a force whose line **cuts** the axis, has no moment about it. In the example above the moment about a fore-aft (x) axis through the trunnion is M_x = −1500 N·m, and the same axis through (5, 0, 0) gives −1500 N·m again.

### 2.4 Couples, force-couple systems and the central axis

Two equal and opposite forces F and −F on parallel lines form a **couple** (4.5). Their resultant is zero, but their moment is (r₁ − r₂) × F about **every** point: a couple is a free vector, a pure turning effect.

Any set of forces F_i at r_i and couples C_j is **statically equivalent** (4.6, 4.7) to a single force R at a chosen point A together with a couple M_A:

R = Σ F_i,   M_A = Σ (r_i − A) × F_i + Σ C_j.   (4)

Move the reference to B and the force is unchanged while the couple gains the moment of R about B:

**M_B = M_A + (A − B) × R.**   (5)

Equation (5) is the single most useful identity in statics: it says that once ΣF = 0, the moment sum is the same about every point, so you may take moments wherever the most unknowns drop out.

Is there a best point? Split M_A into a part parallel to R and a part perpendicular to it. Equation (5) can cancel the perpendicular part by choosing B suitably, but it can never change the parallel part, because (A − B) × R is always perpendicular to R. The set of points where only the parallel part remains is a line, the **central axis**; the point of it closest to the origin and its **pitch** are

p = R × M_O / |R|²,   pitch = R · M_O / |R|²,   (6)

and about any point of the axis the moment is pitch × R. The system is then a **wrench**: a force and a couple about the same line, like a screwdriver pushing and twisting. A **planar** system has R in the plane and M_O perpendicular to it, so R · M_O = 0 and the pitch is zero: every planar system with R ≠ 0 reduces to **one force on one line**. Example: the three forces (−1000, 3000, 500) N at (0, −1.2, 0.3), (0, 0, −400) N at the origin and (200, −100, 0) N at (0.2, 0.5, −0.1) m have R = (−800, 2900, 100) N, M_O = (−1510, −320, −1320) N·m and a pitch of 0.0163 m.

### 2.5 Distributed loads

Air pressure on a wing, the weight of a spar and fuel in a tank are spread over length or area. A line load w(x) in N/m has the **resultant** and line of action (Baker & Haynes 7.8, used again in `P05.09`)

W = ∫ w dx,   x̄ = ∫ x w dx / W.   (7)

Statically, the distributed load is equivalent to W acting at x̄, the centroid of the area under w (lesson `P05.05` generalises this). Closed forms: a uniform load w₀ on a length L gives W = w₀L at L/2; a triangular load rising from 0 to w₀ gives W = w₀L/2 at 2L/3 from the zero end. An **elliptic** spanwise lift on a semispan s, w(y) = w₀√(1 − (y/s)²), has W = πw₀s/4, so w₀ = 4W/(πs), and its centroid is

ȳ = 4s/(3π) = 0.4244 s.   (8)

(Derive it in Problem 5: with y = s sin u the integrals become trigonometric.) For the 1100 kg aircraft at 1 g with a 5 m semispan, each half-wing carries 5 393.66 N at 2.1221 m from the root.

Numerically, (7) is a pair of integrals of sampled data: the trapezoidal rule is exact for a linear w, second order for a smooth one. The elliptic load is not smooth at the tip, where w falls like √(s − y), and the error of x̄ shrinks more slowly than h²: 2.0831 m with 11 points, 2.1175 m with 41, 2.1215 m with 161 and 2.1220 m with 641, against the exact 2.1221 m.

### 2.6 The free-body diagram of the gear leg

Side view, x forward and y up, with the tyre contact point C at the origin:

```
                 E (0.60, 1.10)  airframe pin of the drag brace
                /
  A (0, 1.20)  /      A: trunnion pin, reacts Ax and Ay (no moment)
     o        /
     |       /  brace: two-force member, force T along B -> E
     |      /
     o B (0, 0.70)
     |
     |  leg
     |
     o C (0, 0)  tyre contact: ground force (-D, V)
  ---------------- runway
```

Cut the leg free from the runway, the trunnion and the brace. Each cut becomes a force:

- **Ground force at C**: (−D, V), V up and D aft (spin-up drag).
- **Trunnion pin at A**: a frictionless pin passes a force in any direction but no moment: two unknowns A_x, A_y.
- **Drag brace at B**: the brace is pinned at both ends and nothing loads it in between. Its own free-body diagram has only two forces, one at each pin, and for it to be in equilibrium those two forces must be equal, opposite and **along the line joining the pins** (two forces with different lines would form a couple). So the brace is a **two-force member** and its force on the leg acts along B→E with one unknown, T. Take T > 0 for tension, which pulls B toward E.

Three unknowns, and a planar rigid body has three equilibrium equations, ΣF_x = 0, ΣF_y = 0, ΣM = 0 (OpenStax 12.1). Before solving, note what to expect from the picture: the drag pushes the bottom of the leg aft, which would swing it about A; only the brace has a lever about A to resist that.

### 2.7 Solving it twice

**By one moment equation.** Take moments about A. Both trunnion reactions pass through A and drop out; V passes through A too, because C is directly below A. What is left is the drag and the brace:

Σ M_A = −D ℓ + T d = 0  ⇒  **T = D ℓ / d**,   (9)

with ℓ = 1.2 m the height of A above C and d the perpendicular distance from A to the brace line. The brace direction is u = (0.60, 0.40)/0.7211 = (0.8321, 0.5547), and with B − A = (0, −0.5) m, d = |(B − A) × u| = 0.5 × 0.8321 = **0.4160 m**. So T = 2.8844 D = **10 501 N, tension**.

**By force balance.** ΣF_x = 0: −D + A_x + T u_x = 0 gives A_x = D − T u_x = 3 641 − 8 738 = **−5 097 N**. ΣF_y = 0: V + A_y + T u_y = 0 gives A_y = −V − T u_y = −14 563 − 5 825 = **−20 388 N**. The trunnion pushes the leg forward and down. Its resultant, **21 015 N**, exceeds the vertical ground load, because the brace that holds the drag also pulls the leg up into the trunnion.

**As a linear system.** Each unknown force of known direction u at point r contributes a column [u_x, u_y, (r − O) × u] to a 3 × 3 matrix A; the known loads give b = (ΣF_x, ΣF_y, ΣM_O); equilibrium is A x = −b. This is what `solve_fbd_2d` does, and it is the seed of the general reaction solver of `P05.02`. The matrix is singular when the three unknown lines meet in one point or are all parallel, because then no combination of them can produce a moment about that point (or a force across those lines): the supports are improperly placed, and the solver says so.

**The load case.** V = n × (share of the weight on one main gear) × g. With 90 % of 1100 kg on the two mains, the static load per main gear is 4 854 N; an illustrative landing factor n = 3.0 gives V = 14 563 N, and a spin-up drag of 0.25 V gives D = 3 641 N. Real landing-gear load cases (spin-up, spring-back, braked roll, side load, one-wheel landing) come from the applicable airworthiness regulations and drop tests; they are not reproduced here.

### 2.8 Where the model stops being valid

- **Dynamic loads.** A touchdown lasts a fraction of a second. The spin-up drag peaks and then reverses as the strut springs back ("spring-back"), and the leg vibrates. A static FBD at the peak load is the first step, not the last.
- **Rigid, frictionless pins.** Real pins have friction (`P05.04`) and some fit-up play; a welded or bolted joint transmits moment and is not a pin.
- **The two-force brace.** If anything is attached to the brace between its pins (a retraction actuator, a down-lock link), it is a multi-force member and carries bending too.
- **Small deflections.** The geometry is taken as drawn; a soft gear that deflects changes the levers.

## 3. Build the computational model

Open `starter/main.py`. The helpers `unit`, `to3` (pads 2-D vectors with z = 0), `force_vector`, `ground_force`, the dataclass `GearGeometry` and `main` are given.

### Step 1: moments

```python
def moment_about_point(r, F, P=None):
    r3, F3 = to3(r), to3(F)
    P3 = np.zeros(3) if P is None else to3(P)
    return np.cross(r3 - P3, F3)
```

It works on one force or a stack of shape (n, 3) because `np.cross` broadcasts. `moment_2d` is (2) written out; `moment_about_axis` is (3) with û normalised inside the function, so the caller may pass any non-zero vector along the axis.

### Step 2: equivalent systems

`resultant(points, forces, couples, P)` returns (4); couples may be given as 3-vectors or, for a planar problem, as a single z value. `transport_moment` is (5); `central_axis` is (6) and raises `ValueError` for a pure couple (R = 0), which has no axis.

### Step 3: distributed loads

`distributed_load_resultant(x, w)` evaluates (7) with `np.trapezoid` and raises `ValueError` if W = 0 (a load whose resultant vanishes is a couple, not a force, and has no line of action). `elliptic_lift(y, s, L)` returns w₀√(1 − (y/s)²) with w₀ = 4L/(πs); clip the argument of the square root at zero so that rounding at the tip does not produce NaN.

### Step 4: the free-body solver and the gear

`equilibrium_residual_2d` sums forces and moments about a reference point. `solve_fbd_2d` builds the 3 × 3 matrix of §2.7, checks its condition number (above 10¹⁰ it raises `ValueError`) and solves. `gear_reactions` calls it with the ground force and the three unknowns (A_x, A_y at the trunnion, T along the brace). `brace_force_by_moments` is (9) written generally: T = −M_A(ground)/M_A(unit brace force). `brace_sweep` re-solves the gear for a list of brace attachment heights and returns NaN where the brace line passes through A.

Running `main()` prints the 3-D moment example, the gear solution by both routes with the residual about three points, the zero-drag case, the elliptic-load convergence, and the experiment tables, and writes `landing_gear_fbd.png` (the leg's FBD with force arrows, the brace force against attachment height, the elliptic load with its resultant line).

### Verify

Run `python -m pytest curriculum/phase-05/lesson-01-forces-moments-free-body-diagrams`. The tests check:

1. **Moments**: M ⊥ lever and force; |M| = |r||F| sin θ; zero for a force through P; the 2-D form equals the z component of the 3-D one, anticlockwise positive.
2. **Varignon**: the moment of a sum equals the sum of moments.
3. **Axis moment**: independent of the point on the axis and of the length of û; zero for a force parallel to the axis or cutting it.
4. **Transport and couples**: (5) against a direct recomputation; a couple has the same moment about every point; a free couple adds to M only.
5. **Central axis**: the moment about it is parallel to R; a planar system has pitch 0 and R placed on the axis reproduces M_O.
6. **Distributed loads**: uniform and triangular closed forms; the elliptic total and centroid converge monotonically (order above 1.2 despite the square-root tip).
7. **Solver**: zero residual about three different points; the answer does not depend on the reference point; concurrent or parallel unknowns are refused.
8. **Gear**: the 3 × 3 solve equals (9) for random loads; the FBD is in equilibrium about four points; A_y balances the vertical loads; the brace is in tension under aft drag.
9. **Superposition and limits**: the solution for (V, D) is the sum of those for (V, 0) and (0, D); doubling the loads doubles the answer; with D = 0 the brace carries nothing and A_y = −V.
10. **Brace sweep**: the force grows without bound as B approaches A, with T × d constant; a brace along the leg is refused.
11. `main` writes the figure.

### Compare

The independent check is the one-equation moment route (9) against the general linear solve: they agree to 0 N for the lesson case, because (9) is exact, and to 10⁻¹² relative for random loads in the tests. The residual of the full FBD is below 10⁻¹¹ N about the trunnion, the contact point and the point (3, −2) m. These checks catch the classic errors at once: a brace direction reversed (T comes out negative, compression, with the same magnitude), a moment arm measured to the wrong point (the two routes disagree), and a missing force (the residual is not zero).

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** The results below are computations with the model, not measurements of a real landing gear. Write each prediction down before running it.

1. **Brace attachment height.** Keep E fixed and move B up the leg from h = 0.05 m to 1.15 m. Predict, from (9), how T varies with the lever d(h), and where it blows up.
2. **Drag ratio.** For D/V from 0 to 0.8, predict the brace force and the trunnion resultant. Which one is linear in D/V, and why is the other one not?
3. **Geometric tolerance.** Move B and E each by up to ±5 mm in x and y (2 000 seeded random samples, `default_rng(501)`). Predict the spread of T from the derivative of (9) before you run it.

Raw data format: CSV tables `h_m, brace_N` (1); `drag_ratio, brace_N, trunnion_N` (2); `sample, Bx, By, Ex, Ey, brace_N` (3); plus `landing_gear_fbd.png`.

What you should find (model results): (1) T rises slowly at first and then steeply: 7.66 kN at h = 0.05 m, 9.12 kN at 0.55 m, 11.24 kN at 0.75 m, 18.01 kN at 0.95 m, 29.23 kN at 1.05 m and 87.68 kN at 1.15 m. The product T × d stays at D ℓ = 4 369 N·m; as B approaches the trunnion the brace line passes ever closer to A and d → 0. A brace attached low on the leg is lightly loaded but long and heavy; one attached high is short and heavily loaded. (2) T is exactly linear, 4.20 kN per 0.1 of D/V, and zero at D/V = 0 because V passes through A; the trunnion resultant is not linear (14.56 kN at 0, 22.40 kN at 0.3, 36.99 kN at 0.8) because it is the magnitude of a vector whose components grow at different rates. (3) T stays between 10 359 and 10 655 N, a spread of 2.82 % of the nominal 10 501 N: the brace load is insensitive to millimetre tolerances, as it should be in a sound design.

## 5. Validate and improve

| quantity | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| brace lever d about A (m) | 0.4160 | 0.4160 | 0 | — |
| brace force T (N) | 10 501 (from 9) | 10 501.4 | 0 | exact formula |
| trunnion A_x, A_y (N) | −5 097, −20 388 | −5 097.0, −20 388.0 | 0 | — |
| residual about A, C and (3, −2) (N, N·m) | 0 | below 10⁻¹¹ | round-off | — |
| brace force with D = 0 (N) | 0 | 0.0 | 0 | V passes through A |
| elliptic centroid, 641 points (m) | 2.12207 | 2.12200 | 3 × 10⁻⁵ relative | trapezoid at the √ tip |
| moment of the 3-D example (N·m) | (−1500, −300, −1200) | (−1500, −300, −1200) | 0 | — |

Every quantity with a closed form matches; the only numerical error is the trapezoidal rule at the tip of the elliptic load, which a change of variable removes (★ extension). The model's weaknesses are physical, not numerical: the landing is dynamic, the pins are not frictionless, and the load case is illustrative. The first improvement would be a load case from a drop-test or a dynamic model of the oleo strut (Phase 06), the second a 3-D gear with a side brace (`P05.02`).

## 6. Defend your work

Submit: the completed `main.py` with passing tests; your free-body diagram of the leg, drawn before the equations; the derivations of (5), (6) and (9); `landing_gear_fbd.png`; the three experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) answering the engineering question: the brace force and its sense, the trunnion load, and the sensitivity to the brace geometry.

### Adversarial self-checks

- **What if the load doubles?** Every force doubles (21 003 N in the brace): the problem is linear, and the test of superposition proves it. What does not double is the margin: a brace that was adequate may now buckle or yield, which `P05.06` and `P05.15` quantify.
- **What if the drag goes to zero?** The brace idles, and the trunnion carries exactly V. That is a property of this geometry (C directly below A), not of landing gears: a raked leg with the contact point ahead of the trunnion has a brace load even with no drag (Problem 4).
- **What if the drag reverses?** During spring-back the tyre force points forward and the brace goes into **compression**. A brace must be checked both ways; a member sized as a tie can buckle as a strut.
- **Which assumption fails first in the real world?** The quasi-static load case: dynamic overshoot and spring-back change both the magnitude and the sign of D.
- **What measurement could falsify the model?** A strain gauge on the brace during a drop test, converted to force (`P05.06`), compared with D ℓ/d using the measured drag; a large disagreement would point at pin friction, a moment-carrying joint or a misjudged drag.
- **Does the result survive a change of method?** Yes: one moment equation, the 3 × 3 solve and the force balance give the same three numbers, and the residual vanishes about any point.

## Check yourself

1. A force of 2 kN acts along a line 0.3 m from point P. What is the magnitude of its moment about P, and what happens to it if the force slides along its line?
2. Why is the moment of a force about an axis independent of the point you choose on the axis?
3. State the transport rule for a moment. What does it tell you about ΣM once ΣF = 0?
4. What is a couple, and why is its moment the same about every point?
5. Where does a uniform line load act? A triangular one? An elliptic one on a semispan s?
6. Why does a drag brace pinned at both ends carry its force along the line of its pins?
7. In the gear FBD, why does the trunnion reaction drop out of ΣM_A, and what is left?
8. When does a 3 × 3 free-body solve fail, and what does that failure mean physically?
