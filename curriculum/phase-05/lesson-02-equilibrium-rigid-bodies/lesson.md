# P05.02  Equilibrium of Rigid Bodies

> Six equations, six unknowns, and a rank that tells you whether the supports can do their job.

| | |
|---|---|
| **Learning outcomes** | 1. Model common supports (roller, pin, fixed support, two-force link, ball-and-socket, hinge) by the reaction components they can provide, and write the equilibrium of a rigid body as the linear system A x + b = 0. 2. Classify a support set by counting and by the rank of A as determinate, statically indeterminate (with its degree), partially constrained or improperly constrained, and explain each case physically. 3. Solve 2-D and 3-D reaction problems by hand with well-chosen moment equations (about a point and about an axis) and by a general solver, and show the two agree. 4. Verify a 3-D solution independently by virtual work, with link-length changes measured from the moved geometry. 5. Explain why aerospace mounts are made statically determinate, and quantify how a mount degrades as it approaches an improper constraint. |
| **Prerequisites** | `P05.01` (moments, free-body diagrams), `P01.13` (linear systems), `P01.16` (eigenvalues; the singular values used here generalise them), `P04.01` (Newton's laws). Skill check: you can say what the rank of a matrix is and compute it with `numpy.linalg.matrix_rank`. |
| **Engineering question** | A strut-braced high wing of a 1100 kg aircraft at its limit load factor of 3.8: what do the root pin and the strut carry? An aileron on two hinges and a pushrod: what does each hinge carry? An engine on six links: is the mount determinate, what does each link carry, and what happens if a link is added, removed, or turned until every link cuts the thrust axis? |
| **Assumptions and limits** | Rigid bodies and small displacements (the geometry is taken as drawn); frictionless joints; two-force links pinned at both ends and unloaded between them; static loads (inertia loads included as load factor × weight). Statically indeterminate supports cannot be solved here: they need deformation and compatibility (`P05.03`, `P05.07`, `P05.16`). |
| **Analytical work** | Support models; equilibrium in 2-D (3 equations) and 3-D (6 equations); the equilibrium matrix; counting and rank; the four constraint classes; the strut force from moments about the root; the pushrod force from the moment about the hinge line; the principle of virtual work for a rigid body; why a near-singular mount carries large forces. |
| **Computational deliverable** | `starter/main.py`: `equilibrium_matrix`, `load_vector`, `classify_supports`, `solve_reactions`, `residual`, `strut_wing_reactions`, `strut_force_closed_form`, `aileron_reactions`, `hinge_moment`, `pushrod_force_by_axis_moment`, `link_lengths`, `virtual_work_residual`, `spread_sweep` (the `Unknown` dataclass, the support constructors `roller`, `pin`, `fixed_2d`, `link`, `ball_socket`, `hinge`, `fixed_3d`, the `StrutWing` and `Aileron` dataclasses, `hexapod_mount`, `engine_loads` and `main` are given). Tests prove: matrix columns and reference-point changes; the classification of beams, 3-D supports and mount variants; refusal of non-determinate sets; lever-rule beams; the strut wing against its closed form; the aileron against the hinge-line moment; the mount against virtual work; reversal, superposition and the near-improper blow-up. |
| **Practical activity** | **computational experiment**: a sweep of the strut attachment station of the braced wing, a classification table of support sets, and a sweep of the engine-mount link spread toward an improper constraint with the condition number and the largest link force. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `equilibrium.png`, hand solutions of the strut force and the pushrod force, the classification table, the spread sweep, and a half-page note answering the engineering question. |
| **Assessment** | `assessments/phase-05/grading-rubric.md`; competencies 1–4. Pass: every support is drawn with its reaction components before the equations; the strut and pushrod forces are found by hand with one equation each and match the code; all public tests pass; the classification of every support set is justified physically, not only by the rank. |
| **References** | `@libretexts-baker-haynes-engineering-statics` Ch. 5 "Rigid Body Equilibrium": 5.1 Degree of Freedom, 5.2 Free Body Diagrams, 5.3 Equations of Equilibrium, 5.4 2D Rigid Body Equilibrium, 5.5 3D Rigid Body Equilibrium, 5.6 Stability and Determinacy, 5.7 Equilibrium Examples (required-free: support types, the equations and the determinacy discussion); `@openstax-university-physics-vol1` Ch. 12: 12.1 Conditions for Static Equilibrium, 12.2 Examples of Static Equilibrium (required-free: a gentle review); `@mit-ocw-16-001-unified-materials-structures` lecture notes "Structural Supports and Reactions; Requirements for Structural Equilibrium" (required-free: the same ideas from an aerospace department). The aerospace bodies, the rank classification code and the virtual-work check are original to this lesson. |
| **Extensions** | ★ Add a side brace to the landing gear of `P05.01` and solve it in 3-D with a side load. ★ Write the virtual-work check as the transpose relation: link elongation rates = Aᵀ (rigid-body velocity), and show it numerically. ★ Solve the seven-link mount as a statically indeterminate problem by giving each link a stiffness EA/L (the method of `P05.03`) and show how the redundant link's force depends on the stiffnesses. |
| **Time** | unmeasured |

## 1. The engineering problem

A structure is held by its supports, and the supports are where loads leave it. Before anything can be sized, every reaction must be known: the force in a wing strut, the load on each hinge of a control surface, the force in each link of an engine mount. Equilibrium alone gives them when the supports are **statically determinate**: exactly as many independent reaction components as there are equations of equilibrium. With fewer, the body moves; with more, equilibrium alone cannot share the load; and with the right number badly arranged, the body moves anyway.

This lesson writes equilibrium as a linear system A x + b = 0, with one column of A per unknown reaction, and lets the rank of A classify the supports. It applies the solver to three aerospace bodies. A **strut-braced wing** (the high-wing light aircraft layout) at a load factor of 3.8: the strut carries **40 611 N** in tension and the root pin **38 082 N**, both far more than the 20 496 N of lift on the half-wing, because the strut is inclined and the pin must balance its horizontal pull. An **aileron** on two hinges and a pushrod: the pushrod carries **1 200 N**, which one moment about the hinge line gives directly. And an **engine on six links**: determinate, every link loaded between 0.39 and 4.31 kN, verified by virtual work; add a seventh link and it becomes indeterminate, remove one and it becomes a mechanism, and turn the links until each one cuts the thrust axis and it cannot hold the propeller torque at all.

Success criteria, written first: zero force and moment residual about any point; the strut and pushrod forces equal to their one-equation hand solutions; the mount passing a virtual-work check whose ingredients do not use the equilibrium matrix.

## 2. Understand and derive

### 2.1 Supports and the reactions they provide

A support stops some motions of the body; for each motion it stops it provides one reaction component (Baker & Haynes 5.2). The rule: **a support resists a translation with a force and a rotation with a couple**, in exactly the directions it prevents.

| support | prevents | unknowns |
|---|---|---|
| roller / rocker (2-D) | motion along the surface normal | 1 force along the normal |
| two-force link (2-D or 3-D) | stretching along its line | 1 force along the link |
| pin (2-D) | both translations | 2 forces |
| fixed support (2-D) | both translations and rotation | 2 forces + 1 couple |
| ball-and-socket (3-D) | all translations | 3 forces |
| hinge / journal bearing (3-D) | translations normal to the pin (+ axial for a thrust bearing) | 2 (or 3) forces; 2 couples only if the bearing resists tilting |
| fixed support (3-D) | everything | 3 forces + 3 couples |

A body supported on two well-aligned hinges is normally modelled with **no** hinge couples: each hinge alone could resist tilting only through its short bearing length, and including those couples would make the body indeterminate for no physical gain.

### 2.2 Equilibrium as a linear system

A rigid body is in equilibrium when ΣF = 0 and ΣM_O = 0 about any point O: six scalar equations in 3-D, three in a plane (Fx, Fy, Mz) (5.3). By the transport rule of `P05.01`, once ΣF = 0, ΣM is the same about every point, so O is free.

Each unknown reaction x_k is a force along a unit vector u_k through a point r_k, or a couple about u_k. Its contribution to the six sums is x_k times a fixed column:

force: **[u_k ; (r_k − O) × u_k]**,  couple: **[0 ; u_k]**.   (1)

Collect the columns into A (6 × n) and the applied loads into b = [ΣF ; ΣM_O]. Equilibrium is

**A x + b = 0.**   (2)

In 2-D keep the rows F_x, F_y, M_z. Changing O changes the moment rows of A and b by (−O) × (force rows), the same row operation on both sides, so the solution does not depend on O (a test checks this).

### 2.3 Counting, rank and the four classes

Let n_eq be the number of equations (3 or 6), n the number of unknowns and ρ the rank of A, the number of independent columns (P01.13). Two numbers say everything:

- **n − ρ** = the number of independent **self-stress** states: combinations of reactions that balance each other with no load. If n − ρ > 0, equilibrium cannot decide how much of each the body carries: the body is **statically indeterminate** of that degree.
- **n_eq − ρ** = the number of independent **mechanisms**: rigid motions the supports do not stop. If n_eq − ρ > 0, a general load cannot be balanced.

So (Baker & Haynes 5.6):

| class | condition | example |
|---|---|---|
| determinate | ρ = n_eq = n | beam on a pin and a roller |
| statically indeterminate | ρ = n_eq < n | beam on two pins (degree 1) |
| partially constrained | ρ < n_eq, n < n_eq | beam on two vertical rollers (slides sideways) |
| improperly constrained | ρ < n_eq, n ≥ n_eq | beam on three vertical rollers; three concurrent links |

The last class is the dangerous one, because counting says "enough". Three parallel rollers cannot stop a horizontal slide; three links whose lines meet at one point cannot stop a rotation about that point. In both, the count is three and the rank is two. And since n − ρ = 1 too, an improperly constrained body is also indeterminate: the rank deficit that lets it move is the same deficit that lets two of its reactions fight each other.

`classify_supports` computes ρ from the singular values of A (those above 10⁻⁹ times the largest), and also returns the **condition number**, the ratio of the largest to the smallest singular value. A large condition number means a nearly improper support: the body can be held, but only by large reactions (§2.7).

### 2.4 Example 1: a strut-braced wing (2-D)

Front view of the left half-wing, x outboard along the span and y up, root pinned to the fuselage at the origin. A strut runs from the fuselage at (0, −1.0) m to the wing at y_s = 2.6 m; the semispan is s = 5 m. At the limit load factor n = 3.8 (typical of a light aircraft in the normal category) the lift on one half-wing is L = ½ × 3.8 × 1100 × 9.80665 = **20 496 N**, elliptic, acting at 4s/(3π) = 2.1221 m (`P05.01`). The half-wing's own weight, 60 kg, is multiplied by the same load factor (it is accelerated with the aircraft) and acts down at s/2: n W = **2 236 N**.

Unknowns: the root pin's two components and the strut force T (a two-force member, tension positive). Moments about the root pin remove both pin components:

T × (lever of the strut about the root) = L ȳ − n W s/2.

The strut has length ℓ_s = √(2.6² + 1.0²) = 2.7857 m; its lever about the root is the perpendicular distance from the root to its line, y_s × h/ℓ_s = 2.6 × 1.0/2.7857 = **0.9333 m**. So

**T = (L ȳ − n W s/2) ℓ_s/(y_s h)** = (43 494 − 5 590)/0.9333 = **40 611 N** (tension).   (3)

The force balance gives the pin: the strut pulls the wing toward the fuselage point, with components T(−2.6, −1.0)/2.7857 = (−37 904, −14 578) N, so the pin supplies **(37 904, −3 682) N**: a large outboard push and a small downward pull. A strut is a very efficient way to carry wing bending, but only because the fuselage pin takes a horizontal force almost as large as the strut force.

At the negative limit load factor (−1.52 for this category, 0.4 times the positive one) the lift reverses and the strut is pushed: **T = −16 244 N, compression**. A strut must be designed as a column too (`P05.15`); light-aircraft struts often have a jury strut part-way along to shorten the buckling length.

### 2.5 Example 2: an aileron on two hinges (3-D)

Axes x forward, y spanwise, z up; the hinge line is the y axis. The inboard hinge at y = 0.10 m is a ball-and-socket (3 forces); the outboard hinge at y = 1.40 m lets the aileron slide along the hinge line (2 forces, x and z), so that thermal expansion and wing bending do not load the aileron spanwise. A pushrod is pinned to a horn 0.06 m below the hinge line at y = 0.75 m and runs forward (a link along x). Six unknowns, six equations: determinate.

The aerodynamic load is a normal force N = 900 N (up) acting 0.08 m aft of the hinge line at y = 0.70 m, plus a chordwise force of 40 N aft. Its moment about the hinge line, the **hinge moment**, is M_y = z F_x − x F_z = −(−0.08)(900) = **72.00 N·m**.

**One equation about an axis.** Both hinge reactions pass through the hinge line, so they have no moment about it (`P05.01`, §2.3). Only the pushrod does: a tension P along +x at z = −0.06 m has moment (−0.06)P about y. So

−0.06 P + 72.00 = 0  ⇒  **P = 1 200 N** (tension).   (4)

The general solver gives the same 1 200 N and the remaining five reactions: ball joint (−578.5, 0, −484.6) N, outboard hinge (−581.5, −415.4) N in x and z. The spanwise reaction is zero because nothing loads the aileron spanwise. The rule generalises: **to find one unknown, take moments about an axis that every other unknown cuts or is parallel to**.

### 2.6 Example 3: an engine on six links (3-D)

An engine mount that is statically determinate is called **isostatic**: six links, each pinned at both ends, arranged so that A is square and well conditioned. The mount of `hexapod_mount` has three engine points on a 0.25 m radius around the thrust axis x, and each is joined to two firewall points 0.40 m behind: one at 0.34 m radius turned +40° around the axis, one at 0.16 m radius turned −40°. Loads: a 120 kg engine at a load factor of 4 (4 707 N down) at its centre of gravity 0.25 m ahead of the attachment plane, 3 000 N of thrust and a propeller-torque reaction of −450 N·m about x.

The solver finds the mount determinate (rank 6, condition number 9.98) and the link forces **4 310, 386, −4 265, 3 428, 3 603 and −4 031 N** (tension positive). Why the mount is built this way, and not with more links:

- **Determinate mounts are stress-free under deformation.** If the firewall bends, or the engine grows 1 mm hotter than the airframe, a determinate mount simply moves: there is no self-stress state to excite. An indeterminate one fights the deformation and loads its links even with no external load.
- **Load paths are known.** Each link's force follows from equilibrium alone, so it does not depend on manufacturing tolerances or joint stiffnesses.

The variants show the four classes: with a **seventh link** the rank is still 6, n − ρ = 1: indeterminate of degree 1. With **any one link removed**, ρ = 5 < 6: partially constrained, one free motion. With the link spread set to **0°**, every link lies in a plane containing the thrust axis, so each one cuts or is parallel to that axis and has no moment about it: n = 6 but ρ = 5, improperly constrained. The mount then holds any load except a torque about x, which no combination of links can provide.

### 2.7 Near an improper constraint

What happens between 40° and 0°? Torque about x can only come from the links' tangential components, which scale as sin(spread). To make 450 N·m the links must pull harder as the spread closes: the largest link force under torque alone is 851 N at 40°, 1 482 N at 20°, 2 855 N at 10°, 5 656 N at 5° and 14 102 N at 2°, while its product with sin(spread) settles at about 492 N. The condition number grows the same way, from 9.98 at 40° to 183 at 2°. **A support set can be determinate and still useless**: the condition number, not the count, says whether it is a good design.

### 2.8 An independent check: virtual work

The principle of virtual work states that a rigid body is in equilibrium if and only if the total work of all forces on it vanishes for every small rigid displacement. Move the engine by a small translation εt and rotation εθ. A point moves by ε(t + θ × r), a load F at r does work ε F · (t + θ × r), a couple C does ε C · θ, and a link of tension T_i that lengthens by δℓ_i does work −T_i δℓ_i (the link pulls the engine back toward the firewall). The check

W(t, θ) = −Σ T_i δℓ_i/ε + Σ F · (t + θ × r) + Σ C · θ = 0   (5)

uses **no equilibrium matrix**: δℓ_i is measured by moving the engine points with an exact rotation (Rodrigues' formula), recomputing each link length, and differencing (`link_lengths`). For the lesson's mount, (5) is below 2 × 10⁻⁷ N·m for five random rigid motions, against forces of thousands of newtons on levers of tenths of a metre; a wrong link force (one sign flipped) fails it by orders of magnitude. The duality behind it: the matrix that turns rigid-body velocities into link elongation rates is Aᵀ, the transpose of the equilibrium matrix, which is why a mechanism (a motion with no elongation) and a self-stress (forces with no resultant) are counted by the same rank.

### 2.9 Where the model stops being valid

- **Rigid bodies.** Wings, firewalls and engines deform. For a determinate support set this changes nothing to first order; for an indeterminate one it decides the answer.
- **Pinned links.** Real rod ends have friction and some bending stiffness; elastomeric engine-mount isolators are springs, not pins, and make a mount indeterminate (with intended flexibility).
- **Small displacements.** The geometry is taken undeformed. Near an improper constraint, small deflections change the levers a lot, and the linear answer can be badly wrong.
- **Static loads.** Inertia loads are included as load factor × weight; vibration and gyroscopic moments of a spinning propeller are not.

## 3. Build the computational model

Open `starter/main.py`. The `Unknown` dataclass (a name, a point, a unit direction and a kind, `"force"` or `"couple"`; 2-D points and directions are padded to 3-D), the support constructors, `StrutWing`, `Aileron`, `hexapod_mount`, `engine_loads` and `main` are given.

### Step 1: the equilibrium matrix and the load vector

```python
def equilibrium_matrix(unknowns, dim=3, ref=None):
    O = np.zeros(3) if ref is None else np.append(np.asarray(ref, float), [0.0] * (3 - len(ref)))
    cols = []
    for u in unknowns:
        if u.kind == "force":
            cols.append(np.concatenate([u.direction, np.cross(u.point - O, u.direction)]))
        else:
            cols.append(np.concatenate([np.zeros(3), u.direction]))
    A = np.array(cols).T
    return A[[0, 1, 5]] if dim == 2 else A
```

`load_vector` sums the applied forces and their moments about the same reference, adding any couples (3-vectors, or scalars M_z in 2-D).

### Step 2: classify and solve

`classify_supports(A)` takes the singular values with `np.linalg.svd(A, compute_uv=False)`, counts those above 10⁻⁹ times the largest, and returns the counts, the two degrees, the status and the condition number. `solve_reactions` raises `SupportError` (a `ValueError`) unless the status is determinate, and otherwise solves (2) and returns a dictionary keyed by the unknowns' names. `residual` returns A x + b about any reference point.

### Step 3: the three bodies

`strut_wing_reactions` and `aileron_reactions` call the solver; `strut_force_closed_form` is (3) and `pushrod_force_by_axis_moment` is (4) written generally. `link_lengths(eng, frm, t, θ)` rotates the engine points by Rodrigues' formula and translates them; `virtual_work_residual` is (5) with central differences of step 10⁻⁶. `spread_sweep` returns the largest link force and the condition number over a list of spreads.

Running `main()` prints the three examples, the variants and the sweep, and writes `equilibrium.png` (strut and root-pin forces against the strut station; the aileron reactions; the largest link force and the condition number against the spread).

### Verify

Run `python -m pytest curriculum/phase-05/lesson-02-equilibrium-rigid-bodies`. The tests check:

1. **Matrix columns**: [u; r × u] for a force with u normalised, [0; u] for a couple, rows (F_x, F_y, M_z) in 2-D.
2. **Reference point**: the moment rows change by −O × R and the solution does not change.
3. **Planar classification**: pin + roller determinate; two pins indeterminate; two vertical rollers partially constrained; three vertical rollers and three concurrent links improperly constrained; a fixed support determinate.
4. **3-D classification**: a fixed support, the aileron, two ball joints (free to spin about their line and indeterminate along it), and the four mount variants with their degrees.
5. **Hinges**: forces normal to the axis, one axial component only when asked.
6. **Refusal**: the solver raises `SupportError` for indeterminate and partial sets.
7. **Lever rule**: a simply supported beam under a point load and under a couple.
8. **Strut wing**: the solver equals (3); zero residual about three points; tension at positive and compression at negative load factor; superposition; an idle strut when lift and weight moments balance.
9. **Aileron**: the pushrod equals (4) for random loads; zero residual; no spanwise reaction.
10. **Mount**: the virtual-work residual vanishes for six random rigid motions and is large for a wrong answer; link lengths are unchanged by zero motion and a full turn.
11. **Reversal and limits**: reversed loads reverse every link force; zero load gives zero forces.
12. **Near-improper**: the largest force and the condition number grow as the spread closes; at 0° every load except the torque can be balanced, and the unbalanced residual is exactly the 450 N·m torque.
13. `main` writes the figure.

### Compare

The independent checks are the hand equations (3) and (4), the virtual-work residual (5) and the classification of textbook cases whose answers are known. The strut force agrees with (3) to 10⁻¹² relative, the pushrod with (4) exactly, and the virtual-work residual of the mount is 10⁻¹⁰ of the work scale. The commonest error, a sign slip in a moment column, shows up at once in all three: the residual about a second point is no longer zero.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** All results are computations with the rigid-body model, not measurements. Write each prediction down before running it.

1. **Strut station.** Move the strut's wing attachment from y_s = 1.0 m to 4.6 m. Predict, from (3), how T varies (its lever grows with y_s but saturates as the strut becomes nearly horizontal) and whether the root pin load has a minimum.
2. **Classification table.** For six support sets of a 4 m beam (pin + roller; fixed; two pins; three vertical rollers; two vertical rollers and one horizontal; a pin and a link whose line passes through the pin), predict the class by counting and by looking for concurrent or parallel lines, then run `classify_supports`.
3. **Link spread.** Sweep the mount's link spread from 40° down to 1°. Predict how the condition number and the largest link force scale.

Raw data format: CSV tables `y_s_m, strut_N, root_N` (1); `case, n, rank, status, indeterminate_degree, freedom_degree` (2); `spread_deg, max_link_N, cond` (3); plus `equilibrium.png`.

What you should find (model results): (1) T falls from 53.60 kN at y_s = 1.0 m to 40.61 kN at 2.6 m and 38.79 kN at 4.6 m; the root resultant has a shallow minimum of 37.90 kN near y_s = 2.08 m, where the strut's vertical component exactly balances the net upward load and the pin carries a purely horizontal force. The horizontal pin force itself, 37 904 N, does not depend on y_s at all: it is the root bending moment divided by the strut's depth h, (L ȳ − n W s/2)/h. Moving the strut outboard reduces its force but lengthens it (heavier, and more prone to buckle in compression). (2) Determinate: pin + roller, fixed, two vertical rollers + one horizontal; indeterminate (degree 1): two pins; improperly constrained: three vertical rollers, and a pin with a link through it (each with one free motion and one self-stress). (3) Under the full load the largest link force rises from 4.31 kN at 40° to 7.20, 10.65, 14.37, 23.12 and 37.26 kN at 20°, 10°, 5°, 2° and 1°, and the condition number from 10.0 to 18.7, 36.9, 73.4, 183 and 367. The condition number doubles each time the spread halves, as the 1/sin(spread) argument of §2.7 predicts; the largest link force grows more slowly at first, because most of the load is weight and thrust, which the links carry well at any spread, and approaches the same law as the torque term takes over.

## 5. Validate and improve

| quantity | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| strut force, n = 3.8 (N) | 40 611 (3) | 40 610.7 | 10⁻¹² relative | — |
| root pin (N) | (37 904, −3 682) | (37 903.9, −3 681.6) | 0 | — |
| strut force, n = −1.52 (N) | compression | −16 244 | — | sign check |
| hinge moment (N·m) | 72.00 | 72.00 | 0 | — |
| pushrod force (N) | 1 200 (4) | 1 200.00 | 0 | — |
| mount virtual-work residual (N·m) | 0 | below 2 × 10⁻⁷ | round-off of the difference | — |
| mount rank, 6 / 7 / 5 links / spread 0 | 6 / 6 / 5 / 5 | same | 0 | — |

Every hand result is reproduced, and the mount passes a check that does not share its equations. The model's weaknesses are the rigid-body and pinned-joint idealisations: the strut wing's root is in reality a bolted fitting with some moment capacity (which makes it indeterminate), and the engine mount's links are often attached through elastomeric isolators. The improvement for both is to add stiffness, which is what `P05.03` does for trusses.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; free-body diagrams of the three bodies with every reaction drawn; the derivations of (3) and (4); the classification table with a physical reason for every entry; `equilibrium.png`; the two sweep tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) answering the engineering question.

### Adversarial self-checks

- **What if the load factor doubles?** Every reaction doubles; the problem is linear. The classification does not change, because it depends only on the geometry of the supports.
- **What if the load goes to zero?** Every reaction of a determinate body goes to zero. An indeterminate body may still carry self-stress (from fit-up or temperature), which is exactly why aerospace mounts avoid indeterminacy.
- **Which assumption fails first in the real world?** For the strut wing, the moment-free root pin; for the mount, the rigid firewall and the pinned links. Both make the real structure stiffer and somewhat indeterminate.
- **What measurement could falsify the model?** Strain gauges on the strut in a ground load test at a known applied wing load: the measured force divided by the load should match (3); a lower value would mean the root fitting is carrying moment.
- **Does the result survive a change of reference point or method?** Yes: the solution is identical about any reference, and virtual work agrees without using A.
- **Which decision changes if the uncertainty doubles?** Near an improper constraint, a small geometric error changes the link forces a lot (they scale with the condition number); a mount with a condition number of 10 is robust, one of 300 is not.

## Check yourself

1. How many reaction components does a 2-D pin provide? A 3-D ball-and-socket? A hinge that lets the body slide along its pin?
2. Write the column of the equilibrium matrix for a unit force along u through r, and for a unit couple about u.
3. A body has n = n_eq unknowns but its matrix has rank n_eq − 1. What is it called, and what can it do?
4. Why does adding a seventh link to a determinate engine mount make it statically indeterminate rather than "safer"?
5. Which single equation gives the pushrod force of an aileron on two hinges, and why do the hinge reactions drop out of it?
6. Why is a wing strut in tension in positive-g flight and in compression in negative-g flight?
7. What does a large condition number of the equilibrium matrix tell a designer?
8. State the principle of virtual work for a rigid body and say why it is an independent check on a reaction solution.
