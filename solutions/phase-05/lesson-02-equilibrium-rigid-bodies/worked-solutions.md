# P05.02  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` (and the probe script that calls it) and given to the figures shown.

## 1. Classify, then solve (s = 15: P = 2500 N at x = 1.5 m)

**Assumptions.** Rigid beam, frictionless supports, weight of the beam neglected.

| case | n | geometry | class | degrees (indeterminate, free) |
|---|---|---|---|---|
| (a) pin + vertical roller | 3 | not concurrent, not parallel | determinate | 0, 0 |
| (b) fixed support | 3 | — | determinate | 0, 0 |
| (c) two pins | 4 | — | indeterminate | 1, 0 |
| (d) three vertical rollers | 3 | all parallel | improperly constrained | 1, 1 |
| (e) two vertical rollers + horizontal roller | 3 | not concurrent, not parallel | determinate | 0, 0 |
| (f) pin + link along the x axis | 3 | link line passes through the pin | improperly constrained | 1, 1 |

`classify_supports` agrees with every row.

(a) Moments about the pin: 4 B = 1.5 P, **B = 937.5 N**; vertical sum: **A_y = 1562.5 N**; horizontal: **A_x = 0**.
(b) **A_y = 2500 N, A_x = 0**, and the wall couple balances the load's moment −1.5 P: **A_M = +3750 N·m** (anticlockwise).
(e) The horizontal roller carries nothing (no horizontal load): **C = 0**; the vertical rollers carry **1562.5 N and 937.5 N** as in (a).

(d) can slide sideways and its three vertical reactions share any vertical load in a way equilibrium cannot decide; (f) can turn about the pin, because the link's line passes through it and gives no moment about it.

**Check.** In (a) and (e), A_y + B = P and the lever rule gives B/A_y = 1.5/2.5.

**Common mistakes.** Calling (d) or (f) determinate because the count is three; forgetting the wall couple in (b).

## 2. What the rank counts

(a) Let the three link lines meet at O. Each column of A about O is [u_k; 0], so the moment row is zero and the rank is at most 2 (it is 2 if two of the u_k are not parallel). **Free motion**: a small rotation about O, which changes no link length to first order. **Self-stress**: link forces T_k with Σ T_k u_k = 0, a closed force triangle, which is in equilibrium with no load.

(b) A maps reaction vectors (dimension n) to resultants (dimension n_eq). By rank–nullity, dim null(A) = n − ρ: these are reaction sets with zero resultant, i.e. **self-stresses**; any solution of A x = −b can have any of them added, so equilibrium does not fix x when n − ρ > 0. The range of A has dimension ρ; loads b with a component outside it cannot be balanced, and those components form a space of dimension n_eq − ρ (the null space of Aᵀ): these are the **mechanisms** (the duality of §2.8 in the lesson). A unique solution for every load needs n − ρ = 0 and n_eq − ρ = 0: **ρ = n = n_eq**.

(c) Moments about the root pin: lift L at ȳ (+), weight nW at s/2 (−), and the strut force T at (y_s, 0) along (−y_s, −h)/ℓ_s, whose moment is y_s × (−h/ℓ_s) T = −(y_s h/ℓ_s) T. Setting the sum to zero: **T = (L ȳ − nW s/2) ℓ_s/(y_s h)**. The horizontal component of the pin force balances the strut's: H = T y_s/ℓ_s = **(L ȳ − nW s/2)/h**, independent of y_s: the pin and the strut form a couple of arm h that resists the root bending moment.

**Common mistakes.** In (b), counting mechanisms by n_eq − n instead of n_eq − ρ (which misses improper constraints); in (c), using the strut length ℓ_s instead of the perpendicular lever y_s h/ℓ_s.

## 3. A cargo door held open by a gas strut (s = 15: m = 40 kg)

**Assumptions.** Rigid door; frictionless hinges; the gas strut is a two-force member (pinned at both ends, no load between); static gust load.

(a) Ball joint (3) + sliding hinge (2) + strut (1) = 6 unknowns; rank 6, condition number 6.65: **determinate**.

(b) Moments about the hinge line (x axis), M_x = y F_z − z F_y. Weight W = 392.27 N at y = 0.45: −176.52 N·m; gust 150 N at y = 0.6: −90.00 N·m; total **−266.52 N·m**. The strut runs from (0.6, 0.5, 0) to (0.6, 0, −0.6), length 0.7810 m, u = (0, −0.6402, −0.7682); a unit tension at the door end has M_x = 0.5 × (−0.7682) − 0 = −0.3841 m. So T = −(−266.52)/(−0.3841) = **−693.9 N: compression**. A gas strut pushes, as it should.

(c) `solve_reactions`: ball joint **(0, −222.10, 4.61) N**, sliding hinge **(y −222.10, z 4.61) N**, strut −693.86 N; residual below 10⁻¹⁴.

The strut pushes the door outboard with 444.2 N and up with 533.0 N; the hinges pull it back inboard (2 × 222.1 N) and supply the small remaining 9.2 N of vertical force between them.

(d) No load has a component along x, so the ball joint's x reaction is zero, and the sliding hinge cannot take one anyway. A gust with a fore-aft component, or a strut not perpendicular to x, would load the ball joint along x.

**Check.** Vertical: 533.0 + 2 × 4.61 = 542.3 N = 392.3 + 150. The two hinges share equally because all loads act at x = 0.6 m, midway between them.

**Common mistakes.** Taking the strut's moment arm as its length; forgetting that the sliding hinge has no x component (which would make the door indeterminate).

## 4. Size the strut station (s = 15: m = 1200 kg)

L = ½ × 3.8 × 1200 × 9.80665 = 22 359.2 N on the half-wing.

(a) y_s = 2.6 m: **T = 44 847.1 N** (tension), equal to `strut_force_closed_form`; root pin **(41 857.8, −4 024.1) N**.

(b) Over 1.0–4.6 m, the root resultant is smallest, **41 857.8 N, at y_s = 2.08 m**. The pin's horizontal force is fixed at (L ȳ − nW s/2)/h = 41 857.8 N whatever y_s is (Problem 2c), so the resultant is smallest where the vertical component vanishes: where the strut's vertical component exactly equals the net upward load L − nW.

(c) At n = −1.52 the strut carries **−17 938.8 N: compression**.

(d) A strut attached near the tip is long (heavy, and its buckling load in negative g falls as 1/length², `P05.15`), its force barely drops further (44.85 kN at 2.6 m, 42.84 kN at 4.6 m: the lever saturates as the strut flattens), and it adds drag and interferes with ailerons and flaps.

**Check.** At y_s = 1.0 m, T = 59 195.9 N; at 4.6 m, 42 835.5 N: T × (y_s/ℓ_s) is 41 857.8 N in both, the fixed horizontal force.

**Common mistakes.** Forgetting to reverse the wing's inertia load at negative g; applying the load factor to the lift but not to the wing weight.

## 5. The engine mount near an improper constraint

(a) Torque only, 450 N·m:

| spread | largest link force (N) | × sin(spread) (N) | cond(A) |
|---|---|---|---|
| 40° | 851.2 | 547.1 | 9.98 |
| 20° | 1 481.7 | 506.8 | 18.72 |
| 10° | 2 855.0 | 495.8 | 36.86 |
| 5° | 5 655.9 | 492.9 | 73.43 |
| 2° | 14 102.0 | 492.2 | 183.38 |

The torque is carried by the links' tangential components, which are proportional to sin(spread), so the link forces grow as 1/sin(spread); the product settles at about 492 N once the spread is small and the geometry otherwise stops changing. The condition number doubles each time the spread halves.

(b) The minimum-norm solution of the seven-link mount is (4 574.3, 1 424.1, −3 836.0, 3 021.6, 3 344.6, −3 785.5, −1 406.8) N; it balances the load to 4 × 10⁻¹² N. The null vector of A, (−0.1395, −0.5469, −0.2260, 0.2143, 0.1363, −0.1292, 0.7411), is a **self-stress**: link forces in these proportions balance each other. Any multiple can be added, so equilibrium allows infinitely many answers. The real mount picks one by **compatibility**: the links stretch in proportion to force × length/(EA), and those stretches must fit a single rigid motion of the engine. The minimum-norm answer corresponds to no physical stiffness distribution in particular; the six-link forces (4 309.5, 386.1, −4 265.0, 3 428.4, 3 603.3, −4 030.8) N with the seventh at zero are another valid equilibrium, and the truth lies wherever the stiffnesses put it.

(c) Dropping any one of the six links (k = 0…5) gives rank 5: **partially constrained** in every case. Each link of an isostatic mount is essential.

**Common mistakes.** Treating the `lstsq` answer as "the" solution; forgetting that a near-singular mount fails by excessive forces long before it becomes a mechanism.

## 6. Determinacy in aerospace practice

(a) An instrument on six links (or three bipods, each giving two constraints) is determinate: when the spacecraft structure distorts in launch or heats unevenly in orbit, the instrument moves as a rigid body and no stress reaches its optics. An indeterminate mount would bend the optical bench and spoil the alignment.

(b) A redundant link gives a second load path, so the engine stays attached if one link cracks (fail-safe design). It costs a statically indeterminate analysis: the load sharing depends on the link stiffnesses and on fit-up, thermal mismatch stresses appear, and the structure must also be checked with each link removed in turn.

(c) The beam cannot resist any horizontal load, however small: a slight slope, a vibration or a thermal expansion that pushes it sideways moves it. And for vertical loads it is indeterminate, so the share of each roller depends on how level the rollers are. The classification flags both weaknesses before a test does.

(d) Link forces of order (condition number) × (load / typical lever): a condition number of 300 means forces hundreds of times larger than the loads would suggest, so heavy links. And the forces are as sensitive to geometric errors as the matrix is ill-conditioned: a 1 % change in an attachment position can change some link forces by a large fraction.

## Check-yourself answers

1. Two; three; two (only the forces normal to the pin).
2. Force: [u; (r − O) × u]. Couple: [0; u].
3. Improperly constrained: it has a free rigid motion (a mechanism) under some loads, and at the same time a self-stress state.
4. The rank stays at six, so the seventh unknown cannot be found from equilibrium: the load share depends on stiffness and fit-up, and the mount picks up stress from thermal and structural deformation.
5. Moments about the hinge line: both hinge reactions cut the line and have no moment about it, so the pushrod force times its lever equals the hinge moment.
6. The strut joins the wing to the fuselage below it; positive lift pulls the wing up, stretching the strut, negative lift pushes the wing down, compressing it.
7. That the supports are nearly improper: some loads are resisted only with reactions much larger than the load, and the reactions are sensitive to small geometric errors.
8. A body is in equilibrium if and only if the total work of all its forces is zero for every small rigid displacement. It is independent because it uses the measured change of the link lengths under an actual motion, not the equilibrium matrix that produced the solution.
