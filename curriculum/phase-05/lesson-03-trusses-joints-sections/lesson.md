# P05.03  Trusses: Joints & Sections

> A truss turns bending into tension and compression; solve it three ways and make them agree.

| | |
|---|---|
| **Learning outcomes** | 1. Explain why every member of an ideal pin-jointed truss is a two-force member, and identify zero-force members by inspection. 2. Assemble the method of joints as one linear system B x = −f for any 2-D or 3-D truss, and use counting (m + r against d·j) and rank to decide whether a truss is determinate, indeterminate or a mechanism. 3. Find the force in any member of a planar truss by the method of sections, with the reactions taken from equilibrium of the whole truss, and relate the chord and diagonal forces of a girder to the bending moment and shear of an equivalent beam. 4. Implement the direct stiffness method and show that it reproduces the joint forces of a determinate truss whatever the member stiffnesses, and that it alone can solve an indeterminate one. 5. Use the solver to choose the depth of a girder and to check a welded-tube engine mount, and say when a compression member is governed by buckling rather than by yield. |
| **Prerequisites** | `P05.01` (two-force members, moments), `P05.02` (equilibrium matrix, determinacy by rank), `P03.08` (direct linear solvers), `P03.13` (finite elements: assembly of a stiffness matrix). Skill check: you can assemble a sparse pattern of 2 × 2 blocks into a larger matrix with NumPy slicing. |
| **Engineering question** | A 4.8 m Pratt girder of 4130 steel tube, 0.6 m deep, carries 2 kN at each of its five lower joints. Which members are in tension, which in compression, and how hard is each loaded; what depth would use the least material; and would the most compressed tube yield or buckle first? And for a 3-D welded-tube engine mount: is it determinate, and what does each tube carry under a 4 g landing with full thrust? |
| **Assumptions and limits** | Members straight, pinned at both ends and loaded only at the joints (no member weight, no loads between joints); small displacements; linear elastic material for the stiffness method. Real welded trusses have rigid joints and carry secondary bending; members loaded between joints are beams, not two-force members. |
| **Analytical work** | The two-force member; joint equilibrium; the joint matrix; counting and rank; zero-force member rules; external reactions from the whole truss; the method of sections; chord force = M/h and diagonal force = V/sin θ; member stiffness EA/L and the assembled stiffness matrix; why EA drops out of a determinate truss; a first look at Euler buckling of a compression member. |
| **Computational deliverable** | `starter/main.py`: `member_geometry`, `joint_matrix`, `determinacy`, `solve_joints`, `joint_residual`, `zero_force_members`, `external_reactions_2d`, `solve_section`, `stiffness_solve`, `beam_moment_shear`, `pratt_beam_analogy`, `tube_area`, `tube_second_moment`, `euler_load_pinned`, `volume_index`, `depth_sweep` (the `Truss` dataclass, `TrussError`, the generators `pratt_truss`, `pratt_members`, `pratt_section`, `engine_mount_truss`, `rigid_body_nodal_loads` and `main` are given). Tests prove: member columns are self-equilibrated; the classification of determinate, indeterminate and mechanism trusses; a hand-solved triangle; zero residual at every joint and global balance; symmetry, superposition and the zero-load limit; zero-force members; reactions from the whole body equal those from the joints; sections equal joints in every panel from either side; the beam analogy and the 1/depth scaling; the stiffness method equals joints for any EA and shares load by stiffness when indeterminate; the 3-D mount is determinate and its loads statically equivalent. |
| **Practical activity** | **computational experiment**: a sweep of girder depth against chord force and material volume, a method comparison (joints, sections, beam analogy, stiffness), and a check of the 3-D engine mount. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `truss.png`, a hand solution of one Pratt panel by sections, the depth-sweep table, the method-comparison table, and a half-page note answering the engineering question. |
| **Assessment** | `assessments/phase-05/grading-rubric.md`; competencies 1–4. Pass: the determinacy of each truss is argued by counting and by geometry before the code is run; one panel is solved by sections by hand and matches the code; all public tests pass; the buckling-versus-yield conclusion is stated with numbers. |
| **References** | `@libretexts-baker-haynes-engineering-statics` Ch. 6 "Equilibrium of Structures": 6.1 Structures, 6.2 Interactions between members, 6.3 Trusses, 6.4 Method of Joints, 6.5 Method of Sections (required-free: the two classical methods); `@mit-ocw-3-11-mechanics-of-materials` Roylance module "Trusses": Statics analysis of forces, Statically indeterminate trusses, Matrix analysis of trusses, Stiffness matrix for a single truss element, Assembly of multiple element contributions (required-free: the bridge from joints to the stiffness method); `@mit-ocw-16-001-unified-materials-structures` lecture notes "Truss Analysis" and "Statically Indeterminate Systems", and Lab 01 §3 "Truss Lab" (required-free: an aerospace department's treatment and a measured truss). |
| **Extensions** | ★ Add member self-weight as half-loads at each end joint and see how much it changes the Pratt forces. ★ Solve the 3-D engine mount by sections: a cut through the six struts and the moment about a well-chosen axis. ★ Compute joint deflections with the stiffness method for the girder and compare with the beam formula 5wL⁴/(384EI) using I ≈ 2A_chord(h/2)² (a preview of `P05.12`). |
| **Time** | unmeasured |

## 1. The engineering problem

A beam resists a load by bending: one side stretches, the other shortens, and the material near the middle does little. A **truss** removes that idle middle. Straight members meet at joints and form triangles; loaded only at the joints, each member is simply stretched or squeezed. That is why trusses appear wherever weight matters and depth is available: the welded 4130 steel-tube fuselages of light aircraft and helicopters, engine mounts, the spars of early aircraft, launch-vehicle interstages and the frames of space structures.

The design questions are always the same. Which members pull and which push? How hard? Is the truss stable at all, and can equilibrium alone answer the first two questions? This lesson solves trusses three independent ways: the **method of joints** written as one linear system; the **method of sections**, with the reactions taken from the whole structure; and the **direct stiffness method** of `P03.13`. It applies them to a 4.8 m **Pratt girder** carrying 2 kN at each of five lower joints, whose most compressed tube carries 12.0 kN, and to a 3-D **engine-mount truss**.

Success criteria, written first: the three methods agree member by member to round-off; every joint is in equilibrium; the girder's chord forces equal the bending moment divided by the depth; the conclusions hold when the member stiffnesses are changed at random.

## 2. Understand and derive

### 2.1 The ideal truss

The idealisation has three parts (Baker & Haynes 6.3): members are straight; they are joined by frictionless pins at their ends; and loads act only at the joints. Then each member has forces only at its two ends, and by the argument of `P05.01` those forces are equal, opposite and along the member: **every member is a two-force member** with one unknown, its axial force T. Take T > 0 for **tension**: the member pulls on the joints at its ends. Compression (T < 0) pushes them apart.

The convention matters at the joints. If member k joins joint i to joint j with unit vector u_k pointing from i to j, its tension pulls joint i toward j, a force +T_k u_k on joint i, and pulls joint j toward i, −T_k u_k on joint j.

### 2.2 The method of joints as one linear system

Each joint is a particle: its free-body diagram has the member forces, the applied load and any support reactions, and its equilibrium is d equations (d = 2 in a plane, 3 in space) with no moment equation, because all the forces pass through the joint. For j joints that is d·j equations in the m member forces and the r reaction components. Collect them:

**B x = −f,**   (1)

where x = (T₁, …, T_m, R₁, …, R_r), f stacks the applied joint loads, and B is the **joint matrix**: the column of member k has +u_k in the rows of joint i and −u_k in the rows of joint j (and zeros elsewhere); the column of a reaction component has a single 1 in its row. Each member column sums to zero force: a member is internal to the truss and cannot change its overall equilibrium. The textbook method of joints (6.4), visiting joints one at a time, is Gaussian elimination of (1) done by hand in a lucky order.

### 2.3 Is the truss determinate?

The counting rule says a truss can be determinate only if

m + r = d·j.   (2)

It is necessary but not sufficient. As in `P05.02`, the rank ρ of B decides: **m + r − ρ** is the number of self-stress states (the degree of indeterminacy), and **d·j − ρ** the number of mechanisms. Two examples. A Pratt girder of n panels with vertical end posts has j = 2n + 2 joints, m = 4n + 1 members and r = 3 reactions: m + r = 2j, and the rank is full: **determinate**. Add the second diagonal in every panel and m + r = 2j + n: indeterminate of degree n. Remove one diagonal and m + r = 2j − 1: the panel without a diagonal is a four-bar mechanism that shears. The rank also catches the cases counting misses: three collinear joints joined by two members, the end joints pinned and the middle one free, satisfy the count (m + r = 6 = 2j) but cannot resist a sideways load on the middle joint.

### 2.4 Zero-force members

Two rules find unloaded members by inspection, from the equilibrium of a single joint:

- If **two non-collinear members** meet at an **unloaded** joint with no support, both carry zero force.
- If **three members** meet at an unloaded, unsupported joint and **two are collinear**, the third carries zero force.

Zero-force members are not useless: they brace compression members against buckling and carry load when the loading changes. In the lesson's Pratt girder the two end panels' bottom chords and the middle vertical carry nothing under the symmetric joint loads.

### 2.5 The method of sections

To find one member force without solving the whole truss (6.5): first find the external reactions from the equilibrium of the **whole truss as one rigid body** (three equations in 2-D); then cut through the truss across at most three members whose lines do not all meet at one point, and write the three equilibrium equations of the part on one side of the cut. The cut members' forces appear as external forces on that part. Take moments about the point where two of the cut members meet, and the third member's force comes out of a single equation.

This is a genuinely independent route: it uses the whole-body reactions and one free body, not the joint matrix.

### 2.6 The truss as a beam

Cut a Pratt girder of depth h through a panel. The cut crosses the top chord, the bottom chord and one diagonal. Treat the girder as a simply supported beam with the same loads: at the cut it has a bending moment M(x) and a shear V(x) (both computed from the loads and reactions on one side). The free body to the left of the cut must reproduce them:

- Moments about the joint where the diagonal meets the **top** chord: only the bottom chord has a lever there, h, so **T_bottom = M/h** (tension under a sagging moment).
- Moments about the joint where the diagonal meets the **bottom** chord: **T_top = −M/h** (compression).
- Vertical forces: only the diagonal has a vertical component, so **T_diag sin θ = V**, with θ the diagonal's angle to the horizontal.   (3)

So the chords carry the bending moment as a couple, and the diagonals carry the shear. The lesson's girder: span 4.8 m, six panels of 0.8 m, depth 0.6 m, 2 kN at each of the five interior lower joints. Each reaction is 5 000 N. At x = 0.8 m, M = 5 000 × 0.8 = 4 000 N·m, so the top chord of the first panel carries −4 000/0.6 = **−6 667 N**; the diagonal there (length 1.0 m, sin θ = 0.6) carries V/sin θ = 5 000/0.6 = **8 333 N** in tension. At mid-span, M = 7 200 N·m and the top chord carries **−12 000 N**, the largest compression in the girder.

Equation (3) has a design consequence: the chord forces are **inversely proportional to the depth**. A deeper girder has smaller chords but longer diagonals and verticals. The material needed, if every member is sized to the same allowable stress, is proportional to the **volume index** Σ|T_k| L_k, and that has a minimum at some depth: the experiment of §4 finds it.

### 2.7 The direct stiffness method

Give member k an axial stiffness EA/L: stretched by δ, it carries T = (EA/L)δ. Let q_i be the displacement of joint i. Small displacements stretch member k by δ_k = u_k · (q_j − q_i), and the forces it puts on its two joints are ∓T_k u_k. In matrix form, member k adds to the global stiffness matrix the block

k_k = (EA/L)_k [[u uᵀ, −u uᵀ], [−u uᵀ, u uᵀ]]   (4)

in the rows and columns of joints i and j (Roylance, "Stiffness matrix for a single truss element", "Assembly of multiple element contributions"; `P03.13`). Delete the rows and columns of the supported components, solve K_ff q_f = f_f for the free displacements, and recover each force as T_k = (EA/L)_k δ_k.

For a **determinate** truss, equilibrium (1) already fixes every T_k, so the stiffness method must return the same forces **whatever** the stiffnesses: only the displacements change. That is a powerful test. For an **indeterminate** truss, (1) has infinitely many solutions, and the stiffness method picks the one whose member stretches are **compatible**, i.e. fit together at the joints. There the forces depend on the relative stiffnesses: a stiffer member attracts more load.

### 2.8 Material, stress and a first look at buckling

The girder is built from 4130 steel tube, 19.05 mm (3/4 in) outside diameter with a 0.889 mm (0.035 in) wall, a common size in welded aircraft structure. Its area is π/4 (D² − d²) = **50.72 mm²** and its second moment of area π/64 (D⁴ − d⁴) = **2 096 mm⁴** (`P05.05`). Typical properties of normalized 4130 sheet and tube are E ≈ 200 GPa, a tensile yield strength F_ty ≈ 517 MPa and an ultimate strength F_tu ≈ 655 MPa (the minimum design values long published for this alloy in the US metallic-materials handbook MIL-HDBK-5, now MMPDS; quoted here as typical values for illustration, not as design allowables).

The most loaded tension member, the middle bottom chord at 10 667 N, is stressed to T/A = **210 MPa**, well under yield. The top chord at −12 000 N is stressed to 237 MPa in compression, also under yield, but a slender tube in compression fails by **buckling** first. The Euler load of a pin-ended strut, π²EI/L² (`P05.15`), is **6 465 N** for one 0.8 m bay of this tube: the top chord carries 1.86 times its buckling load. It would buckle at 127 MPa, far below yield. The member forces in this lesson are correct; the tube size is not, and only a compression check finds that.

### 2.9 Where the model stops being valid

- **Welded joints are not pins.** A welded tube truss has rigid joints; the members bend a little as the truss deflects. For slender members (length/diameter above about 10) these secondary moments are small; for stubby members they are not.
- **Loads between joints.** An engine bearer that also supports equipment between its ends, or a fuselage longeron loaded by the floor, is a beam-column, not a two-force member.
- **Eccentric joints.** If the member centre lines do not meet at one point, each joint carries a moment.
- **Buckling and large deflections.** The method gives forces; whether a compression member survives is a separate check.

## 3. Build the computational model

Open `starter/main.py`. The `Truss` dataclass (nodes, members as index pairs, a boolean `fixed` array marking the supported displacement components, and joint loads), `TrussError`, the generators `pratt_truss`, `pratt_members`, `pratt_section` and `engine_mount_truss`, the helper `rigid_body_nodal_loads` and `main` are given.

### Step 1: geometry and the joint matrix

```python
def joint_matrix(t):
    n, d = t.nodes.shape
    m = len(t.members)
    fixed = np.argwhere(t.fixed)
    B = np.zeros((n * d, m + len(fixed)))
    _, U = member_geometry(t)
    for k, (i, j) in enumerate(t.members):
        B[i * d:(i + 1) * d, k] += U[k]
        B[j * d:(j + 1) * d, k] -= U[k]
    for c, (i, comp) in enumerate(fixed):
        B[i * d + comp, m + c] = 1.0
    return B
```

`member_geometry` returns the lengths and unit vectors from the first node to the second; it raises `TrussError` on a zero-length member.

### Step 2: determinacy and the joint solution

`determinacy` computes the singular values of B and the rank (threshold 10⁻¹⁰ of the largest), and returns the counts, degrees and a status. `solve_joints` raises `TrussError` unless the truss is determinate and then solves (1); it returns the member forces and an (n, d) array of reactions. `joint_residual` returns the net force on every joint; `zero_force_members` lists members below 10⁻¹⁰ of a typical force.

### Step 3: sections and the beam analogy

`external_reactions_2d` solves the three whole-body equations for exactly three reaction components. `solve_section(t, cut, side, reactions)` builds the 3 × k system of the free body made of the joints in `side`: each cut member pulls its `side` joint toward the joint across the cut. `beam_moment_shear(t, x)` sums the loads and reactions left of x; `pratt_beam_analogy` applies (3).

### Step 4: the stiffness method

`stiffness_solve(t, EA)` assembles (4) with slices, deletes the fixed components with a boolean mask, checks the condition number of K_ff (above 10¹⁴ it raises `TrussError`: a mechanism has a singular stiffness matrix), solves, and returns displacements, member forces and reactions (K q − f at the supported components).

### Step 5: sizing helpers and the sweep

`tube_area`, `tube_second_moment`, `euler_load_pinned`, `volume_index` and `depth_sweep` are one-liners from §2.6 and §2.8.

Running `main()` prints the Pratt girder by all three methods, the stresses and the buckling check, the crossed-diagonal variant, the depth sweep and the engine mount, and writes `truss.png` (the girder coloured by tension and compression with the forces in kN; the chord force and the volume index against depth/span).

### Verify

Run `python -m pytest curriculum/phase-05/lesson-03-trusses-joints-sections`. The tests check:

1. **Joint matrix**: shape (d·j, m + r); every member column sums to zero force and touches exactly two joints.
2. **Determinacy**: the Pratt girder is determinate; crossed diagonals make it indeterminate of degree n; a missing diagonal makes it a mechanism (refused by both solvers); collinear joints are a mechanism despite the count.
3. **Triangle**: the hand solution (reactions by moments, then two joints).
4. **Equilibrium**: zero net force at every joint and global balance, for two girders and the 3-D mount.
5. **Symmetry, superposition, zero load**; the signs of chords and diagonals under gravity loads.
6. **Zero-force members** in a small truss and the three in the girder.
7. **Reactions** from the whole body equal those from the joints, with an unsymmetric extra load.
8. **Sections** equal joints in every panel of two girders, from either side of the cut.
9. **Beam analogy** equals joints at three depths; the largest chord force scales exactly as 1/depth.
10. A member that does not cross the cut is refused.
11. **Stiffness method** equals joints for uniform and random EA, for the girder and the mount.
12. **Indeterminate truss**: equilibrium holds, the stiffer diagonal attracts load, the external reactions do not change, and uniform scaling of EA scales displacements only.
13. **Mount**: determinate, 12 members, 12 reactions; the nodal loads are statically equivalent to the engine load.
14. **Tube formulas** and their solid-bar limits.
15. `main` writes the figure.

### Compare

Three routes to the same numbers: joints (one linear solve of a 28 × 28 system for the girder), sections (a 3 × 3 solve per panel, with reactions from the whole truss) and stiffness (a 25 × 25 solve of displacements). For the girder, sections and the beam analogy agree with joints to the last printed digit in every panel, and the stiffness method agrees to 7 × 10⁻¹¹ N with uniform EA and 8 × 10⁻¹¹ N with EA varied at random by a factor of 25. A mistake in the sign convention of any one of them (pulling a joint the wrong way) breaks the agreement at once.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** The results below are computations with the pin-jointed model, not measurements. Write each prediction down before running it.

1. **Method comparison.** Solve the girder by joints, by sections (every panel), by the beam analogy and by the stiffness method with uniform and random EA. Predict which agree exactly and which only approximately.
2. **Depth sweep.** Vary the depth from 0.05 to 0.5 of the span. Predict the largest chord force from (3) and the shape of the volume index Σ|T|L against depth; guess where its minimum lies.
3. **Indeterminate girder.** Add the second diagonal in every panel. Predict how the two diagonals of the last panel share the shear with equal EA, and what happens when one is ten times stiffer.
4. **Engine mount.** Solve the 3-D mount (four firewall joints fixed, four ring joints; 120 kg engine at 4 g, 3 kN thrust, −450 N·m torque, CG 0.30 m ahead of the ring). Predict its determinacy by counting before running it.

Raw data format: CSV tables `member, joints_N, sections_N, beam_N, stiffness_N` (1); `depth_over_span, max_chord_N, volume_index_Nm` (2); `member, force_equalEA_N, force_stiff_N` (3); `member, force_N, stress_MPa` (4); plus `truss.png`.

What you should find (model results): (1) all four agree; the forces in kN are, bottom chord 0, 6.67, 10.67, 10.67, 6.67, 0; top chord −6.67, −10.67, −12.00, −12.00, −10.67, −6.67; verticals −5.00, −3.00, −1.00, 0, −1.00, −3.00, −5.00; diagonals 8.33, 5.00, 1.67, 1.67, 5.00, 8.33. (2) The largest chord force is 30.0 kN at depth/span 0.05, 15.0 kN at 0.1, 6.0 kN at 0.25 and 3.0 kN at 0.5, exactly 1/depth. The volume index falls from 243 307 N·m at 0.05 to a broad minimum near 90 000 N·m and rises again to 109 867 N·m at 0.5; a fine sweep puts the minimum at **depth/span = 0.26** (90 133 N·m at 0.25). Below it the chords dominate; above it the long web members do. (3) The indeterminate girder (degree 6) shares the last panel's shear between its two diagonals as 3 997 N and −4 336 N with equal EA; making one diagonal ten times stiffer changes this to 2 463 N and −5 871 N. The external reactions do not change. (4) m = 12, r = 12, 3j = 24, rank 24: determinate. The member forces range from −4 002 N to +4 129 N, i.e. −79 to +81 MPa in the same 19.05 × 0.889 mm tube; the most compressed strut is 0.559 m long with an Euler load of 13 249 N, so it is safe from buckling by a factor of 3.3. The stiffness method agrees to 1.4 × 10⁻¹² N.

## 5. Validate and improve

| quantity | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| reactions (N) | 5 000, 5 000 | 5 000.0, 5 000.0 | 0 | symmetry |
| top chord, panel 1 (N) | −4 000/0.6 = −6 667 | −6 666.7 | 0 | — |
| diagonal, panel 1 (N) | 5 000/0.6 = 8 333 | 8 333.3 | 0 | — |
| mid-span top chord (N) | −7 200/0.6 = −12 000 | −12 000.0 | 0 | — |
| sections vs joints, all panels (N) | equal | equal | below 10⁻⁹ | round-off |
| stiffness vs joints, random EA (N) | equal | 8 × 10⁻¹¹ | round-off | — |
| Euler load of a 0.8 m bay (N) | π²EI/L² = 6 465 | 6 465 | 0 | — |
| optimum depth/span (volume index) | — | 0.26 | — | sweep resolution 0.0025 |

The truss model is exact for its idealisation, and three independent methods confirm it. Its real weakness is the tube: the top chord would buckle. Two improvements follow. Use a larger tube for the top chord (a 22.2 × 0.889 mm (7/8 in) tube has 1.62 times the second moment of area, an Euler load of 10 476 N for the bay, still too little; a 25.4 × 0.889 mm (1 in) tube has 2.46 times, 15 877 N, enough without a safety factor and not with one) or shorten its buckling length with lateral bracing; and in a real welded design, check the secondary bending at the rigid joints.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the counting and rank argument for each truss; a hand solution of one Pratt panel by sections; `truss.png`; the four experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) answering the engineering question, including the buckling conclusion.

### Adversarial self-checks

- **What if the load doubles?** Every member force doubles; the stiffness method's displacements double too. Whether the truss survives depends on the stresses and the buckling loads, which do not scale the same way: the top chord was already beyond its Euler load.
- **What if the load moves?** A single load at one joint changes which diagonals are loaded and can reverse a diagonal from tension to compression; a Pratt girder is efficient only for the loading it was designed for. Problem 3 shows a Warren girder where the middle diagonals carry nothing under a symmetric load.
- **What if a member is infinitely stiff?** In a determinate truss nothing changes in the forces (the test with random EA proves it); in an indeterminate one the stiff member attracts load.
- **Which assumption fails first in the real world?** The pinned joint for stubby members, and the absence of loads between joints for longerons and engine bearers.
- **What measurement could falsify the model?** Strain gauges on a few members of a loaded truss, as in the MIT 16.001 truss lab: measured strains times EA should match the joint forces; a systematic difference in compression members points at bending from rigid joints or initial crookedness.
- **Does the result survive a change of method?** Yes: joints, sections, beam analogy and stiffness agree to round-off.

## Check yourself

1. Why is every member of an ideal truss a two-force member, and what would break that?
2. A planar truss has 9 joints, 15 members and 3 reaction components. Is it necessarily determinate?
3. State the two rules for zero-force members.
4. Why does the method of sections need the reactions first, and how are they found?
5. In a girder, which members carry the bending moment and which the shear? Write the chord force in terms of M and h.
6. Why does the stiffness method give the same member forces as the method of joints for a determinate truss, whatever EA is?
7. A top-chord tube carries less than its yield load. Why might it still fail?
8. What does the volume index Σ|T|L measure, and why does it have a minimum against depth?
