# P05.03  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` (and the probe script that calls it) and given to the figures shown.

## 1. A three-member truss (s = 15: apex at (1.5, 3) m)

**Assumptions.** Ideal truss: pins, loads at joints, weight neglected.

(a) Moments about the pin: 4 R₁ = 1.5 × 1000, **R₁ = 375 N** up at the roller; vertical sum, **R₀ = 625 N** up at the pin; no horizontal load, so the pin's horizontal reaction is 0.

(b) Joint 1 (roller). Member 1–2 runs toward the apex along (−2.5, 3)/3.905 = (−0.6402, 0.7682). Vertical: 375 + 0.7682 T₁₂ = 0, **T₁₂ = −488.1 N** (compression). Horizontal: the base member pulls joint 1 toward joint 0, −T₀₁, so −T₀₁ − 0.6402 T₁₂ = 0, **T₀₁ = 312.5 N** (tension). Joint 0: member 0–2 runs along (1.5, 3)/3.354 = (0.4472, 0.8944); vertical: 625 + 0.8944 T₀₂ = 0, **T₀₂ = −698.8 N** (compression).

(c) `solve_joints` returns (312.5, −698.77, −488.14) N and reactions (0, 625) and (0, 375) N; the joint residual is zero to round-off.

**Check.** Joint 2: the two rafters' vertical components, 0.8944 × 698.8 + 0.7682 × 488.1 = 625.0 + 375.0 = 1000 N, hold the load.

**Common mistakes.** Taking the rafter force on the roller joint with the wrong direction vector (from the apex instead of toward it), which flips the sign of every force.

## 2. Rules and the beam analogy

(a) **Two non-collinear members, unloaded joint.** Resolve forces perpendicular to member A: only member B has a component there, T_B sin α = 0 with α ≠ 0, so T_B = 0; then T_A = 0 too. **Three members, two collinear, unloaded joint.** Resolve perpendicular to the collinear pair: only the third member has a component, so it carries zero.

(b) Take the four-panel Pratt girder and move the diagonal of panel 2 into panel 3 (which now has two). m and r are unchanged, so m + r = 2j still. But panel 2 is now a four-bar mechanism (rank drops by one) and panel 3 has a self-stress (two diagonals that can fight each other). The rank of B detects both: d·j − ρ = 1 mechanism, m + r − ρ = 1 self-stress. Counting cannot.

(c) Cut panel k; keep the left part. Its external loads and reactions give the beam's M and V at the cut. The cut members are the bottom chord (horizontal, lever h about any top joint), the top chord (horizontal, lever h about any bottom joint) and the diagonal. **About the top end of the diagonal**, the top chord and the diagonal pass through the point: T_bottom × h = M there (the bottom chord pulls left on the free body below the moment centre), **T_bottom = M/h**. **About the bottom end of the diagonal**: T_top = −M/h at that point. **Vertical sum**: the chords are horizontal, so T_diag sin θ = ±V, the sign set by the diagonal's slope.

**Common mistakes.** Using the same M for both chords (the moments are taken at two different joints); missing that a diagonal's sign depends on which way it slopes.

## 3. A Warren-girder wing rib (s = 15: P = 450 N at each upper joint)

(a) j = 9, m = 15, r = 3: m + r = 18 = 2j; `determinacy` gives rank 18: **determinate**.

(b) Reactions 900 N at each support. Member forces (N):

| member | joints | force | member | joints | force |
|---|---|---|---|---|---|
| lower 0–1 | 0–1 | 642.9 | diag | 0–5 | −1 106.0 |
| lower 1–2 | 1–2 | 1 285.7 | diag | 5–1 | 553.0 |
| lower 2–3 | 2–3 | 1 285.7 | diag | 1–6 | −553.0 |
| lower 3–4 | 3–4 | 642.9 | diag | 6–2 | 0.0 |
| upper | 5–6 | −964.3 | diag | 2–7 | 0.0 |
| upper | 6–7 | −1 285.7 | diag | 7–3 | −553.0 |
| upper | 7–8 | −964.3 | diag | 3–8 | 553.0 |
| | | | diag | 8–4 | −1 106.0 |

Most compressed: the middle upper chord 6–7, **−1 285.7 N**; most tensioned: the middle lower chords, **1 285.7 N**. The end diagonals (−1 106.0 N) carry the reactions up into the girder.

(c) The section gives **lower 1 285.7 N, upper −1 285.7 N, diagonal 0**, equal to the joint solution. By hand: the left part carries the reaction 900 N up and two loads of 450 N down, so the shear at the cut is **V = 0**: the diagonal, the only member with a vertical component, carries nothing. Symmetric loading makes the middle diagonals zero-force members, although no joint rule identifies them.

(d) M at the lower joint x = 1.0 m: 900 × 1.0 − 450 × 0.75 − 450 × 0.25 = 450 N·m; M/h = 450/0.35 = **1 285.7 N**, equal to the largest chord force.

**Common mistakes.** Taking moments about a joint whose lever to the chord is not h (the upper joints of a Warren girder are midway between lower ones, but the lever is still the depth); forgetting the end diagonals carry the reactions in compression.

## 4. Choose the girder depth (s = 15: span 5.5 m)

(a) The volume index is smallest at **depth/span = 0.26** (depth 1.43 m), 103 189 N·m, where the largest chord force is **5 769.2 N**.

(b) At depth/span = 0.25 (h = 1.375 m): top chord **−6 000.0 N**, bottom chord **5 333.3 N**, diagonal **6 009.3 N**, vertical **−5 000.0 N**. In a deep girder the web members carry forces as large as the chords.

(c) Stresses in 50.72 mm²: top chord **−118.3 MPa**, largest tension (the end diagonal) **118.5 MPa**. A bay is 5.5/6 = 0.917 m long; its Euler load is **4 924 N**, and the top chord carries **1.22 times** that: it buckles, at a stress well below yield. The volume-optimal girder is not automatically a buildable one: the top chord needs lateral bracing at mid-bay or a larger tube.

**Check.** 6 000 N = M_max/h: the reactions are 5 000 N each and M_max = 5 000 × 2.75 − 2 000 × (1.833 + 0.917) = 8 250 N·m, and 8 250/1.375 = 6 000 N.

**Common mistakes.** Sizing compression members for yield only; reading the volume index as the mass (it must be divided by an allowable stress and multiplied by density).

## 5. The 3-D engine mount

(a) Members are numbered: 0–3 the ring sides, 4–7 the struts from each ring joint to the corner behind it, 8–11 the struts to the next corner. Most compressed: **member 7** (ring joint 7 to firewall joint 3), **−4 001.9 N**, −78.9 MPa; most tensioned: **member 5** (ring 5 to firewall 1), **4 128.9 N**, +81.4 MPa. Member 7 is 0.559 m long; its Euler load is **13 249 N**, 3.31 times its force.

(b) With a ring diagonal 4–6 the truss has m + r = 25 > 24 and rank 24: **indeterminate of degree 1**. The stiffness method (uniform EA) puts −22.0 N in the diagonal and changes the other forces by less than 30 N (member 5: 4 146.3 N, member 7: −3 984.5 N); the joint residual is 2 × 10⁻¹². With the diagonal four times stiffer it carries −23.8 N. The ring is already well braced by the struts, so the diagonal attracts little load.

(c) Removing ring member 0 leaves a **mechanism** (one free motion); removing strut 4 also gives a **mechanism**. Every member of a determinate truss is essential.

**Check.** The sum of loads and reactions is zero (printed by `main`), and the stiffness method agrees with joints to 1.4 × 10⁻¹² N in (a).

**Common mistakes.** Forgetting that the ring loads come from a rigid-body transfer of the engine load (`rigid_body_nodal_loads`), not from equal shares.

## 6. When is a welded tube structure a truss?

(a) When the members are slender (length many times their diameter), the joints' centre lines meet at one point and the loads are applied at the joints, the axial forces dominate and the pin-jointed forces are good estimates. The analysis misses the secondary bending caused by rigid joints rotating as the truss deflects, which adds stress near the welds and matters for fatigue.

(b) The bearer becomes a beam-column: it carries bending from the cooler's weight (and its inertia under load factor) on top of its axial force, and its end forces are no longer along its axis. For the rest of the truss, replace the cooler by two statically equivalent joint loads at the bearer's ends; for the bearer itself, add the bending stress and check it as a beam-column (`P05.10`, `P05.15`).

(c) Buckling. The longest member in compression (here the top chord bays or the long compressed struts) is the most likely to fail, at a stress far below yield.

(d) The truss becomes indeterminate: forces now depend on member stiffnesses, fit-up and temperature. It can make some members worse: a stiff added diagonal attracts load and can overload the joints it connects, and a diagonal in compression can buckle. Redundancy is useful (a second load path), but it must be analysed with the stiffness method, not assumed safe.

## Check-yourself answers

1. Pinned at both ends and loaded only at its ends, so its two end forces must be equal, opposite and along the member. A load between the joints, or rigid joints, would break it.
2. Not necessarily: m + r = 18 = 2j is necessary, but a badly arranged truss can have a mechanism and a redundancy at the same time; the rank of the joint matrix decides.
3. Two non-collinear members at an unloaded, unsupported joint: both zero. Three members at an unloaded, unsupported joint with two collinear: the third is zero.
4. The free body on one side of the cut includes the supports on that side, so their reactions must be known; they come from the three equilibrium equations of the whole truss.
5. The chords carry the moment as a couple (T = ±M/h); the diagonals carry the shear (T sin θ = V).
6. Equilibrium alone fixes every member force of a determinate truss, so any solution that satisfies equilibrium, including the stiffness method's, has the same forces; EA changes only the displacements.
7. It can buckle: a slender compression member fails at the Euler load π²EI/L², which can be far below its yield load.
8. Σ|T|L, divided by the allowable stress, is the material volume of a truss with every member fully stressed. Shallow girders need heavy chords (forces ∝ 1/h); deep ones need long web members; the sum has a minimum in between.
