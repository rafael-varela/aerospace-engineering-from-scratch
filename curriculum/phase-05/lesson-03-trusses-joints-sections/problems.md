# P05.03  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer, say which members are in tension and which in compression, and check determinacy before solving.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a three-member truss by hand)

Joints at (0, 0), (4, 0) and the apex (1 + s/30, 3) m. Pin at the first joint, vertical roller at the second; members 0–1, 0–2, 1–2; a load of 1000 N downward at the apex.

(a) Find the reactions by taking moments about the pin.
(b) Find the three member forces by the method of joints, starting at the roller.
(c) Check with `solve_joints` and `joint_residual`.

Answer form: two reactions and three member forces with their sense.

## 2. (Derivation: rules and the beam analogy)

(a) Prove the two zero-force-member rules from the equilibrium of one joint.
(b) Show by an example that m + r = 2j is necessary but not sufficient for a planar truss to be determinate, and explain what the rank of the joint matrix adds.
(c) For a Pratt girder of depth h, show by the method of sections that the bottom chord force is M/h at the right joint, the top chord force −M/h at another, and the diagonal force V/sin θ. Say which joint each moment is taken about.

Answer form: two proofs, an example, a derivation.

## 3. (Unfamiliar scenario: a Warren-girder wing rib)

A Warren girder (no verticals) has five lower joints at x = 0, 0.5, 1.0, 1.5, 2.0 m (y = 0) and four upper joints at x = 0.25, 0.75, 1.25, 1.75 m (y = 0.35 m). Members: the lower chord (4), the upper chord (3), and the diagonals zig-zagging from each lower joint to the upper joint to its right and from each upper joint to the lower joint to its right (8). Pin at x = 0, vertical roller at x = 2.0 m. Each upper joint carries P = 300 + 10s N downward.

(a) Count: is the truss determinate? Confirm with `determinacy`.
(b) Build it as a `Truss` and solve it. Which member is the most compressed, which the most tensioned?
(c) Cut through the lower chord member from x = 1.0 to 1.5 m, the upper chord member from x = 0.75 to 1.25 m and the diagonal from x = 1.0 (lower) to 1.25 m (upper). Solve the cut with `solve_section` and by hand. Explain the force in the diagonal.
(d) Compare the largest chord force with M_max/h.

Answer form: a count, a table of 15 forces, three forces with an explanation, a comparison.

## 4. (Numerical, uses your code: choose the girder depth)

A Pratt girder of six panels spans L = 4 + s/10 m with 2000 N at each interior lower joint.

(a) With `depth_sweep` on depth/span from 0.05 to 0.6 (221 points), find the depth that minimises the volume index and report the largest chord force there.
(b) At depth/span = 0.25, report the largest compression in the top chord, the largest tension, the largest diagonal force and the largest vertical compression.
(c) In 19.05 × 0.889 mm 4130 tube (E = 200 GPa, typical), compute the stresses of (b) and the Euler load of one top-chord bay. Does the top chord buckle?

Answer form: a depth and a force, four forces, stresses and a yes/no with the ratio.

## 5. (Numerical: the 3-D engine mount)

Use `engine_mount_truss()` with its defaults.

(a) Solve it. Which member is the most compressed and which the most tensioned? Compute their stresses in the 19.05 × 0.889 mm tube and the Euler load of the most compressed one.
(b) Add a ring diagonal between ring joints 4 and 6. Classify the truss, solve it by the stiffness method with uniform EA, and compare with (a). Repeat with the diagonal four times stiffer.
(c) Remove ring member 0 (joints 4–5), then instead strut 4 (joints 4–0). Classify each.

Answer form: two members with forces and stresses and a buckling load, a class and a comparison, two classes.

## 6. (Judgment: when is a welded tube structure a truss?)

Answer each in two or three sentences.

(a) A welded 4130 fuselage truss has rigid joints. Under what conditions is the pin-jointed analysis still a good estimate of the member forces, and what does it miss?
(b) An engine bearer also supports an oil cooler at its mid-length. What changes in the analysis of that member and of the truss?
(c) Your truss solver says all members are within yield with a factor of 2.5. Name the check you have not done and the member most likely to fail it.
(d) A colleague adds a diagonal to every panel "for safety". What does that do to the analysis, and can it make some member worse off?

Answer form: four short paragraphs.
