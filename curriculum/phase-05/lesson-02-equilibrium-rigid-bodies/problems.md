# P05.02  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and draw every support with its reaction components before writing any equation.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: classify, then solve)

A rigid beam runs from x = 0 to x = 4 m along the x axis and carries a downward load P = 1000 + 100s N at x = 1.5 m. Six ways of supporting it:

(a) a pin at 0 and a vertical roller at 4 m; (b) a fixed support at 0; (c) pins at 0 and at 4 m; (d) vertical rollers at 0, 2 and 4 m; (e) vertical rollers at 0 and 4 m and a horizontal roller at 2 m; (f) a pin at 0 and a link from x = 4 m to the fixed point (6, 0) m.

For each, predict the class (determinate, indeterminate with its degree, partially or improperly constrained) by counting and by looking for concurrent or parallel reaction lines; then confirm with `classify_supports`. Solve the determinate ones by hand and with `solve_reactions`.

Answer form: a six-row table and the reactions of the determinate cases.

## 2. (Derivation: what the rank counts)

(a) Show that three links whose lines meet at one point, in a plane, give an equilibrium matrix of rank 2, and describe the free motion and the self-stress state.
(b) Using the rank–nullity theorem, show that n − ρ is the number of independent self-stress states and n_eq − ρ the number of independent loads that cannot be balanced. Conclude that "determinate" requires ρ = n = n_eq.
(c) Derive the strut force (3) of the lesson, T = (L ȳ − n W s/2) ℓ_s/(y_s h), and show that the horizontal component of the root pin force does not depend on the strut station y_s.

Answer form: three short derivations.

## 3. (Unfamiliar scenario: a cargo door held open by a gas strut)

Axes x forward, y outboard, z up. A cargo door is hinged along the x axis at the top of the fuselage opening and held open horizontally, sticking outboard. Its hinges are at x = 0.1 m (a ball-and-socket) and x = 1.1 m (a hinge that lets the door slide along x). A gas strut runs from the door at (0.6, 0.5, 0) m down to the fuselage at (0.6, 0, −0.6) m. The door has mass 25 + s kg with its centre of gravity at (0.6, 0.45, 0) m, and a gust pushes down on it with a resultant of 150 N at (0.6, 0.6, 0) m.

(a) Classify the support set.
(b) Find the strut force by one moment equation about the hinge line. Is the strut in tension or compression?
(c) Find all six reactions with `solve_reactions` and check the residual.
(d) Why does neither hinge carry a force along x? What load would change that?

Answer form: a class, a force with its sense, six reactions, two sentences.

## 4. (Numerical, uses your code: size the strut station)

An aircraft of mass m = 900 + 20s kg has the strut-braced wing of the lesson (semispan 5 m, strut from (0, −1.0) m, 60 kg half-wing) at n = 3.8.

(a) Compute the strut force and the root pin reactions for y_s = 2.6 m, and check the strut force with `strut_force_closed_form`.
(b) Sweep y_s from 1.0 to 4.6 m (361 points) and find the station that minimises the root pin resultant. Explain why the minimum is where it is.
(c) At n = −1.52 (lift −0.4 times the positive value, inertia of the wing reversed), what does the strut carry at y_s = 2.6 m?
(d) The strut force falls as y_s grows. Give two reasons why designers do not simply attach it at the tip.

Answer form: three forces, a station and a sentence, a force with its sense, two reasons.

## 5. (Numerical and judgment: the engine mount near an improper constraint)

(a) Apply only the propeller torque, 450 N·m about x, to `hexapod_mount(spread_deg=...)` for spreads of 40°, 20°, 10°, 5° and 2°. Tabulate the largest link force, its product with sin(spread) and the condition number. Explain the trend.
(b) Add the seventh link (`extra_link=True`) and apply the lesson's full load. `solve_reactions` refuses; `numpy.linalg.lstsq` returns the minimum-norm solution. Show that it balances the load, find the self-stress state (the null vector of A), and explain why the minimum-norm answer is not what the real mount carries.
(c) Show that removing any one of the six links leaves the mount partially constrained.

Answer form: a five-row table and a sentence, a vector and a paragraph, a one-line result.

## 6. (Judgment: determinacy in aerospace practice)

Answer each in two or three sentences.

(a) Why are spacecraft instruments (telescopes, optical benches) usually mounted on exactly six links or three bipods?
(b) A fail-safe engine mount deliberately includes a redundant link. What does it buy, and what analysis does it require that a determinate mount does not?
(c) A beam on three vertical rollers is "improperly constrained", yet it carries vertical loads perfectly well. In what sense is the classification still important?
(d) Your mount's condition number is 300. Name one consequence for the link forces and one for the sensitivity of those forces to manufacturing tolerances.

Answer form: four short paragraphs.
