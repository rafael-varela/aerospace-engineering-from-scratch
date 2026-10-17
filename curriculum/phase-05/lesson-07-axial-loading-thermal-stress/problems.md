# P05.07  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical. Material values are the typical values of `MATERIALS` in `starter/main.py`.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a drag-brace tube)

A landing-gear drag brace is a 4130 steel tube, outside diameter 25.4 mm, wall 1.6 mm, length 600 mm, carrying a compressive force of 20 + s kN. Ignore buckling (that is `P05.15`).

(a) Find the stress, the shortening and the ratio of the stress magnitude to the typical yield strength.
(b) The same tube, unloaded, is installed between two rigid fittings and heated by 30 K. What stress does it carry?
(c) What temperature rise alone would bring the restrained tube to its typical yield strength?

Answer form: three numbers with units, one stress, one temperature difference.

## 2. (Derivation: four closed forms of the lesson)

(a) Integrate the elongation of a bar whose area varies linearly from A₀ to A₁ and obtain (2). Show that it tends to PL/(EA₀) as A₁ → A₀.
(b) Derive (5) for the fixed–fixed two-segment bar from equilibrium of the junction, compatibility of the two segments and Hooke's law with thermal strain (3). Show that load and temperature superpose.
(c) Derive the frustum stiffness (7) by partial fractions.
(d) From (8) and (9), derive the bolt and clamp forces (10), the separation load (11), and the thermal preload change (12). Explain why the bolt and the members act in parallel for an external load but in series for a thermal misfit.

Answer form: four derivations and one paragraph.

## 3. (Unfamiliar scenario: a liquid-oxygen feed line)

A 2024-T3 propellant line, outside diameter 50.8 mm, wall 1.25 mm, runs straight between two rigid mounts L = 1.5 + s/30 m apart. It is installed at 20 °C and filled with liquid oxygen at −183 °C. Take α and E as constant (they are not; see (d)).

(a) Without any flexibility, what axial stress and force does the line carry when filled? Compare with the typical yield strength.
(b) A bellows (an axial spring of stiffness k_bel) is inserted in series with the line. Find the largest k_bel for which the line stress stays below one third of the typical yield strength, and the stroke the bellows must absorb.
(c) Check your answer to (b) with `bar_fe_solve`: one element for the line and one for the bellows between two fixed nodes, with the thermal strain as an initial strain on the line only.
(d) Name two reasons the constant-property answer to (a) is inaccurate at −183 °C, and say which way each one moves it.

Answer form: a stress, a force and a comparison; a stiffness and a stroke; an FE stress; two sentences.

## 4. (Numerical, uses your code: can torque control meet the requirements?)

Use the lesson joint, its load case and its temperatures. Draw 100 000 nut factors from N(0.20, 0.025) with `numpy.random.default_rng(s)`.

(a) For tightening torques of 28.8, 32.0, 35.2 and 38.4 N·m (the nominal torques for 18, 20, 22 and 24 kN at K = 0.2), compute the fraction of joints that separate below the 18 kN ultimate load when cold and the fraction whose bolt stress exceeds the typical yield at the 12 kN limit load when hot.
(b) Derive the deterministic preload window: the smallest assembly preload that keeps the cold joint clamped at ultimate, and the largest that keeps the hot bolt below yield at limit. Explain the shape of the trade-off in (a) with it.
(c) Repeat (a) at 32.0 N·m with a nut-factor standard deviation of 0.01, and report the 1st–99th percentile preload range. What preload method could deliver that scatter?

Answer form: a 4 × 2 table and a sentence, two preloads, two fractions, a range and a sentence.

## 5. (Numerical, uses your code: a rotating compressor blade)

Treat a Ti-6Al-4V compressor blade as a uniform bar of length 100 mm, root at radius r₀ = 250 mm, spinning at 12 000 + 200s rpm. Its density is 4430 kg/m³. The centrifugal body force per unit length is q(r) = ρAω²r.

(a) By the method of sections, show that the root stress is σ = ρω²(r₀L + L²/2), independent of the area, and evaluate it. Integrate the strain to find the tip stretch δ = ρω²(r₀L²/2 + L³/3)/E.
(b) With `bar_fe_solve` and consistent nodal loads for a linearly varying q (for an element of length h between loads q_a and q_b: h(2q_a + q_b)/6 to the first node and h(q_a + 2q_b)/6 to the second), compute the tip stretch and the root-element stress with 2, 4, 8, 16 and 32 elements. What do you notice about the tip stretch, and at what order does the root-element stress converge?
(c) The blade also runs 150 K hotter than at assembly. What is its free thermal growth, and what is the total change of tip radius that the casing clearance must allow for (ignore the disk)?

Answer form: a derivation and two numbers; a 5-row table and two sentences; two numbers.

## 6. (Judgment: joint design decisions)

Answer each in two or three sentences, with a number where one is asked for.

(a) Why are preloads typically set at a large fraction (often more than half) of the bolt's yield strength rather than "just enough"?
(b) For the lesson joint, compute C for a bolt with an all-shank grip (no thread in the grip). Is a stiffer bolt better or worse for fatigue of the bolt, and for preload loss in the cold?
(c) The pressure-cone model assumes that the plates extend well beyond the cone. Where on an aircraft fitting does this fail, and which way does it move C?
(d) Name one way to measure the preload in a bolt after installation, and say what it would tell you about the nut factor.

Answer form: four short paragraphs, one number.
