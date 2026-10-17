# P05.11  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical. Draw the free-body diagram (or the flow diagram) before writing any equation.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

Material data (typical, not design allowables): 7075-T6 E = 71.7 GPa, G = 26.9 GPa, ν = 0.33, ultimate shear strength F_su = 331 MPa.

## 1. (Warm-up: a rectangle and a pin)

A solid rectangular beam 25 mm wide and 60 mm deep carries a shear force V = 2000 + 100s N.

(a) Compute the average shear stress V/A, the maximum shear stress and the shear stress a quarter of the depth from the top.
(b) A 20 mm diameter solid round pin carries the same shear as a beam (not in double shear). Compute its average and maximum shear stress.
(c) Show that integrating τb over the depth of the rectangle returns V.

Answer form: five numbers with units, one short integration.

## 2. (Derivation: VQ/(It), the I-beam flows and the channel shear centre)

(a) Draw the free-body diagram of the part of a beam slice dz above a horizontal cut at y₁, with the bending stresses on its ends and the shear stress on the cut, and derive τ = VQ/(It). Which of the phase's three ingredients (equilibrium, compatibility, constitutive law) does the derivation use?
(b) For a thin-walled I-section (centreline depth h_c, flange width b and thickness t_f, web t_w) under S_y, derive the flow in a half-flange and in the web from (4) of the lesson. Show that the web flow integrates to S_y exactly and that the flanges carry no net vertical force.
(c) For a channel (web h, flanges b, thicknesses t_f and t_w), derive e = 3b²t_f/(6bt_f + ht_w). Find its limits as t_w → 0 and as b → 0, and explain each in a sentence.

Answer form: three derivations, two limits with sentences.

## 3. (Unfamiliar scenario: angles and tees)

Many aircraft brackets and stiffeners are angles (L-sections) and tees.

(a) Argue, without computing, where the shear centre of an angle and of a T-section must be. Then confirm with `shear_centre` for an equal angle with legs 40 mm, 3 mm thick, and for a T with a 60 × 4 mm flange and a 50 × 3 mm stem (centreline dimensions).
(b) An equal angle with legs a = 40 + s mm and 3 mm thick is used as a cantilever stiffener with a vertical 500 N load applied through its centroid. Compute the centroid, the second moments and the product of inertia about centroidal axes parallel to the legs, and the torque the load applies about the shear centre.
(c) Using the product of inertia you found, explain why the stiffener also bends sideways under a vertical load, and which bending formula from `P05.10` you would use for its stress.

Answer form: an argument and two confirmations; four section properties and a torque; one paragraph.

## 4. (Numerical, uses your code: the web and the joint of your own spar)

Take the aeroplane and spar of `P05.10` Problem 4 (mass 550 + 12s kg, n = 3.8, semi-span 4.5 + s/20 m, depth h = 0.14 + s/500 m, caps 10 mm thick, web 3 mm), with the cap width rounded to the millimetre from that problem (for s = 15, 46 mm). The root shear is the semi-span lift.

(a) Compute the web shear stress at the neutral axis and at the top of the web with `i_section_shear_stress`, the average V/(h_w t_w), and the fraction of V the web carries.
(b) Build the centreline model with `i_beam_segments(h − t_f, b_f, t_f, t_w)` and compare its peak web stress with (a).
(c) Compute the joint flow VQ_cap/I and the fastener pitch for two rows of fasteners at 4 kN each with the 1.5 factor.
(d) Compute the ultimate shear-strength margin of the web, MS = F_su/(1.5τ) − 1.

Answer form: three stresses and a fraction; a stress and a percentage difference; a flow and a pitch; a margin.

## 5. (Numerical, uses your code: a control-surface channel spar, with and without lips)

An aileron spar is a formed 2024-T3 channel of web depth h = 80 + s mm, flanges 30 mm, thickness 1.6 mm (centreline dimensions). The hinge brackets load it through its web.

(a) Compute the shear centre with `shear_centre` and with (6) of the lesson, and the centroid.
(b) For a hinge reaction of 1000 N through the web, what torque does the spar carry?
(c) Add 10 mm lips turned inward at both flange tips (a lipped channel; build the segment list yourself so the flow starts at a lip's free edge). Recompute the shear centre, the centroid and I_xx. Do the lips move the shear centre toward or away from the web, and why?
(d) Where would you put the hinge line relative to the spar to avoid twisting it, and what else limits that choice?

Answer form: two distances and a centroid; a torque; three properties and a reason; a judgment.

## 6. (Judgment: can the web be thinner?)

A colleague notes that the lesson's spar web has τ = 29 MPa against F_su = 331 MPa and proposes cutting it from 3 mm to 0.5 mm. The web is 160 mm deep between the caps and has vertical stiffeners every 200 mm.

(a) Compute the web stress and the ultimate shear-strength margin for t = 3, 2, 1 and 0.5 mm.
(b) The elastic shear buckling stress of a flat web panel with simply supported edges is τ_cr = k_s π² E/(12(1 − ν²)) (t/b)², with b the shorter side (160 mm) and k_s = 5.35 + 4(b/a)² for a panel of length a ≥ b. Compute τ_cr for each thickness and the ratio τ/τ_cr.
(c) Aircraft webs are often allowed to buckle and carry load as a diagonal tension field. What does the tension field load that an unbuckled web does not? List what else (fastener bearing, minimum gauge, damage tolerance, stiffness) limits the web.
(d) Recommend a web thickness range and the analysis you would require before accepting it.

Answer form: a table of eight numbers; a table of eight numbers; a list; a recommendation.
