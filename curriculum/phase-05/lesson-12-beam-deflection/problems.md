# P05.12  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical. Draw the free-body diagram and the supports before writing any equation.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

Material data (typical, not design allowables): 2024-T3 E = 73.1 GPa, ρ = 2780 kg/m³; 7075-T6 E = 71.7 GPa, G = 26.9 GPa.

## 1. (Warm-up: a pitot boom)

A pitot boom is a 2024-T3 tube, 25 mm outside and 21 mm inside diameter, cantilevered L = 0.6 + s/100 m forward of a wing tip, with a 0.3 kg probe head at its end. Size it for a downward load factor of 6 (a hard landing).

(a) Compute I, EI, the tip load from the probe and the distributed load from the tube's own mass.
(b) Using `cantilever_point` and `cantilever_udl`, find the tip deflection from each load and their sum, and the tip slope. Which result uses superposition?
(c) Compute the root bending moment and the bending stress (`P05.10`). Is the boom limited by stiffness or by strength?

Answer form: four numbers; three deflections and a slope; a moment, a stress and a sentence.

## 2. (Derivation: integration, superposition and the elliptic wing)

(a) Derive the cantilever deflection under a uniform load, w = qx²(6L² − 4Lx + x²)/(24EI), by integrating EI w″ = M twice and applying the boundary conditions. Then derive the same result from EI w⁗ = q with four boundary conditions, and say which route needs the reactions.
(b) Derive the prop force R = 3qL/8 of a propped cantilever by superposition. Then derive it again by integrating EI w″ = M(x; R) with R unknown and using the three kinematic conditions. Which condition did superposition use?
(c) Derive the tip slope and deflection of a uniform cantilever under elliptic lift, equation (5) of the lesson, by superposing point loads. You will need ∫₀¹ t²√(1 − t²) dt = π/16 and ∫₀¹ t³√(1 − t²) dt = 2/15; prove the second with the substitution u = 1 − t².

Answer form: three derivations.

## 3. (Unfamiliar scenario: a strut-braced wing)

A high-wing aeroplane's 5 m semi-span is braced by a strut attached to the spar at a = 2.75 + s/100 m from the root. Treat the wing as a beam with a moment-resisting root fitting (clamped) and the strut as a rigid vertical support at a. The lift is uniform, 13 974.5 N in total, and the spar has EI = 3.479 × 10⁵ N·m².

(a) Why is this beam statically indeterminate? How many redundant reactions does it have?
(b) Find the strut force by superposition of `cantilever_udl` and `cantilever_point`, using the condition at the strut.
(c) Find the root bending moment, the moment at the strut, and the tip deflection, and compare them with the unbraced cantilever.
(d) Check your answer with `fd_beam`: apply the strut force as a point load (on a node) and show the deflection at the strut is zero.
(e) If the root fitting were a pin instead, the beam would be determinate. Find the strut force from statics alone and compare.

Answer form: a sentence and a count; a force; three numbers and two comparisons; a check; a force.

## 4. (Numerical, uses your code: the tip deflection of your own wing)

Take your aeroplane and spar from `P05.10` Problem 4 (mass 550 + 12s kg, n = 3.8, semi-span 4.5 + s/20 m, depth h = 0.14 + s/500 m, 10 mm caps of the width you found (46 mm for s = 15), 3 mm web), elliptic lift, 7075-T6.

(a) Compute EI at the root and the tip deflection and slope by `elliptic_tip`, by `deflection_by_integration` and by `fd_beam` (n = 400). Express the deflection as a percentage of the semi-span.
(b) Repeat with the stiffness tapering linearly to 40 % of the root value at the tip (`tapered_EI(EI_root, 0.6, s)`).
(c) A stiffness requirement limits the tip deflection at limit load to 5 % of the semi-span. By what factor must the root EI grow (with the same taper) to meet it? Why can you answer from (b) without another solve?

Answer form: three deflections, a slope and a percentage; two deflections and a percentage; a factor and a sentence.

## 5. (Numerical, uses your code: testing the finite-difference solver)

(a) Show that `fd_beam` converges at order 2 for a clamped–clamped beam (L = 1 m, EI = 10⁴ N·m², q = 100 N/m), whose exact deflection is qx²(L − x)²/(24EI). Derive that closed form first.
(b) Put a 500 N downward point load on a simply supported beam (same L and EI) at a = 0.3 m (on a node for n = 20–320) and at a = 0.3333 m (between nodes; `fd_beam` snaps it to the nearest node). Compare the maximum error against `simply_supported_point` for both, and their observed orders. Explain the difference.
(c) Verify Maxwell's reciprocal theorem for the propped cantilever with `fd_beam`: compare w(0.7 m) under a unit load at 0.3 m with w(0.3 m) under a unit load at 0.7 m (n = 100). Is the difference what you would expect?

Answer form: a derivation and three orders; two error tables and an explanation; two numbers and a sentence.

## 6. (Judgment: why is the measured tip deflection different?)

In a static test of the lesson's wing loaded to limit load in the elliptic distribution, the measured tip deflection is 15 % below the spar-alone prediction of 243 mm.

(a) Estimate the tip deflection from shear deformation of the web, ≈ M₀/(G A_web) with A_web = 160 × 3 mm². Does it explain the difference, and in which direction?
(b) Estimate the error of the small-slope approximation at the measured slope.
(c) What factor on EI would explain the measurement? Which part of the real wing could supply it?
(d) The root fitting adds deflection, not less. If the fitting stiffness were k_θ = 10⁷ N·m/rad, what EI factor would then be needed to match the measurement?
(e) List two further measurements that would separate these explanations.

Answer form: a deflection, a percentage and a direction; an error; a factor and a source; a factor; two measurements.
