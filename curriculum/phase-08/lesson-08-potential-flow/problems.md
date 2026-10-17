# P08.08  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact for the potential-flow model and which the model cannot be trusted for.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

Conventions as in the lesson: W = φ + iψ, dW/dz = u − iv, Γ positive clockwise, lift ρUΓ upward for flow from left to right. ISA density at 3000 m: 0.90913 kg/m³; at sea level 1.225 kg/m³.

## 1. (Warm-up: a lifting cylinder)

A cylinder of radius a = 0.1 m sits in a stream of U = 20 + s m/s of sea-level air, with a clockwise circulation chosen so that Γ/(4πUa) = 0.3.

(a) Find Γ and the angles of the two stagnation points.
(b) Find the largest and smallest surface pressure coefficients and where they occur. Compare Cp at the top (θ = 90°) and the bottom (θ = 270°) and explain the lift from the difference.
(c) Find the lift per unit span from Kutta–Joukowski and from `force_by_pressure`, and the drag. Express the lift as a coefficient based on the diameter and show that it equals 4π × 0.3.

Answer form: a circulation and two angles; four pressure coefficients and a sentence; two lifts, a drag and a coefficient.

## 2. (Derivation: the cylinder, its lift and Blasius's theorem)

(a) Show that W = U(z + a²/z) has ψ = 0 on |z| = a and tends to uniform flow far away, so it is the flow past a cylinder. Show that adding the vortex (iΓ/2π) log z keeps |z| = a a streamline.
(b) Derive the surface speed u_θ = −2U sin θ − Γ/(2πa) and the pressure coefficient. Integrate the surface pressure by hand, F = −∮ p n ds, to show that the drag is zero and the lift is ρUΓ. Find where the stagnation points are, and the Γ at which they merge.
(c) Derive Blasius's theorem F_x − iF_y = (iρ/2)∮(dW/dz)² dz for a body whose surface is a streamline, starting from p = p∞ + ½ρ(U² − |dW/dz|²) and n ds = −i dz. Then show, using the Laurent series of dW/dz far from any body, dW/dz = U + (m + iΓ)/(2πz) + O(z⁻²), that a closed body (m = 0) has F_x = 0 and F_y = ρUΓ: the Kutta–Joukowski theorem for any shape.

Answer form: three derivations; the Γ for merged stagnation points.

## 3. (Unfamiliar scenario: where to put a static port)

The nose of a pitot-static probe is roughly a Rankine half-body: a uniform stream U plus a source m at the origin.

(a) Show that the body is the streamline ψ = m/2 (with the angle measured from the positive x axis, in [0, 2π)), that its half-width far downstream is b = m/(2U), and that the stagnation point is at x = −m/(2πU) = −b/π.
(b) Trace the surface with U = 1, m = 2π (so b = π and the far diameter D = 2π), compute the surface Cp with `pressure_coefficient` and check the zero-normal-flow condition numerically on your traced points.
(c) Find the minimum Cp and where it occurs, in diameters from the nose. How far behind the nose must a static port be for |Cp| < 0.01, and for |Cp| < 0.001? (The far-field behaviour, Cp ≈ −2b/(πx), explains the answer.)
(d) Real pitot-static probes put their static ports only a few diameters behind the nose. What else in a real probe (and in its installation) changes the pressure there, which the half-body leaves out? Why must every probe be calibrated?

Answer form: a derivation; a check; a Cp, a position and two distances in diameters; two or three sentences.

## 4. (Numerical, uses your code: your own aeroplane's circulation)

Your aeroplane has a mass of 550 + 12s kg, a 10 m span and a 1.5 m chord, and cruises at V = 50 + s/2 m/s at 3000 m.

(a) Treat the lift as uniform along the span (a 2-D strip). Find the lift per span, the circulation and c_l.
(b) Find the angle of attack at which a flat plate of the same chord carries that circulation (Kutta condition), and the stagnation points of the equivalent cylinder of radius c/4. What do you notice?
(c) Check the lift with `force_blasius` on the cylinder and with the Blasius integral in the plate plane (`plate_complex_velocity` on a circle of radius 6a). Resolve the plate-plane force into body axes and wind axes. What is the body-axis force along the plate, and what is it?

Answer form: three numbers; an angle, two angles and a sentence; four forces and a sentence.

## 5. (Numerical, uses your code: the Rankine oval family)

For U = 1, a = 1 and m/(πUa) = 0.1, 0.5, 1, 2, 5, 10, 100:

(a) Tabulate the half-length l, the half-thickness t, the fineness ratio l/t and the maximum surface speed.
(b) Show numerically that the force on every oval is zero to round-off.
(c) Explain the two limits of the maximum surface speed: why it tends to 2U as m grows and towards U as m shrinks. Which oval is the better model of a streamlined fuselage, and why does potential flow still not predict its drag?

Answer form: a seven-row table; seven forces; three sentences.

## 6. (Judgment: what potential flow gets right and wrong)

(a) Potential flow predicts zero drag for every closed body (d'Alembert's paradox). Name the physical mechanism that potential flow leaves out and say why the paradox does not make potential flow useless for wings.
(b) A cylinder of radius 50 mm spins in a 10 m/s stream so that its surface speed is twice the stream speed. If the circulation were simply that of the spinning surface, Γ = 2πa × (surface speed), what lift coefficient (on the diameter) would potential flow predict? Experiments on spinning cylinders (the Magnus effect; MIT 2.25 poses it as a problem) find much less lift. Which assumption of the lifting-cylinder model fails?
(c) For a thin aerofoil at 3° in attached flow, which predictions of this lesson's model would you trust (lift, pressure distribution near the leading edge, drag, the stagnation point) and which not?

Answer form: a mechanism and a sentence; a coefficient and an assumption; four verdicts with a reason each.
