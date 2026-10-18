# P10.03  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. Circulation Γ is positive clockwise (Phase 10 convention), so L′ = ρV∞Γ is positive for lift. Say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a cylinder with circulation)

A cylinder of radius a = 0.5 m sits in a stream V = 10 m/s of sea-level air (ρ = 1.225 kg/m³) with clockwise circulation Γ = g · 4πaV, g = 0.1 + s/100.

(a) Find Γ and the angles of the two stagnation points on the surface.
(b) Find the minimum and maximum surface Cp and where they occur.
(c) Find the lift per metre of span from Kutta–Joukowski and the lift coefficient based on the diameter. Check the lift by integrating the surface pressure with `pressure_force_coefficients`.
(d) Sketch the streamlines for g = 0, g, 1 and 1.5. What happens at g = 1?

Answer form: Γ and two angles; two Cp values with angles; a force, a coefficient and a check; a sketch and a sentence.

## 2. (Derivation: the Kutta condition fixes the lift)

(a) From the complex potential of §2.2, show that the tangential velocity on the circle is u_θ = −2V sin(θ − α) − Γ/(2πa).
(b) The Joukowski or Karman–Trefftz map has dz/dζ = 0 at ζ = b, the image of the trailing edge, which lies on the circle at θ = −β. Explain why the airfoil-plane velocity there is infinite unless the circle-plane velocity vanishes, and derive Γ = 4πVa sin(α + β).
(c) Show that the lift coefficient is c_l = 8πa sin(α + β)/c, and that for a thin, symmetric Joukowski section (a → b, c → 4b, β → 0) it tends to 2π sin α. What do β and the ratio 4a/c separately represent?
(d) Integrate the cylinder's surface pressure analytically (∮ sin θ and ∮ sin² θ dθ) to prove L′ = ρVΓ and D′ = 0 for the cylinder.

Answer form: four derivations.

## 3. (Unfamiliar scenario: a rotor sail)

Some ships carry Flettner rotors: tall vertical cylinders spun by a motor, which develop a side force in the wind (the Magnus effect). Take a rotor 4.0 m in diameter and 25 m tall, a wind of 10 m/s across it, and a surface speed of 30 m/s (spin ratio 3).

(a) Estimate the circulation as the surface speed times the circumference (the flow is dragged round at the surface speed). What is g = Γ/(4πaV)? Where are the stagnation points?
(b) Estimate the side force per metre and on the whole rotor from Kutta–Joukowski, and the lift coefficient based on the diameter.
(c) Real rotors develop much less force than (b) and also have drag. Name the two physical effects that potential flow leaves out here, and say which assumption of (a) is the weakest.

Answer form: Γ, g and a sentence; a force per metre, a force and a coefficient; two effects and a sentence.

## 4. (Numerical, uses your code: a section like yours)

Fit a Karman–Trefftz section with thickness t = 0.08 + s/250 (fraction of chord), camber 0.02 and trailing-edge angle 16° (`fit_karman_trefftz`).

(a) Report the circle centre μ, radius a, angle β, the chord and the chord-line inclination.
(b) Find the zero-lift angle (relative to the chord line), the lift-curve slope as a multiple of 2π, and the angle at which the section gives the trainer's cruise c_l = 0.3309.
(c) At that angle, compute c_l four ways (Kutta–Joukowski; surface pressure with n = 1024; Blasius and momentum on an ellipse of semi-axis 0.8c; the circulation line integral on the same ellipse) and the pressure drag. How close do they agree?

Answer form: five numbers; three numbers; five numbers and a sentence.

## 5. (Numerical, uses your code: how fast does pressure integration converge?)

For a section with t = 0.12 and m = 0.02 at α = 4°, compute the error of the pressure-integrated c_l against the exact c_l for n = 16, 32, 64, 128, 256 surface points, for τ = 0 (Joukowski, cusped trailing edge) and τ = 16°.

(a) Tabulate the errors and the observed orders.
(b) Explain why one case converges "spectrally" (faster than any power of 1/n) and the other only algebraically. Hint: look at dz/dζ near ζ = b.

Answer form: a table and orders; a short argument.

## 6. (Judgment: is the exact solution the right answer?)

The exact potential-flow lift-curve slope of a 12 % symmetric Karman–Trefftz section (τ = 16°) is 1.094 × 2π. The authentic NACA 0012 data of Ladson (NASA TM-4074, `curriculum/phase-10/data/ladson_naca0012_M015.csv`, Table I, free transition) give a slope that you can fit between α ≈ −4° and +6°.

(a) Fit the measured slope at Re = 2.0 × 10⁶ and 5.97 × 10⁶ and express it as a multiple of 2π.
(b) The potential flow is exact for its assumptions. Which assumption explains the difference, and why does the difference shrink as Re grows?
(c) The Kutta condition is sometimes called "a viscous effect imposed on an inviscid flow". In two or three sentences, explain what that means and why it still works.

Answer form: two numbers; two sentences; a short paragraph.
