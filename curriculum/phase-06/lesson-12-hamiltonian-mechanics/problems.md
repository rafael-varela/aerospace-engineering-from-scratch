# P06.12  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a pendulum in phase space)

A pendulum of length l = (0.5 + s/20) m swings with amplitude a = (20 + 2s)° (take a ≤ 170°).

(a) Its ω₀, small-angle period, energy level H = −ω₀² cos a and k² = (E + ω₀²)/(2ω₀²).
(b) Its period by the elliptic integral (8), by the series 2π/ω₀ (1 + a²/16), and by dA/dE with `phase_area`. Which two agree, and to how many digits?
(c) The phase area of the orbit, the largest p on it, and the largest p on the separatrix. Convert the separatrix value to the bob's speed at the bottom (m/s).

Answer form: four numbers; three periods and a sentence; three numbers and a speed.

## 2. (Derivation: the Hamiltonian toolkit)

(a) Starting from the definition (2) and Lagrange's equations, derive Hamilton's equations (3). Where did you use pᵢ = ∂L/∂q̇ᵢ?
(b) Prove T = dA/dE for a one-degree-of-freedom Hamiltonian H = p²/2 + V(q) and a closed orbit.
(c) Prove Liouville's theorem in one degree of freedom from the divergence of the phase velocity. Then show that with ṗ = −V′(q) − cp the area of every region decays as e^{−ct}.
(d) Compute the Jacobian determinant of a kick p ← p + h a(q), of a drift q ← q + h p, and of an explicit Euler step. Conclude which of Euler, symplectic Euler and Verlet preserve area.

Answer form: four derivations.

## 3. (Unfamiliar scenario: a bead on a spinning hoop)

A small mass m slides without friction on a circular hoop of radius R = 0.5 m that is spun about its vertical diameter at a constant rate Ω = (4 + s/5) rad/s, as a deployable boom's hinge mass might on a spinning spacecraft test rig. θ is measured from the bottom of the hoop.

(a) Write L = ½mR²(θ̇² + Ω² sin²θ) + mgR cos θ and find H with `legendre_transform`. Is H equal to T + V? Is H conserved? Is T + V conserved? Explain.
(b) Find the equilibria from Hamilton's equations. Show that the bottom is stable only for Ω < Ω_c = √(g/R), and compute Ω_c.
(c) For your Ω, find every equilibrium, classify it, and give the small-oscillation frequency about each stable one (or the growth rate about an unstable one).
(d) Simulate 10 s from θ = θ₀ + 0.2 rad, θ̇ = 0, about your off-axis equilibrium θ₀ (if it exists; otherwise about the bottom). Report the spread (max − min) of H and of T + V per unit mass along the run, and the range of θ.

Answer form: a Hamiltonian and three yes/no answers with reasons; an inequality and a rate; equilibria with frequencies; two spreads and a range.

## 4. (Numerical, uses your code: your own blob)

With ω₀ = π rad/s, take a blob of 4000 boundary points centred on (q, p) = (0.3, (0.5 + s/50)ω₀) with semi-axes 0.2 rad and 0.2ω₀.

(a) Its energy level at the centre (as a fraction of ω₀²) and its area. Is it entirely inside the separatrix?
(b) Carry it 5 s with the exact flow (`flow_blob`) and report the relative area change. What sets it?
(c) Carry it with RK4 and Verlet for 5 s and with Euler for 1 s, at h = 0.02, 0.04 and 0.08 s. Tabulate each relative area change against the exact flow (Euler against the initial area).
(d) From (c), estimate the power of h in RK4's area loss over a fixed time, and explain the Verlet column.

Answer form: a fraction, an area and a yes/no; a number and a cause; a 3 × 3 table; an exponent and a sentence.

## 5. (Numerical, uses your code: an orbit-determination error)

A satellite is at periapsis of an orbit with a = (6800 + 20s) km and e = 0.001s. Orbit determination has overestimated its periapsis speed by 5 mm/s.

(a) Find the error δa it causes, to first order (from E = v²/2 − μ/r) and exactly (from the perturbed state's energy).
(b) Predict the along-track error after 1 and after 10 orbits with (13).
(c) Check (b) with two `solve_ivp` runs (DOP853, rtol 10⁻¹²): project the separation on the velocity direction.
(d) Compute the radial action J_r of the orbit and check that dJ_r/dE gives its period.

Answer form: two lengths; two lengths; two lengths and a comparison; an action and a period.

## 6. (Judgment: phase space in practice)

Give a short argued answer (one to three sentences) to each.

(a) "Liouville's theorem says the uncertainty of an orbit never grows." What is right and what is wrong in that sentence?
(b) A colleague propagates a cloud of 10 000 orbit states for 1000 orbits with an adaptive RK45 at rtol 10⁻⁶ and finds that its sample covariance determinant has fallen by 30 %. Is that physics or numerics, and what would you change?
(c) For a 90-minute propagation of a single spacecraft state to plan a manoeuvre, does it matter whether the integrator is symplectic? Why or why not?
(d) For a satellite in a 300 km orbit, which assumption of this lesson fails first, and in which direction does it push the along-track error?

Answer form: four labelled items.
