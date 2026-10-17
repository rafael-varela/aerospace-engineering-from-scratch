# P04.04  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a reaction wheel)

A reaction wheel is a uniform solid disk of (2 + s/10) kg, radius 0.10 m and thickness 0.02 m.

(a) Its moment of inertia about its spin axis and about a diameter through its centre.
(b) Its angular momentum and kinetic energy at 6000 rpm. Where does the energy come from, and where does it go when the wheel slows down?
(c) The wheel is mounted with its centre 0.5 m from the spacecraft's centre of mass along body axis 1, spin axis along body axis 3. Its contribution to the spacecraft's inertia tensor about the spacecraft's centre of mass. Are there products of inertia? Why?
(d) Starting from rest, its motor applies 0.1 N·m. How long does it take to reach 6000 rpm, and how fast is the 135.4 kg·m² (about axis 3) spacecraft turning the other way by then? (Neglect the wheel's own share of the spacecraft's inertia.)

Answer form: two moments; two numbers and two sentences; a matrix and a sentence; a time and a rate.

## 2. (Derivation: Euler's equations and their consequences)

(a) Derive Euler's equations (5) from dL/dt = τ and the transport theorem (4). Write the three components in principal axes.
(b) Show from (5) with τ = 0 that the kinetic energy and |L|² are constant, without using the inertial frame. Then show that L_i = R I ω is constant in the inertial frame.
(c) Linearise (5) about steady spin ω₂ about axis 2 and derive the growth rate (9). Show that the same calculation about axis 1 or axis 3 gives an oscillation, and find its frequency.
(d) Derive the steady-precession condition (13) for a heavy symmetric top, then the condition (14) for it to have a real solution. Show that the slow root tends to (12) as ω₃ → ∞.
(e) Derive the turning-point condition (15) from conservation of energy and of L_z for a top released at θ₀ with θ̇ = 0 and no precession.

Answer form: five derivations.

## 3. (Unfamiliar scenario: a pencil-shaped satellite)

A satellite is roughly a uniform solid cylinder of (10 + s/5) kg, radius 0.08 m and length 2.0 m. It is released spinning at 750 rpm about its long axis.

(a) Its principal moments of inertia, its angular momentum and its rotational kinetic energy.
(b) Is spin about the long axis stable in the rigid-body model? Use 2.3.
(c) The satellite carries flexible wire antennas that flex and dissipate energy, but no external torque acts. Show that the state of least kinetic energy at fixed |L| is spin about the axis of largest moment, and find the spin rate (rpm) and the kinetic energy there. How much energy did the antennas dissipate?
(d) A mission designer wants passive spin stabilisation. What rule about the spin axis follows from (c)?

Answer form: three numbers; a yes or no with a reason; a rate, an energy and the energy lost; one sentence.

## 4. (Numerical, uses your code: slewing a spacecraft)

Use the 500 kg box satellite of the lesson, I = diag(260.4, 208.3, 135.4) kg·m².

(a) A rest-to-rest bang-bang slew of (20 + s)° about axis 3 in 60 s: the wheel torque, peak body rate and peak wheel momentum.
(b) Simulate it with `simulate_rigid_body` (`max_step=0.5`, so the integrator sees the torque reversal). Report the final angle and body rate, and verify that body plus wheel momentum stays zero.
(c) Now the satellite also carries a constant-speed momentum-bias wheel along body axis 1 storing h₁ = 2 N·m·s (pass `h0=[2, 0, 0]`; it is not driven). Repeat the same slew. Report the largest body rates about axes 1 and 2, the final tilt of body axis 3 from its starting direction, and the final body rate. Repeat with h₁ = 10 N·m·s. Explain the result with (10).
(d) The axis-3 wheel has I_w = 0.02 kg·m², a top speed of 6000 rpm and a maximum torque of 0.1 N·m. Find the fastest rest-to-rest 90° slew about axis 3 allowed by each limit separately, and say which one binds.

Answer form: three numbers; two numbers and a check; four numbers per bias and two sentences; two times and a sentence.

## 5. (Numerical, uses your code: your own top)

Use the lesson's top (100 g disk, r = 30 mm, 4 mm thick, centre 40 mm from the pivot), tilted 30°, spinning at (40 + s) rev/s.

(a) The fast-top and steady-root precession rates, the precession period, and the simulated precession rate (steady start, 1 s).
(b) The minimum spin for steady precession at tilts of 10°, 30°, 60°, 80° and 89°. Why does it fall to zero at 90°?
(c) Release the top with no precession. Simulate 1 s with 20 001 samples. Report the tilt swing, the turning point from (15), the fast-top estimate of the swing, the nutation frequency (from the times of successive tilt maxima, e.g. `scipy.signal.find_peaks`) against L₃/I₁, and the mean precession rate.
(d) Friction slowly removes spin. At what spin does your top stop precessing steadily at 30°, and what does it do then?

Answer form: four numbers; a five-row table and a sentence; five numbers; a spin and a sentence.

## 6. (Judgment: rotation in practice)

For each, give a short argued answer (one to three sentences).

(a) A space telescope uses reaction wheels to point. A constant solar-pressure torque of 10⁻⁵ N·m acts on it about one axis. Why do the wheels eventually need help, and how long until a wheel with 50 N·m·s of capacity saturates?
(b) A figure skater pulls her arms in and spins faster. Is angular momentum conserved? Is kinetic energy conserved? Where does any extra energy come from?
(c) A quadcopter yaws (turns about the vertical) without any rotor tilting. How?
(d) A colleague's simulation of a torque-free body shows |L| drifting by 1 % over a few minutes. Name three mistakes that could cause it.
(e) A fast top is said to "defy gravity". In what sense does gravity still act on it, and what does gravity's work go into when the top is released without precession?
(f) Why does a spin-stabilised spacecraft that must spin about its long axis carry an active nutation damper?

Answer form: six labelled items.
