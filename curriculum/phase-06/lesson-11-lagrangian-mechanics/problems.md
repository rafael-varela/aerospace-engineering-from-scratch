# P06.11  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical. Draw the coordinates and the free-body diagrams before writing any equation.

Several problems use your personal seed **s**: the day of the month you were born (1–31). All slung-load numbers are illustrative.

## 1. (Warm-up: your own slung load)

A hook block of m₁ = (60 + 2s) kg hangs on a cable of l₁ = (10 + s/5) m; a load of m₂ = (1000 + 20s) kg hangs l₂ = 3 m below the hook. The helicopter hovers (fixed suspension point).

(a) Write M₀ and K₀ of (9) with numbers.
(b) Solve the frequency equation (10) by hand (a quadratic in ω²) for the two natural frequencies and periods, and find the two mode shapes θ₂/θ₁.
(c) Compare the slow period with that of a single pendulum of length l₁ + l₂, and the fast frequency with the hook-mode estimate √(((m₁ + m₂)g/l₁ + m₂g/l₂)/m₁).
(d) Check your numbers with `normal_modes`.

Answer form: two matrices; two frequencies, two periods and two ratios; two comparisons in per cent; a check.

## 2. (Derivation: one system, two routes)

(a) Starting from the positions (1), derive T and V, equation (3), and then the equations of motion (5) by Lagrange's equations. Show that the θ̇₁θ̇₂ terms cancel.
(b) Prove that the energy function h = Σ q̇ᵢ ∂L/∂q̇ᵢ − L is conserved when ∂L/∂t = 0, and that h = T + V for this system.
(c) From the Newton–Euler equations (7), eliminate the tensions: dot the equation of mass 2 with n₂, and the sum of the two mass equations with n₁. Show that you recover (5).
(d) Derive the frequency equation (10) from det(K₀ − ω²M₀) = 0 and find its limits when m₂ → 0 and when m₁ → 0.

Answer form: four derivations.

## 3. (Unfamiliar scenario: a stretchy cable)

Real cables stretch. Model a load m = 1500 kg on a massless elastic cable of unstretched length l₀ = 15 m and stiffness k, hanging from a fixed point, swinging in a plane. Use the generalized coordinates r (the cable's current length) and θ (its angle from the vertical).

(a) Write T and V (gravity and the spring) and derive the two equations of motion with `euler_lagrange`.
(b) Find the static stretch and the equilibrium length r₀. Linearise about it and show that the bounce and swing frequencies are √(k/m) and √(g/r₀).
(c) Choose k so that the bounce frequency is exactly twice the swing frequency (k = 4mg/r₀ with r₀ = l₀ + mg/k: solve for k). Release the load with r = r₀ + (1 + s/50) m, ṙ = 0, and a tiny swing θ = 10⁻³ rad, θ̇ = 0. Simulate 300 s and plot the bounce energy ½mṙ² + ½k(r − r₀)² and the swing energy ½mr²θ̇² + mgr(1 − cos θ) against time. (This cable is far softer than steel; the stiffness is chosen to show the effect.)
(d) What do you observe, and why is a 2 : 1 ratio between a bounce and a swing frequency something a slung-load engineer avoids?

Answer form: two equations; a stretch, a length and two frequencies; a stiffness and a plot; two sentences.

## 4. (Numerical, uses your code: checking and reading a simulation)

Use your slung load from Problem 1.

(a) With `numpy.random.default_rng(s)`, draw 1000 random states (angles in (−π, π), rates in (−3, 3) rad/s) and report the largest relative difference between `lagrange_accelerations` and `newton_euler`.
(b) Release the load from rest with θ₁ = 10°, θ₂ = 0 (not a mode shape). Simulate 120 s. Using the mode shapes, predict how much of each mode this release contains (write q(0) = c₁φ₁ + c₂φ₂ and solve for c₁, c₂). Then find both frequencies in the FFT of θ₁(t) and compare them with Problem 1.
(c) For the same run, report the largest and smallest cable tension T₁ and sling tension T₂ as multiples of their static values. Why does T₂ vary much more than in a mode-1 swing of the same θ₁?
(d) Report the largest energy error at rtol = 10⁻⁶ and 10⁻¹⁰.

Answer form: a number; two coefficients, two frequencies and a comparison; four ratios and a sentence; two numbers.

## 5. (Derivation and scenario: a generalized force — forward flight)

In steady slow forward flight the load has a horizontal drag D (opposite to the flight direction, take it in −x) and the hook block has none. The suspension point moves at constant velocity, so its frame is inertial.

(a) Compute the generalized forces Q₁ = D_vec · ∂r₂/∂θ₁ and Q₂ = D_vec · ∂r₂/∂θ₂ by virtual work.
(b) Find the steady trail angles θ₁, θ₂ from Lagrange's equations with Q (set all rates and accelerations to zero). Show that tan θ₂ = −D/(m₂g) and tan θ₁ = −D/((m₁ + m₂)g).
(c) Check (b) by statics with free-body diagrams of the two masses.
(d) For the lesson's slung load and D = (2 + s/10) kN, find the two angles in degrees and the cable tension. Check with `newton_euler` that the state is an equilibrium (add the drag to the right-hand side of (7) yourself, or argue from (b)).

Answer form: two expressions; two angles in closed form; a check; two angles, a tension and a check.

## 6. (Judgment: what can a simulation of a chaotic system be trusted for?)

(a) A colleague simulates the lab pendulum released from 120° and reports the angle of the lower link at t = 30 s to four decimals. Using the lesson's numbers, estimate what initial-state precision would be needed for that number to be meaningful, and comment.
(b) Which quantities of the same simulation can be trusted at t = 30 s? Give two, with the reason.
(c) Is the slung load in hover chaotic at a 10° swing? Give a test you could run with your code that would show it, and the result you expect.
(d) A flight-control engineer proposes using the lesson's model to design a notch filter for the hook mode. Name the two assumptions in the model that most affect the hook-mode frequency in reality, and what you would measure to pin them down.

Answer form: a precision and a sentence; two quantities; a test and a prediction; two assumptions and two measurements.
