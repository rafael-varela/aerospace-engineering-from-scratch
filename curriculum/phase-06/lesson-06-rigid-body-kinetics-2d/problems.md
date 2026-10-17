# P06.06  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a wheel on a ramp)

A ground-handling dolly wheel of 8 kg, radius 0.25 m and radius of gyration 0.7 R is released on a ramp of (4 + s/2)°. The static friction coefficient is 0.15 and the kinetic 0.12.

(a) Does it roll? Give its acceleration, its angular acceleration and the friction force.
(b) The smallest static friction coefficient that would let it roll on this ramp.
(c) The steepest ramp on which it rolls with μ_s = 0.15. On a 30° ramp, give the acceleration of the centre, the angular acceleration and the rate at which the slip velocity grows.

Answer form: three numbers and a yes or no; a number; an angle and three numbers.

## 2. (Derivation: the spin-up in closed form)

(a) From (5) derive the slip-velocity equation (6) and the effective mass M_eff.
(b) Use impulse–momentum to derive (7) for the friction impulse and the final rolling state. Explain why μ_k does not appear.
(c) Derive the friction energy (9) from ∫F v_s dt and (6).
(d) Starting from (10) with V_f and ω_f from (7), show algebraically that the kinetic-energy loss equals ½ M_eff v_s0². (Hint: write V_f = V₀ − J/M and ω_f R = ω₀R + JR²/I.)
(e) Derive both branches of the spin-up time (8) for a linearly ramped load, and the slip distance for the first branch.
(f) Show that for M → ∞ the friction energy equals the wheel's final kinetic energy, and that pre-spin to a fraction f of V₀/R reduces the energy by (1 − f)².

Answer form: six derivations.

## 3. (Unfamiliar scenario: a canister thrown with backspin)

A thin-walled gas canister (treat it as a cylinder with I = 0.8 m R²) of 3 kg and radius 0.15 m is rolled along a hangar floor at (2 + s/10) m/s with a backspin of (30 + s) rad/s (spinning so that its top moves backwards). The kinetic friction coefficient with the floor is 0.3.

(a) Write (5) for this body. What plays the role of M, and of N?
(b) Will the canister come back towards you? Derive the condition for it to return in terms of V₀, ω₀, R and k² = I/(mR²).
(c) Its final velocity, the time it slides and the distance its centre moves before it rolls.
(d) The energy dissipated, by both routes. What fraction of the initial kinetic energy is it?
(e) Check your answers with `spin_up_closed_form` and `simulate_touchdown` (build a `Wheel(m_wheel=3, R=0.15, k_ratio=√0.8, M=3)`).

Answer form: two sentences; a yes or no and an inequality; three numbers; two numbers and a fraction; a check.

## 4. (Numerical, uses your code: your own touchdown)

Use the lesson's wheel with V₀ = (65 + s/5) m/s, μ_k = 0.35 + s/100, μ_s = 0.8, and a normal load ramping to (100 + s) kN in 0.08 s.

(a) By hand (with a calculator), the friction impulse, which branch of (8) applies, the spin-up time, the final speed and wheel rate, the slip distance, the friction energy and the peak drag load.
(b) Simulate 0.3 s with `simulate_touchdown`. Report the differences from (a) for the spin-up time, the friction energy and the slip distance, and the largest energy-audit residual.
(c) Repeat (b) with rtol = atol = 10⁻⁶. Which quantity degrades most?
(d) Plot the slip velocity against time and mark the end of the ramp. Explain the shape of each piece.

Answer form: seven numbers; four numbers; four numbers and a sentence; a plot and two sentences.

## 5. (Numerical, uses your code: braking at the skid limit)

After the lesson's touchdown (70 m/s, μ_k = 0.55, 0.1 s ramp to 120 kN), the tyre's static friction coefficient is μ_s = 0.7 + s/100.

(a) The largest brake torque without skid from (11), and the simpler estimate μ_s N R. Why is (11) larger?
(b) Apply a constant brake torque of 0.95 and of 1.05 times that limit from t = 0.3 s and simulate to 0.6 s. Report the phase sequence, the final deceleration of the aircraft share and, for the second run, the time from the start of the skid to lock-up (check it against ω I/(T_b − μ_k N R)).
(c) For both runs, the energy absorbed by the brake and the friction heat up to 0.6 s, and the energy-audit residual. Where does most of the energy go in each case?
(d) If each deceleration were held from 70 m/s to rest, how far would the aircraft share travel? What does this say about anti-skid design?

Answer form: two numbers and a sentence; two phase lists, two decelerations, a time; four energies and a sentence; two distances and a sentence.

## 6. (Judgment: touchdown in practice)

For each, give a short argued answer (one to three sentences).

(a) The friction energy does not depend on μ_k. Does that mean the tyre wears the same on a wet and a dry runway?
(b) Engineers have repeatedly proposed spinning the wheels up before touchdown. From (9), what is the benefit at 90 % pre-spin, and what goes wrong if the wheel is spun *faster* than V₀/R?
(c) A colleague takes moments about the contact point and writes I α = −T_b, since friction and the normal force pass through it. What is wrong?
(d) A colleague's simulation shows the wheel's rim speed overshooting the axle speed and oscillating about it after spin-up. What caused it, and how is it fixed?
(e) Why does the spin-up drag load matter for the landing-gear leg even though it lasts only a tenth of a second?
(f) Which input of the lesson's model would you most want measured before trusting its spin-up time, and which output would you trust without it?

Answer form: six labelled items.
