# P10.01  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical. Draw the axes (chord, stream, balance) and the moment reference point before writing any equation.

Several problems use your personal seed **s**: the day of the month you were born (1–31). Air is ISA sea level (ρ = 1.225 kg/m³, ν = 1.461 × 10⁻⁵ m²/s, a = 340.3 m/s) unless the problem says room air (ρ = 1.204 kg/m³, μ = 1.81 × 10⁻⁵ Pa·s, a = 343.2 m/s). Every dataset here is **SYNTHETIC** unless stated otherwise.

## 1. (Warm-up: the trainer, heavier and faster)

The reference trainer (S = 12.0 m², chord 1.2 m) is loaded to m = 700 + 4s kg and cruises at V = 50 + s/2 m/s at sea level.

(a) Compute the dynamic pressure, the cruise lift coefficient C_L = W/(qS), the chord Reynolds number and the Mach number.
(b) A 1:5 scale model (chord 0.24 m) is tested in a room-air tunnel. At what speed would it match the full-scale Reynolds number, and what Mach number would that be? Is that test possible in a low-speed tunnel?
(c) Which of the two similarity parameters would you match, and which would you correct for, if you could match only one? Answer in two sentences.

Answer form: four numbers; a speed, a Mach number and a sentence; two sentences.

## 2. (Derivation: axes, moments and the aerodynamic centre)

(a) From a sketch of the chord, the stream at angle α and the forces N (normal to the chord, up) and A (along the chord, aft), derive L = N cos α − A sin α and D = N sin α + A cos α, and the inverse. Show that the rotation preserves |F|.
(b) Derive the moment-transfer rule M_P = M_O + (x_P − x_O) N − (z_P − z_O) A for nose-up moments with x aft and z up. Check the sign with a picture of an upward force behind P.
(c) Show that if a section has a point (the aerodynamic centre) about which the moment does not change with incidence, then the moment about any other point x_ref satisfies c_m,ref = c_m,ac + (x_ref − x_ac)/c · c_n, so x_ac = x_ref − c · dc_m,ref/dc_n. Why is c_n exact and c_l only approximately right here?
(d) Differentiate D = N sin α + A cos α with respect to α at fixed N and A and show ∂D/∂α = L. Hence show that a standard uncertainty u_α in the measured angle gives u(c_d) ≈ c_l u_α (α in radians). Evaluate it for c_l = 0.64 and u_α = 0.05°.

Answer form: four derivations and one number.

## 3. (Unfamiliar scenario: a balance mounted at an angle)

A 2-D section model (chord 0.25 m, span 0.50 m) sits on a sting balance whose axes are pitched 2.0° nose-down relative to the model's chord line, so that when the chord is at α = 6.0° to the stream the balance axes are at 4.0°. The balance moment centre lies 50 mm aft and 30 mm below the model's quarter-chord point, measured along the balance axes. At q = 600 Pa the balance reads N_b = 74.9 N, A_b = −4.33 N (that is, forward) and M_b = +2.68 N·m (nose-up). These readings are **SYNTHETIC** numbers chosen for this problem.

(a) Compute L, D, c_l, c_d and c_m about the quarter chord.
(b) Repeat (a) with the mount angle forgotten (rotating by 6.0° instead of 4.0°) and with the 30 mm vertical offset forgotten. Which coefficient does each mistake spoil, and by how much?
(c) Why is the axial force negative even though drag is positive?

Answer form: five numbers; two pairs of numbers and two sentences; one sentence.

## 4. (Numerical, uses your code: your own synthetic balance run)

Generate `synthetic_balance_data(seed=s)` (it is **SYNTHETIC**: the generator invents a NACA 2412-like section and adds instrument noise) and reduce it with `reduce_balance`.

(a) Fit the lift line c_l = a₀(α − α_L0) and the moment line c_m,c/4 = c_m,ac + k c_n. Report a₀, α_L0, x_ac/c and c_m,ac with the 1σ uncertainties from the fit covariance (`np.polyfit(..., cov=True)`).
(b) Repeat for seeds 1 to 20 and compute the standard deviation of x_ac/c across seeds. Compare it with the mean of the stated fit uncertainties. What does the agreement (or disagreement) tell you?
(c) The generator's truth is in `TRUE_SECTION`. Is your seed's x_ac within 2σ of it?

Answer form: four numbers with uncertainties; two numbers and a sentence; yes/no with the number of σ.

## 5. (Numerical, uses your code: where does the drag uncertainty come from?)

For the seed-42 dataset, take the first run at nominal α = 0°, 4° and 8°.

(a) Use `propagate_uncertainty` to find u(c_d) and the variance shares of the five inputs (N, A, M, α, q) at each angle. Check one of them with `monte_carlo_uncertainty`.
(b) Explain why the share of α grows with α and why q hardly matters for c_d.
(c) At α = 8°, what standard uncertainty in α would make the α share equal to the A share? What is u(c_d) then? Could any improvement in the α measurement alone bring u(c_d) below 0.0002?

Answer form: a table; two sentences; an angle, a number and a sentence.

## 6. (Judgment: what transfers from the tunnel to the aeroplane?)

The synthetic tunnel model runs at Re = 7.0 × 10⁵ and M = 0.10; the trainer cruises at Re = 4.5 × 10⁶ and M = 0.16. For each of a₀, α_L0, x_ac, c_m,ac, c_d,min and c_l,max, say whether you expect the tunnel value to transfer to the aeroplane unchanged, with a small correction, or not at all, and give the physical reason in one line each. Then name the one measurement you would add to the tunnel campaign to reduce the largest of these risks.

Answer form: a six-row table and one sentence.
