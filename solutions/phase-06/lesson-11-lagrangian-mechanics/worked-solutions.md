# P06.11  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for **s = 15**; numerical results are from `solution.py` and the probe scripts that call it, given to the figures shown. Every "simulated" value is an output of the model, not a measurement of a helicopter. All slung-load numbers are illustrative.

## 1. Your own slung load (s = 15: m₁ = 90 kg, l₁ = 13 m, m₂ = 1300 kg, l₂ = 3 m)

**Assumptions.** Point masses, rigid massless links, fixed hover point, small angles.

(a) M₀ = [[(90 + 1300) × 13², 1300 × 13 × 3], [1300 × 13 × 3, 1300 × 3²]] = **[[234 910, 50 700], [50 700, 11 700]] kg·m²**; K₀ = diag(1390 × 9.80665 × 13, 1300 × 9.80665 × 3) = **diag(177 206, 38 246) N·m/rad**.

(b) a = m₁l₁l₂ = 3510, b = −(m₁ + m₂)g(l₁ + l₂) = −218 100, c = (m₁ + m₂)g² = 133 677 (SI). ω² = 0.61908 and 61.5176 s⁻², so **ω₁ = 0.78682 rad/s (period 7.986 s)** and **ω₂ = 7.84332 rad/s (period 0.801 s)**. Shapes from the first row, θ₂/θ₁ = (K₁₁ − ω²M₁₁)/(ω²M₁₂): **+1.0124** (cable and sling together) and **−4.5765** (hook against load).

(c) A single 16 m pendulum has a period of 8.026 s; mode 1 is **0.50 % shorter**. The hook-mode estimate gives 7.6725 rad/s, **2.2 % below** ω₂ (it treats the load as fixed, which makes the system softer than it is).

(d) `normal_modes(DoublePendulum(90, 1300, 13, 3))` returns the same frequencies and shapes to all printed digits.

**Check.** m₂/m₁ = 14.4 is large, so mode 1 should be close to the single pendulum of length l₁ + l₂ (2.6): it is, to 0.5 %.

**Common mistakes.** Forgetting that M₁₁ carries (m₁ + m₂), not m₁; taking √ of the quadratic's roots in the wrong units (the roots are ω², in s⁻²).

## 2. One system, two routes

(a) ẋ₁ = l₁θ̇₁ cos θ₁, ẏ₁ = l₁θ̇₁ sin θ₁, ẋ₂ = ẋ₁ + l₂θ̇₂ cos θ₂, ẏ₂ = ẏ₁ + l₂θ̇₂ sin θ₂. Then v₁² = l₁²θ̇₁² and v₂² = l₁²θ̇₁² + l₂²θ̇₂² + 2l₁l₂θ̇₁θ̇₂ cos Δ, which gives T in (3); V = m₁gy₁ + m₂gy₂ gives V in (3). For θ₁:
∂L/∂θ̇₁ = (m₁ + m₂)l₁²θ̇₁ + m₂l₁l₂θ̇₂ cos Δ,
d/dt(∂L/∂θ̇₁) = (m₁ + m₂)l₁²θ̈₁ + m₂l₁l₂θ̈₂ cos Δ − m₂l₁l₂θ̇₂ sin Δ (θ̇₁ − θ̇₂),
∂L/∂θ₁ = −m₂l₁l₂θ̇₁θ̇₂ sin Δ − (m₁ + m₂)gl₁ sin θ₁.
Subtracting, the −m₂l₁l₂θ̇₁θ̇₂ sin Δ terms cancel and **(m₁ + m₂)l₁²θ̈₁ + m₂l₁l₂ cos Δ θ̈₂ = −m₂l₁l₂θ̇₂² sin Δ − (m₁ + m₂)gl₁ sin θ₁**. The θ₂ equation is the same calculation with the sign of sin Δ reversed in the derivative of cos Δ: **m₂l₁l₂ cos Δ θ̈₁ + m₂l₂²θ̈₂ = m₂l₁l₂θ̇₁² sin Δ − m₂gl₂ sin θ₂**. That is (5).

(b) Multiply the i-th Lagrange equation (with Q = 0) by q̇ᵢ and sum: Σq̇ᵢ d/dt(∂L/∂q̇ᵢ) − Σq̇ᵢ∂L/∂qᵢ = 0. Add and subtract Σq̈ᵢ∂L/∂q̇ᵢ: d/dt(Σq̇ᵢ∂L/∂q̇ᵢ) − [Σq̇ᵢ∂L/∂qᵢ + Σq̈ᵢ∂L/∂q̇ᵢ] = 0. When ∂L/∂t = 0 the bracket is dL/dt, so **dh/dt = 0 with h = Σq̇ᵢ∂L/∂q̇ᵢ − L**. V does not depend on q̇ and T = ½q̇ᵀMq̇ is homogeneous of degree 2 in q̇, so Σq̇ᵢ∂T/∂q̇ᵢ = 2T (Euler's theorem) and **h = 2T − (T − V) = T + V**.

(c) With nᵢ · eᵢ = 0, nᵢ · nⱼ = cos(θᵢ − θⱼ), e₁ · n₂ = sin Δ and e₂ · n₁ = −sin Δ:
Mass 2 dotted with n₂ (T₂ drops out): m₂[l₁θ̈₁ cos Δ − l₁θ̇₁² sin Δ + l₂θ̈₂] = −m₂g sin θ₂. Multiplied by l₂, this is the second row of (5).
Sum of the two mass equations, m₁a₁ + m₂a₂ = −T₁e₁ + (m₁ + m₂)g (T₂ cancels), dotted with n₁ (T₁ drops out): (m₁ + m₂)l₁θ̈₁ + m₂l₂θ̈₂ cos Δ + m₂l₂θ̇₂² sin Δ = −(m₁ + m₂)g sin θ₁. Multiplied by l₁, this is the first row of (5). The other two components (along e₂ and e₁) give the tensions.

(d) det(K₀ − ω²M₀) = ((m₁ + m₂)gl₁ − ω²(m₁ + m₂)l₁²)(m₂gl₂ − ω²m₂l₂²) − ω⁴m₂²l₁²l₂². Dividing by m₂l₁l₂: (m₁ + m₂)(g − ω²l₁)(g − ω²l₂) − ω⁴m₂l₁l₂ = 0, which expands to **m₁l₁l₂ω⁴ − (m₁ + m₂)g(l₁ + l₂)ω² + (m₁ + m₂)g² = 0**. Limits: m₂ → 0 gives m₁(l₁ω² − g)(l₂ω² − g) = 0, two independent pendulums; m₁ → 0 kills the ω⁴ term and leaves ω² = g/(l₁ + l₂), a single pendulum of the total length.

**Check.** Units: each term of (10) is kg·m²·s⁻⁴. The `test_sympy_mass_matrix_and_forcing_match_the_hand_forms` and `test_lagrange_equals_newton_euler_at_random_states` tests check (a) and (c) numerically.

**Common mistakes.** Dropping the (θ̇₁ − θ̇₂) factor from the derivative of cos Δ, which leaves a spurious θ̇₁θ̇₂ term; getting the sign of e₂ · n₁ wrong.

## 3. A stretchy cable

**Assumptions.** Massless elastic cable, linear spring, planar motion, fixed suspension.

(a) x = r sin θ, y = −r cos θ, T = ½m(ṙ² + r²θ̇²), V = −mgr cos θ + ½k(r − l₀)². `euler_lagrange` gives **m r̈ − m r θ̇² − m g cos θ + k(r − l₀) = 0** and **m r (r θ̈ + 2 ṙ θ̇ + g sin θ) = 0**: the radial equation has the centripetal term r θ̇², the angular one the Coriolis term 2ṙθ̇ (`P06.01`).

(b) At rest with θ = 0: k(r₀ − l₀) = mg, a static stretch **mg/k** and **r₀ = l₀ + mg/k**. Linearising: m δr̈ + k δr = 0 (bounce, **√(k/m)**) and r₀θ̈ + gθ = 0 (swing, **√(g/r₀)**); to first order they are uncoupled.

(c) √(k/m) = 2√(g/r₀) ⇒ k r₀ = 4mg ⇒ k(l₀ + mg/k) = 4mg ⇒ **k = 3mg/l₀ = 2942.0 N/m**, r₀ = **20.0 m** (a 5.0 m static stretch), bounce 1.4005 rad/s, swing 0.7002 rad/s. Released 1.3 m below r₀ with θ = 10⁻³ rad: the bounce starts with 2486 J and the swing with almost none. The swing energy grows slowly (2.3 J at 50 s, 74 J at 100 s), then fast; at **171.7 s the swing holds 2486 J and the bounce 0.1 J**, the load swinging ±7.43°; the energy then flows back (bounce 2485 J, swing 1.6 J at 300 s). The total energy is held to 1.5 × 10⁻⁵ J.

(d) The bounce and the swing exchange their energy completely and periodically: an **autoparametric (internal) resonance**. The bounce modulates the pendulum length at twice the swing frequency, which is exactly the frequency at which a length modulation pumps a pendulum. A slung load whose vertical bounce sits near twice its pendulum frequency can turn a vertical jolt (a gust, a collective input) into a large swing a minute later, with no obvious cause. The linear model of (b) cannot show it: the coupling is in the nonlinear terms.

**Check.** At the moment of full transfer the swing energy equals the initial bounce energy to four figures, as energy conservation requires.

**Common mistakes.** Measuring the potential from the unstretched length but the swing about r₀ inconsistently; using l₀ instead of r₀ in the swing frequency.

## 4. Checking and reading a simulation (s = 15)

(a) 1000 random states from `default_rng(15)`: largest relative difference **2.3 × 10⁻¹⁴**: round-off.

(b) Solving Φc = (10°, 0) with Φ = [[1, 1], [1.0124, −4.5765]]: **c₁ = 8.189°, c₂ = 1.811°** (θ₁ components). The FFT of a 120 s run has its largest peak at **0.785 rad/s** (one bin is 0.052 rad/s; the linear value 0.7868 is within it) and a second peak, about a fifth as high, at **7.12 rad/s**, **9 % below** the linear 7.843 rad/s. That is not a bug: a pure mode-2 release has frequency 7.840 rad/s at θ₁ = 0.1°, 7.586 at 1° and 7.091 at 1.81°. The hook mode softens with amplitude far faster than mode 1, because its relative link angle Δ is large (about 10° here) and the coupling cos Δ, sin Δ terms of (5) are nonlinear in Δ.

(c) Over the run, T₁ ranges **0.684–1.472** times (m₁ + m₂)g and T₂ **0.673–1.503** times m₂g. In a pure mode-1 swing with the same θ₁ amplitude, T₂ ranges only 0.985–1.031. The hook mode drives the large swings in tension: the hook accelerates sideways at about l₁c₂ω₂² ≈ 13 × 0.0316 × 61.5 ≈ 25 m/s², and the links' tensions supply that acceleration. A release that is not a mode shape loads the hook far more than its angle suggests.

(d) Largest energy error: **0.271 J** at rtol 10⁻⁶ and **3.9 × 10⁻⁶ J** at rtol 10⁻¹⁰, against a swing energy of 2692 J (1.0 × 10⁻⁴ and 1.5 × 10⁻⁹ relative).

**Common mistakes.** Reading the 7.12 rad/s peak as a numerical error instead of checking the mode-2 frequency at that amplitude; normalising the energy error by the total energy (dominated by the large, constant −(m₁ + m₂)gl₁) rather than by the swing energy.

## 5. A generalized force: forward flight (s = 15: D = 3.5 kN)

**Assumptions.** Steady flight at constant velocity (an inertial suspension frame), drag only on the load, horizontal and constant.

(a) ∂r₂/∂θ₁ = (l₁ cos θ₁, l₁ sin θ₁), ∂r₂/∂θ₂ = (l₂ cos θ₂, l₂ sin θ₂). With D_vec = (−D, 0): **Q₁ = −D l₁ cos θ₁, Q₂ = −D l₂ cos θ₂**.

(b) At rest, (2) is ∂V/∂θᵢ = Qᵢ: (m₁ + m₂)gl₁ sin θ₁ = −Dl₁ cos θ₁ and m₂gl₂ sin θ₂ = −Dl₂ cos θ₂, so **tan θ₁ = −D/((m₁ + m₂)g), tan θ₂ = −D/(m₂g)**. Both links trail backwards; the sling trails more because only the load's own weight resists the drag.

(c) Mass 2: the sling tension along −e₂ balances m₂g down and D backwards: tan θ₂ = −D/(m₂g), T₂ = √((m₂g)² + D²). Mass 1: the cable balances m₁g plus the sling's pull, which is (m₂g down, D backwards): tan θ₁ = −D/((m₁ + m₂)g), T₁ = √(((m₁ + m₂)g)² + D²). Same angles.

(d) For the lesson's load (W = 15 494 N): **θ₁ = −12.73°, θ₂ = −13.38°, T₁ = 15 885 N = 1.025 W**, T₂ = 15 121 N. At these angles f(q, 0) + Q = (0, 0) to round-off (`generalized_forces` plus Q). `newton_euler` without the drag returns non-zero accelerations (0.0096 and 0.539 rad/s²) at the same state, as it should: the state is an equilibrium only with the drag.

**Common mistakes.** Taking Q₂ at the hook instead of at the load; forgetting that D also tilts the upper cable.

## 6. What can a simulation of a chaotic system be trusted for?

(a) Four decimals means an error Δ ≈ 10⁻⁴ rad at 30 s. With λ ≈ 1.4 s⁻¹, δ₀ = Δ e^{−λt} = 10⁻⁴ × e^{−42} ≈ **6 × 10⁻²³ rad**. Double precision carries about 16 significant digits (round-off near 10⁻¹⁶ rad for an angle of order 1), so the number is not reproducible by any double-precision simulation, let alone by a real pendulum whose release angle is known to perhaps 10⁻³ rad. Extrapolating the lesson's horizons (15.1 s from 10⁻¹⁰ rad, about 2.7 s per factor of 100) gives the same verdict. The reported value is noise.

(b) Two things that can be trusted: the **energy** (conserved to the solver tolerance, 10⁻⁷ J at rtol 10⁻¹⁰) and **statistical properties** of the motion, such as λ itself, the time-averaged kinetic energy, or the fraction of time the lower link spends above the pivot, which converge as the run lengthens and are insensitive to the exact trajectory.

(c) No. `lyapunov_estimate` for the s = 15 load released at θ₁ = 10°, θ₂ = 0 gives **0.101 s⁻¹ over 60 s and 0.026 s⁻¹ over 240 s**: the estimate falls roughly like (ln t)/t, the signature of regular motion with λ = 0. A chaotic motion would settle at a positive value.

(d) (1) The **hook mass** and what moves with it (the hook mode is a light mass between two taut lines, so ω₂ ∝ √(tension/m₁)): the cable's own mass, the hook's swivel and any rigging change the effective m₁; weigh them, or tap the cable in hover and measure ω₂. (2) The **fixed, rigid suspension and cable**: the helicopter's own motion and the cable's elasticity add degrees of freedom and shift the mode; a flight test with an accelerometer on the hook and on the airframe pins it down. Lesson 2.6's model gives the starting value; a notch must be designed with its uncertainty, and Problem 4 shows that the hook mode also shifts with amplitude.

**Common mistakes.** Treating energy conservation as proof that a chaotic trajectory is right; quoting a chaotic state at a long time without a horizon estimate.

## Check-yourself answers

1. Two (four Cartesian coordinates minus two constraints). The tensions do no work on any motion that the coordinates θ₁, θ₂ allow, because those motions keep the link lengths fixed, so they drop out of Hamilton's principle (or of d'Alembert's virtual work).
2. L = ½ml²θ̇² + mgl cos θ; d/dt(ml²θ̇) + mgl sin θ = 0, i.e. θ̈ + (g/l) sin θ = 0.
3. It is the inertial coupling between the links: the lower mass's velocity has a part due to the upper link's rotation, and their dot product carries cos Δ. It vanishes when the links are perpendicular (Δ = ±90°) and when m₂ = 0.
4. Because the Lagrangian does not depend explicitly on time and T is quadratic in the rates, so h = T + V is a first integral. A moving suspension point, friction at a joint, or aerodynamic drag would break it.
5. m₂ → 0: two independent pendulums, √(g/l₁) and √(g/l₂) (the upper mass does not feel the lower one). m₁ → 0: one root, √(g/(l₁ + l₂)) (the massless hook just lines up the two links into one pendulum).
6. The hook (m₁) moves most; the heavy load almost stays put. The hook is a small mass between two taut lines, so its frequency is roughly √((T₁/l₁ + T₂/l₂)/m₁), set by the hook mass and the tensions.
7. The average exponential rate λ at which two nearby motions separate, δ(t) ≈ δ₀e^{λt}. A simulation is useful until the error reaches a tolerance Δ, after t_h ≈ (1/λ) ln(Δ/δ₀): improving δ₀ helps only logarithmically.
8. No. For a small, regular motion an energy error of 10⁻¹⁰ goes with a small trajectory error that grows only slowly, so the trajectory is right. For a chaotic motion the trajectory error grows exponentially even when energy is held to 10⁻¹⁰ (experiment 4: lost after 18 s at rtol 10⁻¹²); energy conservation is necessary but not sufficient.
