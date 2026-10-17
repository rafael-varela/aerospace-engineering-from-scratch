# P06.03  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe script that calls it, given to the figures shown. All gear and aircraft numbers are illustrative.

## 1. Your aircraft on the lesson's gear (s = 15: m = 21 t, v₀ = 62.5 m/s)

**Assumptions.** Particle on a straight deck; no thrust, no wheel drag; force law (3).

(a) KE = **41.02 MJ**, p = **1.3125 MN·s**. Ideal 100 m: F = **410.2 kN**, a = 19.53 m/s² = **1.992 g**, t = 2 × 100/62.5 = **3.200 s**.

(b) Runout **100.17 m**, time **3.168 s**, peak **423.2 kN (2.055 g) at the stop**, efficiency **0.968**. The gear was set for a lighter, faster aircraft: for this one the c v₀² term at engagement (390.6 kN) is below the relaxation value km/(2c) = 443.6 kN, so the force rises along the stroke and peaks at the end.

(c) Distance average ½mv₀²/x_s = **409.5 kN** (work–energy); time average m v₀/t_s = **414.3 kN** (impulse–momentum). Both lie below the peak, as they must.

(d) `simulate_arrestment` stops at 100.168 m and 3.1680 s; the impulse is 1.3125 MN·s (relative error 1.2 × 10⁻⁹). Energy: damper **48.3 %** (19.82 MJ), k x term **51.7 %** stored (21.20 MJ).

**Independent check.** Three routes give the same runout and time to 10⁻⁹.

**Common mistakes.** Using the lesson's 65 m/s; reporting the peak at engagement out of habit.

## 2. Derivations

(a) m dv/dt = −F; multiply by v: m v dv/dt = −F v, and since v = dx/dt, m v dv/dx = −F, i.e. d(½mu)/dx = −F with u = v². Substituting F: (m/2)u′ = −F₀ − kx − cu, so **u′ + (2c/m)u = −(2/m)(F₀ + kx)**. Integrating factor e^{αx}, α = 2c/m: (u e^{αx})′ = −(2/m)(F₀ + kx)e^{αx}. A particular solution a + bx needs bα = −2k/m ⇒ b = −k/c, and b + αa = −2F₀/m ⇒ a = −F₀/c + km/(2c²). With u(0) = v₀²: **u = (v₀² − a)e^{−αx} + a + bx**. Substitution check as in 2.2 of the lesson.

(b) F(x) = F₀ + kx + cu(x) ⇒ F′ = k + cu′ = k + c(−αu − (2/m)(F₀ + kx)) = k − (2c/m)(F₀ + kx + cu) = **k − αF**. Solution F = km/(2c) + (F(0) − km/(2c))e^{−αx}: an exponential approach to km/(2c), so F is monotonic on [0, x_s] and its maximum is at x = 0 or x = x_s.

(c) F constant ⇔ F′ = 0 ⇔ **k = αF = 2cF/m**. With F = ½mv₀²/x_s: k = 2c(½mv₀²/x_s)/m and c = F/v₀² (F₀ = 0), so k = 2F²/(m v₀²) = **m v₀²/(2x_s²)**. A spring alone stops the aircraft in x_s when ½kx_s² = ½mv₀², i.e. k = mv₀²/x_s²: twice as stiff.

(d) Near the stop u(x) ≈ u′(x_s)(x − x_s) = |u′(x_s)|(x_s − x), so 1/v ~ (x_s − x)^{−1/2}, whose integral converges (∫ε^{−1/2}dε = 2√ε). With x = x_s(1 − s²), x_s − x = x_s s² and dx = −2x_s s ds, so dx/√u ≈ 2x_s s ds/(s√(|u′|x_s)): finite at s = 0.

**Common mistakes.** Treating the c v² term as a function of t; claiming an interior maximum because the force "first falls then rises" (it cannot).

## 3. Catching a drone on a bungee

**Assumptions.** Particle; bungee linear, massless; damper force always opposes motion.

(a) k = m v₀²/x² = 25 × 625/225 = **69.44 N/m**; peak force kx = **1041.7 N = 4.25 g** at the stop; time (π/2)√(m/k) = **0.942 s**. With no dissipation all 7.81 kJ is stored and returned: the drone passes the engagement point going backwards at **25.0 m/s** (simulated: 25.000 m/s at 1.885 s).

(b) With k = 34.72 N/m, `brentq` gives **c = 0.8333 N·s²/m²**, which is exactly `constant_force_design(25, 25, 15)`: halving the spring at the same runout *is* the constant-force condition (Problem 2c). Peak **2.12 g**, constant along the stroke.

(c) Integrating out and back with the sign-correct damper: the drone returns to x = 0 at **12.85 m/s** after 2.682 s, carrying **26.4 %** of its original kinetic energy. On the way out the damper dissipated 50 % and the bungee stored 50 %; on the way back the bungee's 50 % is again shared with the damper, which takes nearly half of it (23.6 % of the original).

(d) Stop the bungee giving energy back: a ratchet or one-way brake that locks the line at the stop, or replace most of the k x term with a dissipative brake (F₀). Either removes the rebound; the brake also keeps the peak at the constant-force value.

**Common mistakes.** Using c v² instead of c v|v| on the return (the damper would then push the drone backwards).

## 4. One gear, many aircraft (s = 15: m = 21 t, v₀ = 62.5 m/s)

(a)

| mass (t) | 57.5 m/s | 62.5 m/s | 67.5 m/s |
|---|---|---|---|
| 16.8 | 81.83 m, 2.099 g, stop | 87.00 m, 2.371 g, engagement | 91.98 m, 2.766 g, engagement |
| 21.0 | 94.08 m, 1.930 g, stop | 100.17 m, 2.055 g, stop | 106.03 m, 2.212 g, engagement |
| 25.2 | 105.30 m, 1.800 g, stop | 112.24 m, 1.919 g, stop | 118.93 m, 2.033 g, stop |

Light and fast aircraft peak at engagement (the damper term dominates), heavy and slow ones at the stop.

(b) The runout reaches 105 m at **22.64 t** at 62.5 m/s.

(c) For your aircraft at 100 m: **c = 105.0 N·s²/m², k = 4101.6 N/m**, against 100 and 4225. The damper is 5 % stronger (the same force at a lower speed needs a larger c) and the k term 3 % weaker.

## 5. Choosing a force law

(a) Lowest peak at design: the constant-force design or the pure brake (2.154 g, efficiency 1). Most tolerant of speed errors: the c v² + kx gear (a 70 m/s engagement overruns by 5.6 m against 16.0 m for the brake). Weight errors: neither; a 20 % heavier aircraft overruns by 12.1 m (c v² + kx) and 20 m (brake).

(b) At the stop the k x term holds ½kx_s², half the energy, and would pull the aircraft back up the deck. A real gear avoids it by making the absorber dissipative throughout (fluid through a controlled valve) and holding the cable at the end of the stroke; the "spring" behaviour then comes from a valve schedule, not from stored energy.

(c) The finite rise of the engagement force and the cable's wave dynamics: the damper-heavy design asks for its largest force at the very moment the cable is least able to deliver it smoothly, so the real peak on the hook would be a transient spike above 818.6 kN.

(d) Use a speed-sensitive, nearly constant-force law (a damper with a payout-scheduled term, as in the constant-force design), made dissipative so nothing is returned. Set it for each aircraft's weight before the landing, because no fixed law tolerates a ±20 % weight range within a 105 m deck. Accept speed variation without resetting, because the c v² term absorbs ±5 m/s within a few metres.

## Check-yourself answers

1. Work–energy gives the stopping distance and the force along the stroke (the distance average ½mv₀²/x_s); impulse–momentum gives the stopping time and the time average mv₀/t_s.
2. Because written in terms of x, v appears only as v²; the work–energy equation becomes linear in u = v².
3. Then v = v₀e^{−cx/m}: the force falls with v², so the deceleration vanishes as the aircraft slows and the speed never reaches zero in a finite distance.
4. Only at engagement or at the stop, because F′ = k − (2c/m)F makes F monotonic along the stroke.
5. A constant force; efficiency 1 (peak equals the distance average).
6. Distance average ½mv₀²/(100 m) from work–energy; time average mv₀/(4 s) from impulse–momentum.
7. Because the stored ½kx² would push the aircraft back up the deck once it stops.
8. Weight changes the energy in proportion and a fixed law overruns badly; speed variations are absorbed by the c v² term, which grows with v₀² and so adapts automatically.
