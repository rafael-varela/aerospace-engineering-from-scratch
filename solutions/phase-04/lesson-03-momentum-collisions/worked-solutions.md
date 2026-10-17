# P04.03  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for **s = 15**; numerical results are from `solution.py` and given to the figures shown.

## 1. A docking (s = 15: m = 9500 kg, v = 0.065 m/s)

**Assumptions.** One-dimensional, head-on contact; the station initially at rest in the frame used; no thrusters fire during contact.

(a) `collide_1d(9500, 0.065, 420e3, 0, 0)`: common velocity **1.438 mm/s**; energy absorbed **19.62 J** of the vehicle's 20.07 J.

(b) With e = 0.2: vehicle **−11.27 mm/s** (it rebounds), station **+1.725 mm/s**; energy absorbed **18.84 J**. The crew sees the station push them back at about a centimetre per second after contact: the capture latches must close during the brief contact or the vehicle drifts away and the approach is repeated.

(c) In the centre-of-mass frame (3) gives ΔK = ½μv² for e = 0, and the vehicle's kinetic energy is ½mv², so ΔK/K₀ = μ/m = M/(m + M) = 420/429.5 = **0.978**. "Within 3 %" requires M/(m + M) > 0.97, i.e. **M > 32.3 m**.

**Check.** Momentum: 9500 × 0.065 = 617.5 kg·m/s = 429 500 × 0.001438 to four figures. Limiting: M → ∞ gives v′ → 0 and ΔK → K₀.

**Common mistakes.** Using the closing speed in the station frame for an energy that must be computed in a single frame consistently; forgetting that a stiff mechanism (e > 0) makes the vehicle rebound.

## 2. From pellets to the rocket equation

(a) In the frame where the cart moved at v before the throw: momentum before mv; after, the cart (m − dm) at v + dv and the pellet at v + dv − u, since u is relative to the cart *after* the throw. Equating: mv = mv + m dv − dm u (the terms in v cancel; dm·dv cancels too), so **dv = u dm/m exactly**. If u is relative to the cart *before* the throw, the pellet moves at v − u: mv = (m − dm)(v + dv′) + dm(v − u) ⇒ dv′ = u dm/(m − dm) = (u dm/m)(1 + dm/m + …): the same to first order in dm, so the continuum limit is unaffected; for finite pellets the conventions differ (Problem 3c).

(b) With Δm = (m₀ − m_f)/n and mᵢ = m₀ − iΔm, v_n = uΣ_{i=0}^{n−1} Δm f(mᵢ) with f = 1/m: each term uses the value of f at the *upper* end of its interval [mᵢ₊₁, mᵢ] of mass, the end where f is smallest (call it the left sum in the order the throws happen). The trapezoid sum T = Σ Δm (f(mᵢ) + f(mᵢ₊₁))/2 has error O(Δm²) per interval, O(1/n²) overall. Left sum minus trapezoid = (Δm/2)Σ(f(mᵢ) − f(mᵢ₊₁)) = (Δm/2)(f(m₀) − f(m_f)), a telescoping sum. So u ln(m₀/m_f) − v_n = u(Δm/2)(1/m_f − 1/m₀) + O(1/n²) = **u(m₀ − m_f)(1/m_f − 1/m₀)/(2n)**. At n = 3: formula 1103.2 m/s, `pellet_convergence` **891.3 m/s**: the O(1/n²) term (which has the opposite sign) is 19 % of the total at n = 3; at n = 100 the formula is within 0.6 %, at n = 10 000 within 0.006 %.

(c) A rocket leaks propellant with zero relative speed: each parcel leaves at the rocket's own velocity v. The pellet balance (a) with u = 0 gives dv = 0: nothing changes, in every frame. The naive F_ext = d(mv)/dt = m dv/dt + v dm/dt with F_ext = 0 gives dv/dt = −(v/m)(dm/dt), which is non-zero whenever v ≠ 0 and has a magnitude set by v, i.e. by the observer's choice of frame. That violates Galilean invariance, so it is wrong. The correct statement applies momentum conservation to the *fixed* system "rocket at t plus the mass that leaves in dt": the momentum change in dt is F_ext dt = [(m − dm)(v + dv) + dm(v + dv − u)] − mv = m dv − u dm, hence **m dv/dt = ṁu + F_ext** with ṁ = dm/dt of the leaving mass.

(d) z_b = ∫₀^{t_b} [u ln(m₀/m(t)) − gt] dt. With m = m₀ − ṁt, dt = −dm/ṁ: ∫₀^{t_b} ln(m₀/m) dt = (1/ṁ)∫_{m_f}^{m₀} ln(m₀/m) dm = (1/ṁ)[m ln(m₀/m) + m]_{m_f}^{m₀} = (1/ṁ)(m₀ − m_f − m_f ln(m₀/m_f)) = t_b − (m_f/ṁ) ln(m₀/m_f). So **z_b = u t_b − u(m_f/ṁ) ln(m₀/m_f) − ½g t_b²**. For the lesson's case, 441 299 − 203 924 − 110 325 = 127 050 m, as the ODE gives.

**Common mistakes.** Measuring u relative to the ground; treating the left sum as an overestimate (f = 1/m decreases with m, but the throws run from large m to small, so each term uses the smaller value); differentiating (10) to check it and forgetting that m depends on t.

## 3. An astronaut without a tether (s = 15: tools of 3.5 kg, m₀ = 130.5 kg)

**Assumptions.** One dimension; throws are impulsive; the station is infinitely massive.

(a) All at once: 120v + 10.5(v − 5) = 0 ⇒ v = 52.5/130.5 = **0.4023 m/s**, 20 m in **49.7 s**. One at a time: `pellet_rocket_velocities(130.5, 120, 5, 3)[-1]` = **0.4136 m/s**, **48.4 s**.

(b) Tsiolkovsky: 5 ln(130.5/120) = **0.4194 m/s**. It is the upper bound: v_n increases with n toward it (each pellet thrown from a lighter body) and never exceeds it (2.3); the one-at-a-time strategy captures 98.6 % of it, all-at-once 95.9 %.

(c) Relative to the initial frame, momentum gives 120v = 10.5 × 5, **v = 0.4375 m/s**, larger because the tools are credited with 5 m/s relative to the original frame, not relative to her recoiling body: she would have to throw harder than 5 m/s relative to herself. A real throw is limited by how fast the arm can move relative to the body, so "relative to herself after the throw" is the realistic convention.

(d) Kinetic energy at 0.4136 m/s: ½ × 120 × 0.4136² = **10.3 J** < 30 J: **yes**, her arms can absorb it, with a safety factor of about 3.

**Check.** Momentum: after the one-at-a-time throws the tools carry −120 × 0.4136 = −49.6 kg·m/s in total. Units: (kg × m/s)/kg.

**Common mistakes.** Using 5 m/s relative to the initial frame for every throw; forgetting the tools' mass in m₀.

## 4. Choosing the thrust-to-weight ratio (s = 15: k = 0.289 kg/m)

(a)

| T/W₀ | burnout (m/s) | gravity loss (m/s) | drag loss (m/s) | max-q (kPa) |
|---|---|---|---|---|
| 1.2 | 2142.2 | 1838.7 | 97.5 | 17.8 |
| 1.5 | 2490.7 | 1471.0 | 116.8 | 27.9 |
| 2.0 | 2824.2 | 1103.2 | 151.1 | 47.3 |
| 2.5 | 3011.8 | 882.6 | 184.1 | 67.7 |
| 3.0 | 3126.8 | 735.5 | 216.2 | 88.7 |
| 4.0 | 3247.6 | 551.6 | 279.2 | 132.4 |
| 5.0 | 3296.5 | 441.3 | 340.7 | 178.0 |
| 6.0 | 3311.6 | 367.7 | 399.2 | 225.3 |

(b) Maximum burnout speed at **T/W₀ = 6.35**, **3312 m/s**: the gravity loss falls as 1/(T/W) while the drag loss grows, and they balance only at a very high thrust. A real vehicle cannot fly there because max-q would be about 240 kPa (structure) and the engines for six times the weight would be far heavier.

(c) Max-q = 40 kPa at **T/W₀ = 1.816**, burnout **2725 m/s**: the structural limit costs 3312 − 2725 = **587 m/s** against the unconstrained optimum.

(d) The model treats the engine as massless (in reality engine mass grows roughly with thrust, so a higher T/W eats the dry-mass budget), and it ignores acceleration limits (crew, payload) and the transonic rise of C_D, which raises the drag loss of a fast ascent. All three push the choice down toward the 1.2–1.5 of real first stages.

**Check.** With ṁ = (T/W₀)m₀g/u, the gravity loss is g t_b = g(m₀ − m_f)/ṁ = (m₀ − m_f)u/((T/W₀)m₀) = 7500 × 2942/(10 000 × T/W₀) = 2206.5/(T/W₀) m/s: 1838.7 at 1.2 and 1471.0 at 1.5, as tabulated.

**Common mistakes.** Optimising without the max-q constraint and reporting 6.35 as a design; using the drag-free burnout to compute q.

## 5. Staging with two different stages (s = 15: payload 350 kg)

(a) Best first-stage fraction **f = 0.914**, Δv = **10 185 m/s** (4934 + 5251); at f = 0.5, 8390 m/s.

(b) u₁ ln(1/ε₁) = 2800 ln 12.5 = **7072 m/s**; u₂ ln(1/ε₂) = 3400 ln(8.33) = **7209 m/s**. Each limit applies to a stage that must carry its own empty tanks to the end; staging drops stage 1's tanks after 4.9 km/s, so stage 2 starts its own logarithm afresh from a much lighter vehicle. The velocities add, and so the total can exceed either limit.

(c) The best split gives exactly 9400 m/s with a payload of **578.9 kg** (f = 0.889); gross lift-off mass = payload + 30 000 kg propellant + structures = **33 353 kg**, so **1.74 %** of the lift-off mass reaches the target velocity.

(d) Not always. In (11), the upper stage's mass ratio is capped by its own ε (and here ε₂ is the larger, so its tanks are relatively heavier), and every kilogram added to stage 2 is dead weight that stage 1 must accelerate through its whole burn. The optimum in (a) still puts 91 % of the propellant in the lower-u first stage; the higher exhaust speed is worth most where the vehicle is lightest, which is why high-u engines (hydrogen, vacuum nozzles) are upper-stage engines.

**Check.** The two-stage optimum exceeds both single-stage limits, as (b) says it can; the split is a maximum (neighbouring f give less).

**Common mistakes.** Using the total vehicle mass as m₀ for stage 2; forgetting stage 2's structure in stage 1's dead weight.

## 6. Momentum in practice

(a) m_p = M(1 − e^{−Δv/u}): cold gas **66.6 kg**, electric **1.66 kg**. The electric thruster's thrust is tiny (millinewtons per kilowatt of electrical power), so the manoeuvre takes days to weeks instead of seconds, and it needs a power system; for a quick collision-avoidance manoeuvre or a small satellite without spare power, the cold gas is the right choice.
(b) The mass falls during the burn, so the velocity change is u ln(m₀/m_f), not J/m₀; dividing the impulse by the initial mass underestimates Δv (by the factor ln R/(1 − 1/R)), and by the final mass overestimates it.
(c) e = v_up/v_down = √(h_up/h_down) = √0.64 = **0.8**; not elastic (36 % of the kinetic energy is lost at each bounce).
(d) A vent through a single nozzle throws momentum in one direction and changes the velocity (and so the orbit) by u ln(m₀/m_f); a non-propulsive vent with two opposite nozzles throws equal momenta both ways, so the net momentum given to the spacecraft is zero (only imperfect symmetry leaves a residual).
(e) No. Orbit at 127 km needs about 7.8 km/s *horizontally*; the vertical ascent has no horizontal speed at all and falls back. The missing pieces are a gravity turn and roughly three times the velocity.
(f) A 1 % drift means a modelling error, not a numerical one (the solver holds 10⁻¹¹): a thrust term not equal to ṁu, an exhaust momentum written as −ṁu instead of ṁ(v − u) (wrong frame), or the d(mv)/dt form of the equation of motion.

## Check-yourself answers

1. J = ∫F dt = Δp. For a half-sine of peak F₀ and duration τ, J = 2F₀τ/π = 12.7 N·s.
2. Momentum is conserved; kinetic energy is not. The energy of the relative motion, ½μv_rel², goes into deformation, heat and sound; in a docking it goes into the mechanism's dampers.
3. The energy lost is ½μv_rel², and μ ≈ m when M ≫ m (μ/m = 0.972 for 12 t on 420 t), so almost all of the vehicle's kinetic energy is relative-motion energy.
4. mv = (m − dm)(v + dv) + dm(v + dv − u) ⇒ m dv = u dm ⇒ dv = u dm/m.
5. Each pellet is thrown from a lighter vehicle, so it gives a larger velocity increment u dm/m; one large throw uses the heaviest mass for all of the propellant.
6. R = e¹ = 2.72; Δv = 3000 ln 10 = 6908 m/s.
7. A rocket is not a fixed set of particles; d(mv)/dt ignores the momentum carried away by the exhaust and gives a frame-dependent result (a leak with zero relative speed would change the velocity).
8. g t_b = 9.81 × 100 = 981 m/s. Higher thrust shortens the burn but drives the vehicle faster through the dense lower atmosphere, raising the drag loss and the dynamic pressure; the structure and the engine mass set the limit.
