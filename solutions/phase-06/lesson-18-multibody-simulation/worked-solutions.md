# P06.18  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and probe scripts that call it, given to the figures shown. All gear parameters are illustrative.

## 1. The gas spring (s = 15: M = 1400 kg)

**Assumptions.** Polytropic gas, the lesson's strut, rake 8°.

(a)

| n | F(0) [N] | F(0.15 m) [N] | F(0.30 m) [N] | U(0.30 m) [J] |
|---|---|---|---|---|
| 1.0 | 10 985 | 18 055 | 46 335 | 6294 |
| 1.1 | 10 985 | 18 961 | 53 335 | 6769 |
| 1.4 | 10 985 | 21 953 | 81 252 | 8492 |

Isothermal check: p₀V₀ ln(V₀/V) − p_atm A_a s = 1.5 × 10⁶ × 3.1416 × 10⁻³ × ln 4 − 101 325 × 7.854 × 10⁻³ × 0.3 = **6294 J** ✓. At full extension the force does not depend on n; at the end of stroke the adiabatic gas force is 75 % higher than the isothermal.

(b) Gas force = 1400 × 9.80665 × cos 8° = 13 596 N at **s = 66.5 mm**, where the stiffness is **dF/ds = 47.5 kN/m** (= n A_a² p/V).

(c) fₙ = (1/2π)√(47 475/1400) = **0.93 Hz**. A landing gear must absorb the aircraft's sink energy over a fraction of a metre of stroke without bottoming and give a soft ride on the ground, so it is soft; the IMU isolator must stay above the guidance bandwidth and carries only a 0.35 kg box. The gas spring's nonlinearity lets the gear be soft near static and stiff near the end of stroke.

**Common mistakes.** Using gauge pressure for p₀ in the polytropic law; forgetting the atmospheric pressure on the other side of the piston.

## 2. Derivations

(a) x₂ = −(ℓ₀ − s) sin β, y₂ = y₁ − (ℓ₀ − s) cos β, so ẋ₂ = ṡ sin β, ẏ₂ = ẏ₁ + ṡ cos β and T = ½Mẏ₁² + ½m(ṡ² sin²β + ẏ₁² + 2ẏ₁ṡ cos β + ṡ² cos²β) = ½(M + m)ẏ₁² + m cos β ẏ₁ṡ + ½mṡ², so **M_r = [[M + m, m cos β], [m cos β, m]]**. Virtual work: gravity −Mg δy₁ − mg δy₂, lift L δy₁, tyre F_t δy₂, strut −F δs, with δy₂ = δy₁ + cos β δs. Collecting: **Q_y₁ = −(M + m)g + L + F_t**, **Q_s = (F_t − mg) cos β − F**. The strut force is internal: a δy₁ moves both of its ends together and does not change its length, so it does no virtual work on y₁. The joint and guide forces act perpendicular to every allowed virtual displacement and do no virtual work at all.

(b) U = ∫₀ˢ A_a p₀V₀ⁿ(V₀ − A_a s′)⁻ⁿ ds′ − p_atm A_a s = p₀V₀ⁿ[(V₀ − A_a s)^{1−n} − V₀^{1−n}]/(n − 1) − p_atm A_a s = **p₀V₀/(n − 1)[(V₀/V)^{n−1} − 1] − p_atm A_a s**. As n → 1, (x^{n−1} − 1)/(n − 1) → ln x, giving **p₀V₀ ln(V₀/V) − p_atm A_a s**.

(c) Φ = (x₁, θ₁, θ₂ − θ₁, d·n(θ₁)), d = r₂ − r₁. Rows of Φ_q: (1, 0, 0, 0, 0, 0); (0, 0, 1, 0, 0, 0); (0, 0, −1, 0, 0, 1); (−n_x, −n_y, d·n′, n_x, n_y, 0), with n′ = dn/dθ₁. Differentiating Φ₄ twice with ṅ = θ̇₁n′ and ṅ′ = −θ̇₁n: Φ̈₄ = d̈·n + θ̈₁ d·n′ + 2θ̇₁ ḋ·n′ − θ̇₁² d·n. The first two terms are Φ_q q̈, so **γ₄ = −2θ̇₁ ḋ·n′ + θ̇₁² d·n** (zero along this motion, where θ̇₁ = 0).

(d) The power of the constraint forces is q̇ᵀΦ_qᵀλ = λᵀ(Φ_q q̇) = λᵀΦ̇ = 0 for any motion that satisfies the constraints. So the Newton–Euler model has the same energy balance as the Lagrange model: only applied forces appear. If the constraints drift, Φ̇ ≠ 0 and the constraint forces start doing spurious work, which is one way drift shows up in an energy check.

**Common mistakes.** Including the strut force in Q_y₁; dropping the θ̈₁ d·n′ term from the Jacobian.

## 3. A crushable lander leg

**Assumptions.** Vertical motion, no lift, equal legs, constant crushing force over the full 0.15 m.

(a) E = ½ × 300 × 2.5² + 300 × 3.711 × 0.15 = 937.5 + 167.0 = **1104.5 J**. Constant force over 0.15 m: **7363 N in total, 1841 N per leg**.

(b) a = F/m − g_Mars = 24.54 − 3.71 = 20.83 m/s² = **2.12 Earth g** of deceleration; load factor F/(m g₀) = **2.50**.

(c) A constant force over the stroke is the ideal: **η = 1**. A linear spring absorbing 1104.5 J over 0.15 m needs a peak force of 2E/x = 14 727 N: **twice** the crushable's. Crushables are light, simple and ideal for a single landing with no rebound (the energy goes into plastic deformation); aircraft land thousands of times and need a reusable absorber whose force adapts to the sink rate, which the orifice does.

**Common mistakes.** Forgetting the work of gravity over the stroke; quoting deceleration and load factor as the same number (they differ by g_Mars/g₀ in the normalisation).

## 4. Choose the orifice (s = 15: v = 3.25 m/s)

(a) On the 0.1 cm² grid the minimum is at **A_o = 4.2 cm²**: load factor **2.230**, stroke **245.8 mm**, efficiency **0.898**, max energy residual **1.1 × 10⁻⁷ J**, no bottoming. The lesson's 4.0 cm² orifice gives **2.282** and **242.1 mm**: 2.3 % more load for 4 mm less stroke.

(b) The curve has a V shape (the lesson's experiment 2 at 3.05 m/s shows it): to the left the orifice force dominates and peaks early when the stroke rate is highest; to the right the orifice cannot stop the aircraft quickly and the gas spring, stiffening towards the end of stroke, sets a late peak. The minimum is where the two peaks are about equal, which is also near the highest efficiency.

(c) Not exactly. C_d is uncertain by about ±0.1, which is ±14 % in the effective area C_d A_o; the left side of the curve is steeper than the right, so a design slightly to the right of the minimum (a slightly larger orifice) loses little if C_d is high and is protected if C_d is low. Real gears go further and use a metering pin.

**Common mistakes.** Choosing the orifice by stroke alone; comparing runs whose energy budgets were not checked.

## 5. When does it bottom?

(a) Bottoming at 0.30 m: **4.48 m/s with lift = weight** (load factor 3.63), **3.02 m/s with no lift** (load factor 3.62).

(b) At bottoming the stored energy is almost the same in both cases (8648 J against 8642 J: gas at 300 mm plus tyre at 68 mm). With lift, the input is almost all kinetic: 15.35 kJ, of which 6.70 kJ is dissipated in the orifice. Without lift, gravity adds 5.42 kJ of work over the stroke and tyre deflection, so a much smaller kinetic energy (7.04 kJ) reaches the same stroke; and the slower sink makes a smaller orifice force, so only 3.82 kJ is dissipated. The limits differ by more than the gravity work alone suggests because the orifice's share of the work depends on the speed.

(c) In sink rate, 4.48/3.05 = **1.47**; in energy, (4.48/3.05)² = **2.16**. The energy ratio is the honest margin for an energy absorber (and the lift assumption must be stated with it: without lift the strut bottoms at 3.02 m/s, *below* the requirement).

**Common mistakes.** Looking for bottoming with a run too short to reach the peak stroke; treating the margin in velocity as the margin in energy.

## 6. The drop test disagrees

(a) Peak tyre force and its time after touchdown: **C_d = 0.6**: 35.9 kN at 40.0 ms; **k_t = 1.2 × 10⁶ N/m**: 34.5 kN at 21.5 ms; **k_t = 1.6 × 10⁶ N/m**: 36.6 kN at 15.7 ms; baseline 32.2 kN at 33.4 ms.

(b) The **tyre stiffness**. A stiffer tyre transmits the impact sooner, before the strut has built up its stroke rate, and k_t ≈ 1.6 × 10⁶ N/m (twice the assumed value) reproduces both the force and the time. A lower C_d also matches the peak force (C_d = 0.599 gives 36.0 kN) but moves the peak *later*, to 40 ms, the opposite of what was measured: the timing rules it out as the sole explanation.

(c) (1) The tyre deflection against tyre force, from a separate static tyre test or from the axle height and ground force in the drop: it should show the stiffer slope directly. (2) The stroke and stroke-rate records: if the tyre is the cause, the stroke history should be close to the model's (229 mm peak stroke against 230.6 mm with the stiff tyre), whereas a low C_d would show a shorter stroke (217 mm).

(d) **No**, not to match one number. Tuning C_d to the peak force would give the right peak at the wrong time and a wrong stroke, and would make every other prediction worse. Fix the parameter the evidence points to (measure the tyre), and keep C_d at its measured or handbook value; if a residual difference remains, report it.

**Common mistakes.** Using only the peak force to compare model and test; tuning the most convenient parameter.

## Check-yourself answers

1. Gas spring (pneumatic): grows with stroke, polytropically, A_a[p₀(V₀/V)ⁿ − p_atm]. Orifice (hydraulic): grows with the square of the stroke rate and opposes it. Top-out stop: acts only beyond full extension and holds the strut extended.
2. A virtual change of y₁ moves both ends of the strut together without changing its length, so the strut's internal force does no virtual work on y₁.
3. A multiplier is the force (or moment) a joint exerts to enforce its constraint; λ₄, conjugate to the constraint (r₂ − r₁)·n = 0, is the force across the sliding joint, the bearing side load.
4. Along this motion θ₁ stays exactly zero, so all constraints are linear in q, and explicit Runge–Kutta methods preserve linear invariants exactly. It would not with nonlinear constraints, for example an airframe free to pitch or a revolute joint.
5. Φ̈ + 2αΦ̇ + β_B²Φ = 0: the gains make any constraint violation decay like a damped oscillator instead of persisting or growing linearly.
6. The work done by the strut force over the compression stroke divided by the peak force times the stroke; 0.5 for a linear spring and 1 for a constant force.
7. With lift = weight, essentially the kinetic energy ½(M + m)v²; without lift, also the work of the whole weight over the stroke and the tyre deflection, which made the lesson's strut bottom.
8. The orifice area (experiment 2); the hydraulic force depends on C_d and A_o only through their product C_d A_o, so changing C_d is exactly equivalent to changing the area.
