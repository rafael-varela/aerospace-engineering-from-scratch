# P07.02  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15; numbers come from `solution.py` and the probe scripts that call it. N₂: R = 296.803 J/(kg·K), cv = 742.008, cp = 1038.811 J/(kg·K). He: R = 2077.264, cv = 3115.897, cp = 5193.161 J/(kg·K). All results are computed; no pressure vessel is built, filled or vented at home in this course.

## 1. Closed systems (s = 15: T₁ = 325 K)

**Assumptions.** Ideal gas, constant cv, quasi-static processes.

(a) v₁ = R T₁/p₁ = 0.096461 m³/kg; v₂ = 3v₁.

| path | p₂ (kPa) | T₂ (K) | w (kJ/kg) | Δu (kJ/kg) | q (kJ/kg) |
|---|---|---|---|---|---|
| n = 1 | 333.33 | 325.00 | 105.973 | 0 | 105.973 |
| n = 1.2 | 267.58 | 260.89 | 95.139 | −47.569 | 47.569 |
| n = 1.4 | 214.80 | 209.43 | 85.755 | −85.755 | 0 |

The end volume, and so v₂, is the same for every path; p₂, T₂, Δu, w and q are not. Δu depends only on the end states (u is a property), but the end states differ because the paths end at different pressures. Notice that w per kilogram does not depend on p₁ at all: p₁v₁ = R T₁, so only T₁ and the volume ratio enter. (With T₁ = 300 K these are the lesson's 97.82, 87.82 and 79.16 kJ/kg, scaled by 325/300.)

(b) m = p₀V/(R T₀) = 20 × 10⁶ × 0.01/(2077.26 × 300) = **0.32093 kg**. Rigid and closed: p₂ = p₀ T₂/T₀ = **21.0 MPa**. Q = m cv ΔT = **15.0 kJ**. Proof: m cv = (p₀V/(R T₀)) × R/(γ − 1) = p₀V/((γ − 1)T₀) = 20 × 10⁶ × 0.01/(0.6667 × 300) = 1000 J/K; times 15 K = 15 kJ.

**Check.** For n = γ, q = 0 to round-off: the adiabatic path is the polytropic one with n = γ.

**Common mistakes.** Using the isothermal formula with n = 1.2; taking Δu = 0 for every expansion "because it is an ideal gas" (Δu = 0 only if T does not change).

## 2. Derivations

(a) CV = the inside of the rigid tank. No shaft or boundary work (rigid); one exit, no inlet; uniform state, so the exit gas carries the tank's h; kinetic and potential energy of the contents negligible, so E_cv = m u. Energy: d(mu)/dt = Q̇ − ṁh. Mass: dm/dt = −ṁ. Expanding, m du/dt + u(−ṁ) = Q̇ − ṁh, so m du/dt = −ṁ(h − u) + Q̇. With du = cv dT and h − u = p v = R T: **m cv dT/dt = −ṁ R T + Q̇**.

(b) Energy along the isentropic streamline from the tank (at rest): cp T₀ = cp T + V²/2, so T₀/T = 1 + (γ − 1)M²/2 with M = V/√(γRT). Isentropic: ρ/ρ₀ = (T/T₀)^(1/(γ−1)). The mass flux is ρV = ρ₀(T/T₀)^(1/(γ−1)) M √(γ R T); at M = 1, T* = 2T₀/(γ + 1), so ρ*V* = ρ₀ (2/(γ + 1))^(1/(γ−1)) √(γ R T₀) (2/(γ + 1))^(1/2) = (p₀/(R T₀)) √(γ R T₀) Γ = **Γ p₀ √(γ/(R T₀))**, with Γ = (2/(γ + 1))^((γ+1)/(2(γ−1))). The throat reaches M = 1 when p_back ≤ p* = **p₀(2/(γ + 1))^(γ/(γ−1))**. Multiply by Cd A.

(c) With Q̇ = 0, every element of gas still in the tank has undergone only quasi-static expansion with no heat and no dissipation (the dissipation happens in the jet, outside the CV), so its entropy is constant and p ∝ ρ^γ, T ∝ ρ^(γ−1). Then V dρ/dt = −Cd A Γ p √(γ/(R T)) = −Cd A Γ √(γ R T₀) ρ₀ r^((γ+1)/2), with r = ρ/ρ₀; that is dr/dt = −r^((γ+1)/2)/τ with τ = V/(Cd A Γ a₀). Separating: ∫ r^(−(γ+1)/2) dr = −t/τ ⇒ (2/(γ − 1))(r^(−(γ−1)/2) − 1) = t/τ ⇒ r = (1 + (γ − 1)t/(2τ))^(−2/(γ−1)) and **p/p₀ = r^γ = (1 + (γ − 1)t/(2τ))^(−2γ/(γ−1))**. And ṁ₀ = Cd A Γ ρ₀ a₀, so m₀/ṁ₀ = Vρ₀/(Cd A Γ ρ₀ a₀) = **τ**.

(d) d(U + H_out − Q_in)/dt = (−ṁh + Q̇) + ṁh − Q̇ = 0 for every state, so the combination is a linear invariant. A Runge–Kutta step forms y_{n+1} = y_n + Δt Σ bᵢ kᵢ, where every stage value kᵢ is a right-hand-side evaluation and so satisfies cᵀkᵢ = 0 for c = (0, 1, 1, −1); hence cᵀy_{n+1} = cᵀy_n exactly, apart from floating-point rounding. Agreement with (c) depends on how well the nonlinear part (m, U) is integrated, which is the truncation error controlled by rtol.

**Common mistakes.** Using u instead of h for the exiting gas (that drops the flow work and gives no cooling at all); treating the whole tank's gas as having lost entropy (the irreversibility is in the jet, outside the CV).

## 3. Filling a tank

**Assumptions.** Adiabatic, rigid tank, uniform state, constant line state, ideal gas with constant cv.

(a) With m₁ = 0: m₂u₂ = h_line m₂, so cv T₂ = cp T_line and **T₂ = γ T_line**, for any final pressure. N₂: **420 K**; He: **500 K**. The flow work done by the line on the incoming gas becomes internal energy.

(b) m = 30 × 10⁶ × 0.05/(296.803 × 420) = **12.033 kg**. Cooled at constant volume to 300 K: p = 30 × 300/420 = **21.43 MPa**. An isothermal fill to 30 MPa holds 16.846 kg, **4.81 kg (40 %) more**.

(c) m₁ = 10 × 10⁶ × 0.05/(296.803 × 300) = 5.6154 kg. Energy: cv T₂ m₂ − cv T₁ m₁ = cp T_line (m₂ − m₁), with m₂ = p₂V/(R T₂), so cv p₂V/R − cv T₁ m₁ = cp T_line p₂V/(R T₂) − cp T_line m₁. Solving (the probe used `brentq`): **T₂ = 370.59 K**, **m₂ = 13.637 kg**.

(d) Filling heats the gas; the hot gas reads a high pressure that falls as it cools, so the bottle is under-filled unless it is topped up after cooling, and the composite and liner have temperature limits. Slow or staged fills let the wall carry the heat away.

**Check.** The 16.Unified notes' "Transient filling of a tank" (2.3.3) reaches the same T = γ T_line from the same balance. Limiting case: γ → 1 gives no heating.

**Common mistakes.** Writing the balance with u_line instead of h_line; assuming the final temperature depends on the final pressure.

## 4. Your orifice (s = 15: d = 2.25 mm)

(a) A = 3.976 × 10⁻⁶ m²; ṁ₀ = **218.97 g/s**; τ = **76.933 s**. Adiabatic: t = (2τ/0.4)(0.1^(−1/7) − 1) = 384.66 × 0.3895 = **149.82 s**. Isothermal: τ ln 10 = **177.14 s**.

(b)

| mode | t to 3 MPa (s) | T (K) | m left (kg) | energy residual |
|---|---|---|---|---|
| adiabatic | 149.82 | 155.38 | 3.2525 | 1.2 × 10⁻¹⁶ |
| isothermal | 177.14 | 300.00 | 1.6846 | 2.5 × 10⁻¹⁶ |

The event times agree with (a) to the integrator's tolerance (149.824 368 against 149.824 368).

(c) p = 3 × 300/155.38 = **5.79 MPa**.

(d) τ(2.25 mm)/τ(2 mm) = 76.93/97.37 = **0.790 = (2/2.25)²**: the area ratio, as expected, because τ = V/(Cd A Γ a₀). The end states at 3 MPa did not change: they depend only on p/p₀ and γ.

**Common mistakes.** Using the diameter instead of the area; expecting a bigger orifice to strand less gas in the adiabatic limit.

## 5. Nitrogen against helium

(a)

| gas | hA (W/K) | t to 3 MPa (s) | T (K) | m left (kg) |
|---|---|---|---|---|
| N₂ | 0 | 189.62 | 155.38 | 3.2525 |
| N₂ | 10 | 205.71 | 183.16 | 2.7592 |
| N₂ | 100 | 232.20 | 282.61 | 1.7883 |
| He | 0 | 60.89 | 119.43 | 0.6046 |
| He | 10 | 65.13 | 132.38 | 0.5455 |
| He | 100 | 84.20 | 250.11 | 0.2887 |

(b) m₀cv = p₀V/((γ − 1)T₀): 12 500 J/K for N₂, 7500 J/K for He. Thermal time constants: N₂ **1250 s** (10 W/K) and **125 s** (100 W/K); He **750 s** and **75 s** (hA = 0: infinite). Crossover hA = m₀cv/τ: N₂ 12 500/97.37 = **128.4 W/K**; He 7500/34.70 = **216.1 W/K**.

(c) Helium blows down 2.8 times faster (τ = 34.7 s against 97.4 s) because its speed of sound is 2.89 times higher, while its heat capacity per tank is only 40 % lower. The wall has less time per joule of cooling to act, so helium needs a larger hA (216 against 128 W/K) to reach the same fraction of its isothermal limit; at 100 W/K nitrogen ends 17 K below the wall temperature, helium 50 K below.

(d) The stage carries the bottle, not just the gas: per unit of ullage volume pressurised, helium's low molar mass means a much lighter gas load for the same pressure, which outweighs its larger cooling. And helium does not condense or dissolve appreciably in cryogenic propellants (liquid oxygen, liquid hydrogen), where nitrogen would liquefy or dissolve; its very low boiling point makes it the only practical pressurant for liquid hydrogen.

**Check.** Each row lies between its gas's adiabatic and isothermal limits in T and m.

**Common mistakes.** Comparing the gases at equal mass instead of at equal tank; reading the time column as monotonic in hA (it is not, for nitrogen at large hA: lesson, Experiment 1).

## 6. Judgment

(a) **The equation of state.** At 30 MPa the ideal gas mispredicts the mass stored by more than 1 %; use a real-gas property table or equation of state (lesson 03's Z method), and verify by weighing the bottle.

(b) **The uniform state.** In a 5 s vent the gas near the outlet and the core cool differently and the wall barely responds; the well-mixed adiabatic model gives a reasonable average but not the outlet temperature, which needs a measured or computed stratification (CFD, or test data).

(c) **Fine, near the isothermal limit.** An hour-long draw is far slower than any thermal time constant here, so the isothermal model applies; the remaining issue is the real-gas mass at high pressure.

(d) **The single-phase ideal gas.** The model's gas reaches 93 K at 0.5 MPa, at nitrogen's condensation temperature; below that the "gas" would start to condense and the model is wrong; in reality wall heat keeps the gas warmer, which is one reason the adiabatic limit is never reached.

(e) **Conservative for the gas temperature, not for everything.** The adiabatic gas temperature is the coldest the gas can get for a given pressure ratio, so it is a sound worst case for a cold seal; but the time to a pressure is not bounded by either limit (Experiment 1), and a seal that sees a rapid temperature swing also needs the rate, which the adiabatic model overstates.

## Check-yourself answers

1. ΔU = Q − W, with Q positive into the system and W positive when done by the system.
2. A flowing stream carries u and also does (or has done on it) flow work p v per kilogram as it crosses the boundary; u + p v = h.
3. m cv dT/dt = −ṁ R T + Q̇. The remaining gas does flow work p v on the gas it pushes out, at the expense of its internal energy.
4. 0.5283 (p_back/p₀). The throat is sonic, so the mass flow depends only on the tank's p₀ and T₀ (ṁ ∝ p₀/√T₀), not on the back pressure.
5. Its pressure falls both because mass leaves and because the remaining gas cools; it reaches the limit pressure while the gas is still cold and dense, so more mass is left inside.
6. p/p₀ = exp(−t/τ), τ = V/(Cd A Γ a₀) = m₀/ṁ₀: the time the initial mass flow would take to empty the tank.
7. U + ∫ṁh dt − ∫Q̇ dt has zero time derivative identically; Runge–Kutta steps preserve linear invariants exactly, so only rounding remains.
8. The gas warms at constant volume and the pressure rises, to 5.79 MPa from 3 MPa in the lesson's case: p ∝ T at fixed m and V.
