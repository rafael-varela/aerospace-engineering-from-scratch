# P04.02  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for **s = 15**; numerical results are from `solution.py` with the lesson's trainer (5000 kg, S = 20 m², C_D0 = 0.022, k = 0.060, T₀ = 16 kN, thrust ∝ ρ, M ≤ 0.75) and given to the figures shown.

## 1. The energy accounts (s = 15: h = 3500 m, V = 135 m/s)

**Assumptions.** Point mass, flat Earth, constant g; in (b) no thrust and no drag during the zoom.

(a) h_E = 3500 + 135²/(2 × 9.807) = 3500 + 929.2 = **4429.2 m**. K = ½ × 5000 × 135² = **45.56 MJ**; U = mgh = 5000 × 9.807 × 3500 = **171.62 MJ**. The altitude account holds almost four times the speed account; it is the speed account that the pilot can change quickly.

(b) The top altitude h_t satisfies h_t = h_E − (1.2V_s(h_t))²/2g with V_s(h) = √(2W/(ρ(h)SC_Lmax)). `brentq` gives **h_t = 4111.6 m** with V_top = 78.93 m/s (1.2 × 65.78): a **611.6 m** gain.

(c) P_s(3500 m, 135 m/s) = **18.99 m/s**. With dV/dh = 0.02 s⁻¹: ḣ = 18.99/(1 + 135 × 0.02/9.807) = 18.99/1.2753 = **14.89 m/s**.

**Check.** Without the stall constraint a zoom to V = 0 would reach h_E = 4429 m; the answer is below it by (78.93)²/2g = 317.6 m, as it must be. Dimensional: V/g × dV/dh is (m/s)/(m/s²) × 1/s, dimensionless.

**Common mistakes.** Using the sea-level stall speed at the top (the stall speed rises with altitude, 63.7 m/s at 3500 m against 65.8 at 4112 m); forgetting that a zoom conserves h_E, not h.

## 2. The special speeds and the phugoid

(a) The forces are T and D along the velocity, L perpendicular to it, and W. The power of the net force is (T − D)V + L·v + W·v; L·v = 0 because lift is *defined* perpendicular to the flight path, and W·v = −mg ḣ = −d(mgh)/dt because gravity is conservative. So d(½mV²)/dt = (T − D)V − d(mgh)/dt, i.e. d(½mV² + mgh)/dt = (T − D)V. Dividing by W = mg: **d(h + V²/2g)/dt = (T − D)V/W = P_s**.

(b) With q = ½ρV², D(q) = qSC_D0 + kW²/(qS). dD/dq = SC_D0 − kW²/(q²S) = 0 ⇒ q²S²C_D0 = kW² ⇒ (W/(qS))² = C_D0/k ⇒ **C_D0 = kC_L²**: the zero-lift and induced drags are equal. Then C_D = 2C_D0, C_L = √(C_D0/k) and L/D = C_L/C_D = √(C_D0/k)/(2C_D0) = **1/(2√(kC_D0))**; in level flight L = W so **D_min = W/(L/D)max**. Trainer: 13.76, 3.56 kN.

(c) With q ∝ V², P = DV = aV³ + b/V with a = ½ρSC_D0 and b = 2kW²/(ρS). dP/dV = 3aV² − b/V² = 0 ⇒ 3aV⁴ = b ⇒ 3 × qSC_D0 = kW²/(qS) ⇒ **3C_D0 = kC_L²**. Since C_L ∝ 1/V² at fixed W, V_mp²/V_md² = C_L,md/C_L,mp = √(C_D0/k)/√(3C_D0/k) = 1/√3, so **V_mp = V_md/3^{1/4}** (61.8 against 81.3 m/s at sea level). There C_D = C_D0 + 3C_D0 = 4C_D0 and L/D = √(3C_D0/k)/(4C_D0) = (√3/2) × 1/(2√(kC_D0)) = **(√3/2)(L/D)max** = 11.92 for the trainer (the code: 11.918).

(d) Linearise about level flight at V₀, h₀ with γ small. Energy (no thrust, no drag): ½V² + gh = const ⇒ **V₀δV = −g δh**. Lift at constant C_L: L = ½ρV²SC_L, so δL/L₀ = δρ/ρ + 2δV/V₀ with L₀ = W; with δρ/ρ = −βδh (linear in δh), **δL = W(−βδh + 2δV/V₀)**. Vertical motion: m d²(δh)/dt² = L cos γ − W ≈ δL (cos γ ≈ 1 to first order). Substituting δV = −gδh/V₀: m δḧ = W(−β − 2g/V₀²)δh, and with W/m = g, **δḧ = −(2g²/V₀² + gβ)δh**. At 150 m/s and 3 km, β = 1.0297 × 10⁻⁴ m⁻¹: ω² = 8.548 × 10⁻³ + 1.010 × 10⁻³ = 9.558 × 10⁻³ s⁻², **T = 64.27 s**; `simulate_flight` with the ISA gives **64.33 s**; at constant density the simulation gives 67.96 s = √2π × 150/9.807.

**Check.** (b) and (c): the numerical minima of D and DV in the tests land on these speeds to 0.02 m/s. (d): β → 0 recovers Lanchester.

**Common mistakes.** In (c), setting dD/dV = 0 instead of d(DV)/dV = 0; in (d), forgetting the density term or the factor 2 from V².

## 3. Storing solar energy as altitude (s = 15: m = 55 kg)

**Assumptions.** The 150 W is power available beyond what level flight needs at the start; no change in solar input or propeller efficiency with altitude; quasi-steady climb.

(a) W = 55 × 9.807 = 539.4 N; ḣ = 150/539.4 = **0.278 m/s**; in 6 h (21 600 s), **Δh = 6007 m**, from 1000 to 7007 m. A fair *first* estimate because the climb rate is small and the change in power required is a correction; but see (c): the minimum power required grows from 372 W to 511 W between 1000 and 7007 m (P_r ∝ 1/√ρ at fixed C_L), which would eat 139 W of the 150 W spare, so the real gain is much smaller. The estimate is an upper bound.

(b) E = 150 × 21 600 = **3.24 MJ = 900 Wh**; at 250 Wh/kg, **3.6 kg** of cells.

(c) (L/D)max = 23.57. At 1000 m: V_mp = **14.07 m/s** (stall 13.66: only 3 % above, too close for a real aircraft in gusts), P_r = 371.9 W, **sink 0.689 m/s**. At 7007 m: V_mp = 19.33 m/s (stall 18.77), P_r = 510.9 W, **sink 0.947 m/s**. Integrating dt = dh/sink(h) from 7007 down to 1000 m: **2.08 h**.

(d) No, not on these numbers: the stored altitude carries the aircraft through about 2 hours of a 12-hour night, the climb estimate of (a) is optimistic by most of the spare power, and the 3.6 kg of batteries it replaces is a small fraction of 55 kg. To decide you need the night length, the solar power curve through the day, the battery mass actually carried, and the power required at the altitudes considered; the idea is used in practice as a supplement to batteries, not a replacement.

**Check.** Energy: WΔh = 539.4 × 6007 = 3.24 MJ = PΔt, consistent. Glide time against a crude estimate Δh/mean sink = 6007/0.82 = 7300 s = 2.0 h.

**Common mistakes.** Treating the 150 W as total power rather than excess power; flying the night glide at V_md (sink 0.786 m/s at 18.5 m/s: V_md minimises the sink *per metre flown*, V_mp the sink *per second*, and endurance is per second).

## 4. The minimum-time climb (s = 15: h_E = 10 500 m)

(a) From h_E = 1147 m to 10 500 m: **643.4 s**, ending at **V = 168.1 m/s, h = 9059 m** (P_s = 7.09 m/s there).

(b) Accelerate-then-constant: **130 m/s: 735.9 s** (+14 %), **150 m/s: 654.6 s** (+1.7 %), **160 m/s: 645.1 s** (+0.3 %). Too slow a constant speed sits far below the P_s ridge at every altitude; 150 and 160 m/s straddle the optimal schedule (153–168 m/s), and the P_s curve is flat near its maximum (a quadratic: a 5 % speed error costs a fraction of a per cent in P_s), so they lose little.

(c) ∫₀^{9059} dh/P_s,max(h) = **611.6 s**, shorter than (a). It is wrong because it is the time to reach 9059 m at the best-climb speed of each altitude, which ignores the kinetic energy: the aircraft starts at 150 m/s and must end at 168 m/s, and the best-climb speeds rise along the climb, so energy must be spent on acceleration. The energy method charges for it; the steady estimate does not (equivalently, it ignores the acceleration factor).

(d) With 1.1 × mass: **746.0 s**, **+16 %**. The mass enters P_s twice: the induced drag grows as W² (through C_L), and the excess thrust is divided by a larger W. The excess T − D shrinks while the divisor grows, so P_s falls by more than 10 %, most of all near the ceiling where T − D is small.

**Check.** (a) ≤ every time in (b), as the optimality argument requires; (c) is not a schedule from the same initial and final *energy* states, so it does not contradict optimality.

**Common mistakes.** Comparing (a) with a climb that ends at a different speed; integrating dh rather than dh_E.

## 5. What a hard turn costs (5 km, 150 m/s)

(a)

| n | bank | P_s (m/s) |
|---|---|---|
| 1 | 0° | 15.60 |
| 1.5 | 48.2° | 12.27 |
| 2 | 60.0° | 7.61 |
| 3 | 70.5° | −5.71 |

The induced drag is qSkC_L² with C_L = nW/(qS), so it grows as n²; the zero-lift drag does not change. P_s therefore falls quadratically in n, and between n = 2 and 3 it changes sign.

(b) P_s(n) = 0 at **n = 2.618** (bank **67.5°**); turn rate g√(n² − 1)/V = 9.807 × 2.420/150 = 0.1582 rad/s = **9.06°/s**; radius V²/(g√(n² − 1)) = **948 m**.

(c) The sustained load factor is largest at **172.3 m/s** with **n = 2.70**: a faster aircraft has more dynamic pressure to make lift cheaply (lower C_L, lower induced drag) until the zero-lift drag growing as V² takes over.

(d) q = ½ × 0.7364 × 150² = 8285 Pa, n_max = qSC_Lmax/W = 8285 × 20 × 1.4/49 033 = **4.73**. At n = 4, P_s = −24.37 m/s; in 5 s, **Δh_E = −121.8 m**. Holding altitude, all of it comes from speed: √(150² − 2g × 121.8) = 141.8 m/s, a loss of **8.2 m/s** (the linear estimate gΔh_E/V gives 8.0).

**Check.** At n = 1, (a) reproduces the P_s map at 5 km; the sustained n lies between the n = 2 and 3 rows where P_s changes sign.

**Common mistakes.** Confusing the instantaneous (C_Lmax-limited) and sustained (P_s = 0) load factors; using the sea-level density at 5 km.

## 6. Trusting the energy method

(a) **True**: h_E = h + V²/2g contains no mass, so the ideal zoom altitude is the same for any aircraft at the same h and V; the mass matters only through drag during the real, non-ideal zoom.
(b) **False in general**: P_s = (T − D)V/W with D containing kW²/(qS); at fixed T/W the thrust term is the same, but the induced drag per unit weight, kW/(qS), grows with W (at fixed S), so the heavier aircraft has a lower ceiling. True only if the wing loading W/S is also kept the same.
(c) **True only if** the kinetic energy at the start and end is the same and the speed schedule is flat; in general the minimum-time path maximises P_s on each *energy* level, not ḣ on each altitude, and for supersonic aircraft the two are very different (Problem 4c shows the difference even for the trainer).
(d) **True only if** drag is negligible during the manoeuvre: an ideal zoom conserves h_E, but a real dive and pull-up lose energy to drag (Problem 5d: 24 m/s of energy height per second at n = 4), so the altitude regained is less than the altitude lost unless the engine makes up the difference.
(e) **First half true, second half true only at a given C_L**: in a glide T = 0 so P_s = −DV/W < 0; the glide *angle* depends only on L/D, which at a given C_L does not depend on weight, but the glide *speed* scales as √(W/S), so a heavier glider flies the same angle faster and arrives sooner.
(f) **True**: 10⁻⁸ m is ten orders of magnitude below the uncertainty of the thrust model (±10 % moves the ceiling by 600 m). Numerical accuracy is not the limit of this model; the physics is.

## Check-yourself answers

1. h_E = h + V²/(2g) [m]: mechanical energy per unit weight. P_s = (T − D)V/W [m/s]: excess power per unit weight. dh_E/dt = P_s.
2. Lift is perpendicular to the velocity, so L·v = 0. Drag is antiparallel to the velocity, so D·v = −DV < 0 at every instant.
3. Gravity is conservative (constant here, the gradient of mgz), so its work depends only on the endpoints. Drag depends on the direction of motion and is always opposed to it, so its work is minus its magnitude times the arc length: path-dependent.
4. h_E = 1000 + 100²/19.61 = 1509.9 m. Zoom to 50 m/s: 1509.9 − 127.5 = 1382.4 m.
5. At C_D0 = kC_L² (zero-lift drag equals induced drag). (L/D)max = 1/(2√(0.02 × 0.05)) = 15.81.
6. √2π × 100/9.807 = 45.3 s. Shorter: the density gradient adds gβ to ω².
7. n = 1/cos 60° = 2. P_s falls because the induced drag grows as n² while the thrust does not change.
8. The minimum-time schedule maximises P_s along each energy level. Without a transonic drag rise, P_s along an energy level peaks close to the best-climb speed, and the gain from exchanging speed and height is small. A supersonic fighter has a drag peak near Mach 1. Its P_s is highest on the far side of that peak, so the fastest path dives through the drag rise, trading height for speed, and then climbs at high speed.
