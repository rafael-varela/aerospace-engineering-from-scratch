# P06.02  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe script that calls it, given to the figures shown. Aircraft numbers are illustrative.

## 1. A pull-up and a turn at the same g (s = 15: V = 165 m/s)

**Assumptions.** Particle, constant speed, g = 9.80665 m/s².

(a) R = V²/(g(n − 1)) = 165²/(9.80665 × 5) = **555.2 m**. Over the top at the same R and V: n = V²/(gR) − 1 = 5 − 1 = **4.0** (bottom minus top is always 2 g).

(b) φ = arccos(1/6) = **80.41°**; R = V²/(g√35) = **469.3 m**; ω = V/R = 0.3516 rad/s = **20.15 °/s**.

(c) Seat: n m g = 6 × 85 × 9.80665 = **5001 N**. Neck: 6 × 6 × 9.80665 = **353 N**: the neck carries the head and helmet at six times their weight, which is why helmet-mounted displays are weighed in grams.

(d) n_top = V²/(gR) − 1 ≥ 0 ⇒ V ≥ √(gR) = **99.0 m/s** for R = 1000 m.

**Independent check.** At n = 6 the level turn is tighter than the pull-up by √35/5 = 1.183: 555.2/469.3 = 1.183 ✓.

**Common mistakes.** Using n instead of n − 1 for the pull-up radius; forgetting that the turn's radius uses √(n² − 1).

## 2. Path equations and the dive recovery

(a) Along e_t: T − D − W sin γ = m V̇ (thrust assumed along the velocity: a thrust line inclined at ε would add T sin ε to the normal equation and multiply T by cos ε in the tangential one). Along e_n: L − W cos γ = m V γ̇, since the velocity direction turns at γ̇ and a_n = V γ̇.

(b) An accelerometer measures the specific force f = a − g_vec (the push of its mounting, per unit mass). Along e_n: f·e_n = V γ̇ + g cos γ = (L − W cos γ)/m + g cos γ = L/m, so f·e_n/g = L/W = n. Level flight: **1**. Zero-g parabola: **0**. Level turn at 60° bank: **2** (1/cos 60°; the accelerometer is along the lift, tilted with the aircraft).

(c) Constant V and n: dh/dγ = ḣ/γ̇ = V sin γ / (g(n − cos γ)/V) = (V²/g) sin γ/(n − cos γ). Since d(n − cos γ)/dγ = sin γ, ∫_{γ₀}^{0} = (V²/g)[ln(n − cos γ)]_{γ₀}^{0} = (V²/g) ln((n − 1)/(n − cos γ₀)): a negative change of height, so the loss is **Δh = (V²/g) ln((n − cos γ₀)/(n − 1))**. Time: dt = V dγ/(g(n − cos γ)); with u = tan(γ/2), cos γ = (1 − u²)/(1 + u²) and dγ = 2du/(1 + u²), so ∫ dγ/(n − cos γ) = ∫ 2du/((n − 1) + (n + 1)u²) = (2/√(n² − 1)) atan(√((n + 1)/(n − 1)) u), which is (4). Small angle: n − cos γ₀ ≈ (n − 1)(1 + γ₀²/(2(n − 1))), so Δh ≈ (V²/g) γ₀²/(2(n − 1)) = **R γ₀²/2** with **R = V²/(g(n − 1))**, the drop of a circular arc of that radius through angle γ₀.

**Common mistakes.** Integrating dh/dt with γ treated as constant; losing the sign of sin γ for a dive.

## 3. A zero-g research flight

(a) R = 160²/(9.80665 × 0.8) = **3263 m**.

(b) With n = 0, L = 0; with T = D the only force is the weight, so ẍ = 0, ḧ = −g: a projectile parabola. Entering at 160 m/s and 45°, the horizontal speed 113.14 m/s is constant and the vertical speed goes from +113.14 to −113.14 m/s (γ = −45° when the two are equal and opposite). Duration **2 V sin 45°/g = 23.07 s**; height gained to the apex **(V sin 45°)²/(2g) = 652.6 m**; horizontal distance **V²/g = 2610.5 m**; exit speed **160 m/s** (the same height as entry, so the same speed by energy).

(c) `solve_ivp` on `flight_path_rhs` with n = 0, `speed_hold=False`, `excess_thrust=0` and an event on γ = −45°: 23.07354 s, x = 2610.47 m, exit speed 160.0000 m/s, agreeing with (b) to better than 10⁻⁸ relative.

(d) At the apex V = 113.14 m/s and the whole of g is normal to the path: ρ = V²/g = **1305 m**. n = 0 means the *lift* is zero, not the acceleration: the path still curves under gravity.

**Common mistakes.** Thinking n = 0 means a straight line; using the entry speed at the apex.

## 4. Your own dive recovery (s = 15: γ₀ = −45°)

(a) Δh = **303.76 m**, t = **3.283 s**.

(b) 2 g/s: **539.4 m** (extra 235.6 m; estimate (5) 265.2 m). 6 g/s: **386.2 m** (extra 82.5 m; estimate 88.4 m). The estimate is 7–13 % high, as in the lesson, because the path starts to turn during the ramp.

(c) Thrust equal to drag at 6 g/s: **409.8 m**, 23.6 m more: the aircraft accelerates down the dive while the g builds, and a faster aircraft needs a bigger radius.

(d) The load factor recovered from the positions agrees with the 6 g/s ramp to **3.3 × 10⁻⁵** away from the ends and the kink at 1 s; the largest error anywhere, at the ends where `np.gradient` is one-sided, is 5.6 × 10⁻³.

**Common mistakes.** Using degrees in `np.cos`; starting the onset ramp at 0 g instead of 1 g.

## 5. Setting a recovery trigger

(a) At 6 g/s from 60°: 501.9 m at 225 m/s and 725.9 m at 275 m/s against 608.8 m: **−18 % and +19 %** for ±10 % speed. At 6 g instead of 7 g: 688.8 m, **+13 %**.

(b) Larger than predicted: a lift-limited pull (n not available), a slower real onset or pitch lag, thrust left high (the aircraft speeds up), a vertical gust or downdraft. Smaller: speed falling below the start value (speed brakes, drag), a higher n than planned, a dive that is shallower than estimated.

(c) n W = ½ ρ V² S C_L,max ⇒ V = √(2 n W/(ρ S C_L,max)) = **215 m/s** for 7 g (244 m/s for 9 g). Below 215 m/s the aircraft cannot pull 7 g with this wing; the recovery is flown at a lower n and needs more height than (3) with n = 7, so a trigger computed with a fixed 7 g would fire too late at low speed. The trigger must use the n actually available at the current speed.

(d) One defensible choice: predicted loss at the available n and a 4 g/s onset, plus 20 % for a 10 % speed uncertainty, plus a fixed floor (say 50 m) for terrain-database error. Speed and onset dominate the uncertainty, so the margin should scale with V² and with 1/ṅ rather than be a fixed number.

## Check-yourself answers

1. m V̇ = T − D − W sin γ; m V γ̇ = L − W cos γ.
2. n = L/W. In a free-fall parabola it reads 0.
3. R = V²/(g(n − 1)) = 200²/(9.80665 × 4) = 1019.7 m.
4. At constant speed and radius n = cos γ + V²/(gR); cos γ is +1 at the bottom and −1 at the top.
5. The level turn: its radius is V²/(g√(n² − 1)) against V²/(g(n − 1)) for the pull-up; in the pull-up one g of lift holds up the weight, in the turn the weight is balanced by the vertical component and more is left for turning.
6. As V² (at fixed n and dive angle), from (3).
7. During the ramp the aircraft is still diving with less than the full load factor; roughly it acts like a delay of half the ramp time, costing about V sin|γ₀| (n − 1)/(2ṅ).
8. V²/2 + gh; the lift is perpendicular to the velocity, so it does no work.
