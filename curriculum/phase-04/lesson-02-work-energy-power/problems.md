# P04.02  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31). The aircraft is the lesson's trainer (`Aircraft()`) unless stated.

## 1. (Warm-up: the energy accounts)

The trainer flies at h = (2000 + 100s) m and V = (120 + s) m/s.

(a) Its energy height, kinetic energy and potential energy (from sea level). Which account holds more?
(b) It zooms (thrust and drag neglected) until the speed is 1.2 times the stall speed *at the altitude it reaches*. Find that altitude (one equation in one unknown: `brentq`) and the height gained.
(c) Its P_s at the starting point (`specific_excess_power`). If it climbs while accelerating at 2 m/s per 100 m of height, what is its rate of climb?

Answer form: three numbers with units and a sentence; two numbers; two numbers.

## 2. (Derivation: the special speeds and the phugoid)

(a) From d(½mV² + mgh)/dt = (T − D)V, derive dh_E/dt = P_s. Which force does no work, and why?
(b) For level flight with D = qSC_D0 + kW²/(qS), show that the drag is minimum where C_D0 = kC_L² and that D_min = W/(L/D)max with (L/D)max = 1/(2√(kC_D0)).
(c) Show that the power required DV is minimum where C_D0 = kC_L²/3, that V_mp = V_md/3^{1/4}, and that L/D there is (√3/2)(L/D)max.
(d) Derive ω² = 2g²/V₀² + gβ for the small-amplitude, drag-free phugoid at constant C_L, with β = −(1/ρ)dρ/dh. State every linearisation. Evaluate the period at 150 m/s at 3 km with β from `density_scale_gradient` and compare with `simulate_flight`.

Answer form: four derivations; two periods.

## 3. (Unfamiliar scenario: storing solar energy as altitude)

A solar-powered UAV of mass (40 + s) kg (S = 4 m², C_D0 = 0.015, k = 0.03, C_Lmax = 1.3) has 150 W of electrical power left over, after flying level, for 6 hours in the middle of the day. Instead of charging a battery, it climbs.

(a) Its rate of climb, and the altitude gained in 6 hours from 1000 m (neglect the change of P_s with altitude for this part; say why it is a fair first estimate).
(b) The energy stored, in MJ and Wh, and the mass of lithium-ion cells at 250 Wh/kg that would store the same.
(c) At dusk the motor stops and the UAV glides at its minimum-power speed. Using `replace(Aircraft(), mass=..., S=4.0, CD0=0.015, k=0.03, T0=0.0, CLmax=1.3)`, find V_mp and the sink rate at 1000 m and at the top altitude, check that V_mp is above the stall speed, and integrate the glide time from the top back to 1000 m.
(d) In two sentences: is "altitude as a battery" a good idea for this aircraft? What would you need to know to decide?

Answer form: two numbers; three numbers; five numbers and a time; two sentences.

## 4. (Numerical, uses your code: the minimum-time climb)

(a) With `rutowski_path` and `time_to_energy`, find the minimum time from (0 m, 150 m/s) to h_E = (9000 + 100s) m, and the speed and altitude at the end of the schedule.
(b) Compare with `climb_time_accelerate_then_constant` at 130, 150 and 160 m/s. Explain the ranking and the size of the differences.
(c) A "steady best-climb" estimate integrates dh/P_s,max(h) from sea level to the final altitude of (a), ignoring the kinetic energy. Compute it and explain why it is *shorter* than the energy-method time, and why it is wrong.
(d) Repeat (a) with 10 % more mass. By how much does the time grow, and why is the growth larger than 10 %?

Answer form: a time and two numbers; three times with two sentences; a time and two sentences; a time and a sentence.

## 5. (Turns: what a hard turn costs)

At 5 km and 150 m/s:

(a) P_s at load factors 1, 1.5, 2 and 3, with the bank angles. Explain the trend from the drag formula.
(b) The maximum sustained load factor (P_s = 0: `brentq` on n), the bank angle, the turn rate g√(n² − 1)/V and the turn radius V²/(g√(n² − 1)).
(c) Find the speed at which the sustained load factor at 5 km is largest, between 110 m/s and the Mach cap, and its value.
(d) The aerodynamic limit is n = qSC_Lmax/W. Compute it at 150 m/s. A pilot pulls to n = 4 for 5 s: how much energy height is lost (P_s at n = 4 times 5 s), and how much speed if the altitude is held?

Answer form: a four-row table and a sentence; four numbers; two numbers; three numbers.

## 6. (Judgment: trusting the energy method)

For each statement, say whether it is true, false, or "true only if…", in one or two sentences each.

(a) "The zoom altitude does not depend on the aircraft's mass."
(b) "A heavier aircraft with the same T/W has the same ceiling."
(c) "The minimum-time path is the path of maximum rate of climb."
(d) "An aircraft can always regain the altitude lost in a dive by zooming back up."
(e) "Drag does negative work, so P_s is always negative in a glide, and the glide ratio therefore does not depend on the weight."
(f) "If the simulated energy height drifts by 10⁻⁸ m over 300 s in a drag-free flight, the integration is accurate enough for performance work."

Answer form: six labelled items.
