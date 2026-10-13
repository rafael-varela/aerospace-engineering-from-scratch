# P00.12  Worked solutions

Each solution states: assumptions, units, the derivation, the result with uncertainty where applicable, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes.

## 1.

**Assumptions.** Sea-level air 1.2 kg/m³ (±10 %); 100 kg per passenger with baggage (±20 %); jet fuel 43 MJ/kg as given.
**Derivation.** (a) $90 \times 90 \times 30 = 2.43 \times 10^5$ m³; × 1.2 = $2.9 \times 10^5$ kg. (b) $180 \times 100 = 1.8 \times 10^4$ kg. (c) $\tfrac12 \times 7 \times 10^4 \times 230^2 = 0.5 \times 7 \times 10^4 \times 5.29 \times 10^4 = 1.85 \times 10^9$ J; divided by $4.3 \times 10^7$ J/kg gives 43 kg of fuel.
**Result.** (a) $\sim 10^5$ kg; $3 \times 10^5$ kg (factor 1.3). (b) $\sim 10^4$ kg; $2 \times 10^4$ kg (factor 1.3). (c) $\sim 10^9$ J; $1.9 \times 10^9$ J (factor 1.1, since mass and speed are given). The airliner's cruise kinetic energy equals the chemical energy of only about 40 kg of fuel, roughly 0.2 % of a typical fuel load; almost all the fuel goes into overcoming drag, not into getting up to speed.
**Check.** (a) The air outweighs the aircraft's empty mass (about 40 t) by a factor of 7: surprising but right; dimensions m³ × kg/m³ = kg. (c) Units: kg·m²/s² = J.
**Common mistakes.** Forgetting to square the speed; using 1.2 g/L and dropping a factor of 1000; quoting kinetic energy "per passenger" when total was asked.

## 2.

**Assumptions.** Parabolic polar with constant $C_{D0}$, $e$, $AR$; level flight; ISA.
**Derivation.** $\frac{d}{dC_L}(C_{D0}/C_L + kC_L) = -C_{D0}/C_L^2 + k = 0 \Rightarrow C_L^* = \sqrt{C_{D0}/k}$. Substituting: $(T/W)_{\min} = C_{D0}/\sqrt{C_{D0}/k} + k\sqrt{C_{D0}/k} = \sqrt{C_{D0}k} + \sqrt{C_{D0}k} = 2\sqrt{C_{D0}k}$, and $(L/D)_{\max} = 1/(2\sqrt{C_{D0}k})$. Speed: $C_L^* = 2W/(\rho V^{*2} S)$ gives
$$V^* = \sqrt{\frac{2}{\rho}\,\frac{W}{S}\,\frac{1}{C_L^*}} = \sqrt{\frac{2\,(W/S)}{\rho\sqrt{C_{D0}/k}}}.$$
Scaling: $V^* \propto \sqrt{W/S}$ (a heavier aircraft on the same wing flies its best L/D faster); $V^* \propto \rho^{-1/2}$ (higher altitude, faster); $V^*$ is independent of weight at fixed wing loading and independent of $C_L$'s value except through the polar constants. With $k = 1/(\pi \times 0.8 \times 9) = 0.04421$: $C_L^* = \sqrt{0.02/0.04421} = 0.672$; $(T/W)_{\min} = 2\sqrt{0.02 \times 0.04421} = 0.0595$; $(L/D)_{\max} = 16.8$. $W/S = 5600$ N/m², $\rho = 0.364$ kg/m³: $V^* = \sqrt{2 \times 5600/(0.364 \times 0.672)} = \sqrt{45\,790} = 214$ m/s.
**Result.** $V^* = 214$ m/s (Mach 0.73 at 11 km), below the 230 m/s cruise speed, as expected: airliners cruise somewhat faster than minimum-drag speed because the $T/W$ curve is flat there and time matters.
**Check.** The test `test_thrust_to_weight_minimum_at_analytic_optimum` finds the same $C_L^*$ by sweeping speed. Limiting case: $C_{D0} \to 0$ makes $C_L^* \to 0$, meaning fly as fast as possible when there is no parasite drag, which is the right direction.
**Common mistakes.** Dropping the factor 2 in $(T/W)_{\min}$; using sea-level density; confusing $W/S$ in N/m² with kg/m².

## 3.

**Assumptions.** $c$ in kg/(N·s) as stated; $W$ in newtons.
**Derivation.** $[V/c] = (\text{m/s}) / (\text{kg}/(\text{N·s})) = \text{m·N}/\text{kg} = \text{m}\cdot\text{kg·m/s}^2/\text{kg} = \text{m}^2/\text{s}^2$: not a length, so (a) fails. $[V/(cg)] = (\text{m}^2/\text{s}^2)/(\text{m/s}^2) = \text{m}$: (b) and (c) pass dimensionally. Limiting case $W_i = W_f$: (b) gives $\ln 1 = 0$, range zero, correct; (c) gives $R = V L/(c g D)$, a nonzero range with no fuel burned, which is impossible. Only (b) survives (it is the Breguet range equation; Phase 16). Evaluate: $V/(cg) = 230/(1.7 \times 10^{-5} \times 9.807) = 1.380 \times 10^6$ m; × 17 = $2.346 \times 10^7$ m; × $\ln 1.3 = 0.2624$ gives $6.16 \times 10^6$ m.
**Result.** $R \approx 6200$ km (to about 10 %, given the one-significant-figure inputs). Plausible: a narrow-body airliner burning 23 % of its weight in fuel flies roughly 5000–7000 km.
**Check.** Order of magnitude: $10^6$ m × 10 × 0.3 ≈ $3 \times 10^6$ to $10^7$ m. Another way to see that $V/(cg)$ is a length: $cg$ has units (kg/(N·s))(m/s²) = (kg·m/s²)/(N·s) = 1/s, so $1/(cg)$ is a time and $V/(cg)$ is speed × time.
**Common mistakes.** Treating $c$ as dimensionless; accepting (c) because its units are right (dimensions are necessary, not sufficient); forgetting the logarithm.

## 4.

**Assumptions.** Seed-dependent; shown for $s = 90$: mass 82 500 kg, wing area 132.5 m², 11 km, default polar.
**Derivation.** $W/S = 6106$ N/m². Sweep: $C_L$ falls below 0.8 above $V = \sqrt{2 \times 6106/(0.364 \times 0.8)} = 205$ m/s and below 0.3 above 334 m/s; Mach exceeds 0.92 above 271 m/s and falls below 0.6 under 177 m/s; $T/W$ stays within 0.04–0.12 for $C_L$ between about 0.19 and 2.4, so it is not binding here. Passing interval: about 205–271 m/s. Minimum $T/W$ from the sweep at the speed where $C_L = 0.672$: $V^* = \sqrt{2 \times 6106/(0.364 \times 0.672)} = 223$ m/s; the sweep's argmin agrees to within the 1 m/s step.
**Result.** (a) ≈ 205–271 m/s for this seed; (b) 223 m/s both ways, difference < 1 m/s; (c) plot with the minimum marked at $T/W = 0.0595$.
**Check.** The lower bound is set by $C_L = 0.8$ and the upper by Mach 0.92, as the hint said. The interval must contain $V^*$ for the dataset to be a sensible airliner; if it does not, the wing loading is unusual.
**Common mistakes.** Sweeping at sea-level density; reading the minimum from a coarse grid and reporting a spurious 5 m/s difference; forgetting `matplotlib.use("Agg")`.

## 5.

**Assumptions.** ISA at 500 m: $\rho = 1.167$ kg/m³, $a = 338$ m/s.
**Derivation.** $W = 245$ N; $W/S = 204$ N/m²; $q = 0.5 \times 1.167 \times 400 = 233$ Pa; $C_L = 0.88$; Mach 0.059; with the airliner polar, $C_D = 0.02 + 0.0442 \times 0.77 = 0.0542$, $T/W = 0.062$. Flags: mass (too light), altitude (too low), Mach (too slow), wing loading (too low), $C_L$ (above 0.8). $T/W$ passes only because the airliner polar flatters the drone.
**Result.** Five flags; every one is a true statement about the drone and none is a data error; the checker's ranges describe a different vehicle class. Proposed UAV table: mass 1–50 kg (hand-launch to small-runway class); altitude 0–3000 m (regulatory and battery limits); Mach 0.03–0.15 (10–50 m/s); wing loading 50–400 N/m² (foam gliders to heavy cargo drones); $C_L$ 0.4–1.0 (slow flight, high-lift sections); $T/W$ 0.08–0.30 ($L/D$ of 3–12 for small, draggy airframes). The drone passes every row.
**Check.** Dimensionless numbers transfer: a drone and an airliner can be compared on $C_L$, and the drone's 0.88 is high but attainable with a cambered section and no flaps; its Mach says compressibility is irrelevant.
**Common mistakes.** Treating the flags as errors in the drone's data; proposing a UAV table with airliner $T/W$ values.

## 6.

**Assumptions.** ISA; airliner polar for both (a limitation itself).
**Derivation.** Concorde-class: $\rho(17\,\text{km}) = 0.1413$ kg/m³, $a = 295$ m/s, Mach 2.03 (flagged: beyond the transonic drag rise); $W/S = 5068$ N/m² (passes); $q = 25\,600$ Pa, $C_L = 0.198$ (flagged: very low); altitude 17 km (flagged: above airliner ceiling); $T/W = 0.02/0.198 + 0.0442 \times 0.198 = 0.110$ (passes, barely, and for the wrong reason: the subsonic polar does not apply). Turboprop: $\rho(7.5\,\text{km}) = 0.557$, $a = 311$, Mach 0.48 (flagged: slower than jet cruise); $W/S = 3537$ (passes); $q = 6266$, $C_L = 0.56$ (passes); altitude 7.5 km (passes, inside 6–14); mass 22 t (passes). One genuine flag, Mach, and it is a correct description of a turboprop.
**Result.** No flag on either aircraft is a data error. Defence of option (ii), vehicle class as a parameter: the ranges encode knowledge about a class and are only meaningful when the class is known; asking the user to name the class costs one keyword and keeps the airliner ranges tight enough to catch the km/h mistake (Mach 1.5), which option (i) would let through (a widened Mach range to 2.1 passes it). Cost of (ii): a user who picks the wrong class, or a vehicle between classes (a supersonic business jet), gets wrong verdicts; the mitigation is option (iii)'s habit of reporting flags as questions. Example error still caught under (ii): mass in pounds, through $C_L$ and wing loading. Example missed: any 15 % input error.
**Check.** Concorde's low $C_L$ is physically right: supersonic cruise at high $q$ needs little lift coefficient, which is why the ogival delta wing is small relative to its mass.
**Common mistakes.** Concluding the datasets are wrong; choosing (i) without naming the mistake it lets through.

## 7.

**Assumptions.** From the mishap report: the small-forces file carried angular-momentum-desaturation thruster impulse in lbf·s, read as N·s, a factor of 4.45; the effect accumulated over many firings during cruise.
**Result.** The check: for each desaturation event, compute the expected velocity change $\Delta v = F\,\Delta t / m$ from the thruster's rated thrust (newtons, from the propulsion specification), the commanded firing duration (seconds, from the command log), and the current spacecraft mass (kilograms, from the mass budget). Compare it with the $\Delta v$ the navigation team's orbit-determination filter needed to fit the tracking data. The expected value is a reference magnitude; the fitted value is the measurement. A ratio of 4.45 between them, the same ratio on every event, is far outside the few-percent uncertainty of either and is the fingerprint of a constant conversion factor rather than a physics error. A dimensional check could not have caught it because lbf·s and N·s have the same dimensions (force × time); only a *magnitude* comparison against an independently derived expected value distinguishes unit systems. The contributing causes in the report say in effect that the discrepancy was seen but not pursued to this comparison.
**Check.** 1 lbf = 4.448 N, so the ratio predicted by the hypothesis "file in lbf·s" is 4.45, matching the observed under-prediction of the trajectory perturbation.
**Common mistakes.** Proposing a dimensional check as the fix; proposing "use SI" without a comparison that would detect a violation.

## Check-yourself answers

1. Wing loading $\sim 10^3$ to $10^4$ N/m² (nearest power of ten for 5600 is $10^4$); density at 11 km $\sim 10^{-1}$ kg/m³ (0.36 rounds to $10^0$ by the nearest-power rule, which is why one-significant-figure values are better than orders of magnitude for numbers near $\sqrt{10}$); speed of sound $\sim 10^2$ m/s.
2. (kg/m³)(m²/s²) = kg/(m·s²) = N/m² = Pa, a pressure.
3. $C_L \propto 1/V^2$: $1/1.1^2 = 0.826$, so about 17 % lower.
4. The ISA relations written in the file are only valid to 20 km; above that the temperature profile changes and the formulas would return a wrong number silently. A meaningless input is refused; an implausible one is flagged.
5. $C_L^* = \sqrt{\pi e\,AR\,C_{D0}}$; it does not depend on weight or altitude. Those set the speed at which $C_L^*$ occurs.
6. Mass in pounds entered as kilograms; speed in knots entered as m/s (too slow); or wing area in a smaller unit. Physically, a clean wing's maximum lift coefficient is about 1.2–1.5, so cruising at 1.3 would leave no margin to a stall and would require far more thrust than any airliner carries.
7. lbf·s and N·s are both force × time; dimensions match. A comparison of each firing's effect against its expected magnitude from thrust × duration / mass would have shown a constant factor of 4.45.
8. An error that lands inside the plausible band (a 15 % mass error), and a vehicle of the wrong class whose correct numbers happen to fall inside the airliner ranges (a large business jet).
