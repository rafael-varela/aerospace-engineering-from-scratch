# P00.12  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units and uncertainty on every numerical answer.

## 1. (Warm-up: three Fermi estimates)

Estimate, to the nearest power of ten and then to one significant figure, (a) the mass of air inside a hangar 90 m × 90 m × 30 m at sea level, (b) the total passenger-plus-baggage mass on a 180-seat airliner, (c) the kinetic energy of that airliner at 70 t and 230 m/s, and compare (c) with the chemical energy in 1 kg of jet fuel (take 43 MJ/kg). For each, state the reference values you used and your expected error factor. Answer form: three numbers with units and an error factor each, plus the comparison ratio.

## 2. (Derivation: minimum thrust and the speed where it occurs)

Starting from $T/W = C_{D0}/C_L + kC_L$, derive $C_L^*$, $(T/W)_{\min}$, and $(L/D)_{\max}$. Then, using $C_L = W/(qS)$, derive the speed $V^*$ at which $C_L^*$ occurs, and show how $V^*$ scales with weight, altitude (through density), and wing loading. Evaluate $V^*$ for the A320-class dataset at 11 km with $C_{D0} = 0.02$, $e = 0.8$, $AR = 9$. Answer form: derivation, three scaling statements, one speed in m/s.

## 3. (Dimensional reasoning: find the broken formula)

A colleague proposes three expressions for the range $R$ of a jet aircraft in cruise, where $V$ is speed (m/s), $c$ is the thrust-specific fuel consumption in kg of fuel per newton of thrust per second (kg/(N·s)), $L/D$ is dimensionless, $g$ is in m/s², and $W_i$, $W_f$ are initial and final weights in N:

(a) $R = \dfrac{V}{c}\,\dfrac{L}{D}\,\ln\dfrac{W_i}{W_f}$, (b) $R = \dfrac{V}{c\,g}\,\dfrac{L}{D}\,\ln\dfrac{W_i}{W_f}$, (c) $R = \dfrac{V}{c\,g}\,\dfrac{L}{D}\,\dfrac{W_i}{W_f}$.

Determine by dimensions alone which can be a range in metres. For the survivors, test a limiting case ($W_i = W_f$, no fuel burned) to eliminate further. Then evaluate the survivor for $V = 230$ m/s, $c = 1.7 \times 10^{-5}$ kg/(N·s), $L/D = 17$, $W_i/W_f = 1.3$, and say whether the result is a plausible airliner range. Answer form: a dimensional table, the surviving formula, a range in km, and a plausibility verdict.

## 4. (Numerical, uses your code: the plausible band)

Take a seed from the last two digits of your birth year, $s$. Set mass $= (60 + s/4) \times 10^3$ kg and wing area $= 110 + s/4$ m², at 11 000 m. Sweep speed from 100 to 400 m/s in 1 m/s steps with `check_airliner`. Report (a) the speed interval over which no flag trips, (b) the speed of minimum $T/W$ from the sweep and from equation (4) of the lesson, and their difference, (c) a plot of $T/W$ against speed with the analytic minimum marked, saved to a file. Answer form: an interval, two speeds, a plot.

## 5. (Unfamiliar scenario: a cargo drone)

A fixed-wing cargo drone has mass 25 kg, wing area 1.2 m², cruise speed 20 m/s at 500 m. Run `check_airliner` on it and list every flag with its reason. Then explain, for each flag, whether it indicates an error in the data or a limitation of the checker. Finally, propose a `PLAUSIBLE` table for small fixed-wing UAVs (ranges for mass, altitude, Mach, wing loading, $C_L$, $T/W$) with one line of justification per row, and show the drone passing it. Answer form: two tables and the flags.

## 6. (Judgment: thresholds and false alarms)

Run `check_airliner` on Concorde-class numbers (185 000 kg, 358 m², 600 m/s, 17 000 m) and on a regional turboprop (22 000 kg, 61 m², 150 m/s, 7 500 m). Both are real, consistent datasets. Which quantities are flagged for each, and is any flag a genuine data error? Propose and defend one of: (i) widening the ranges so both pass, (ii) keeping the ranges and adding a vehicle-class parameter, (iii) keeping the ranges and reporting flags as questions rather than errors. Give the cost of your choice in terms of missed errors or false alarms, with one example of each. Answer form: flag lists, a choice, and a two-paragraph defence.

## 7. (Judgment: the sanity check that was not run)

Read the root cause and contributing causes in the Mars Climate Orbiter mishap report (`@nasa-mco-mib`). The small-forces file reported impulse in lbf·s where N·s was expected, a factor of 4.45. Describe one sanity check, in the style of this lesson (a reference value, an expected magnitude, a comparison), that the navigation team could have run on each thruster firing and that would have flagged the file within the first weeks of cruise. Say what reference value it needs, where that value would have come from, and why a dimensional check alone would not have found the problem. Answer form: a described check in at most ten sentences.
