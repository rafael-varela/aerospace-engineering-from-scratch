# P00.12  Engineering Judgment & Sanity Checks

> Before you trust a number, ask what size it ought to be. Before you trust a formula, ask what units it produces.

| | |
|---|---|
| **Learning outcomes** | 1. Produce an order-of-magnitude estimate for an aerospace quantity from a few remembered reference values, and state its expected error. 2. Check a formula by dimensions and by limiting cases before using it. 3. Derive the lift coefficient, wing loading, Mach number, and thrust-to-weight ratio an airliner needs in cruise from mass, wing area, speed, and altitude, using the ISA relations. 4. Implement and test a sanity checker that flags implausible values with the reason, and explain which errors it cannot catch. 5. Decide, for a given set of aircraft numbers, whether they are mutually consistent and which one is most likely wrong. |
| **Prerequisites** | `P00.03` (the ISA table), `P00.05` (units and the Mars Climate Orbiter); skill check: can rearrange $L = \tfrac12 \rho V^2 S C_L$ for any one variable. |
| **Engineering question** | Given four numbers for an airliner (mass, wing area, cruise speed, altitude), can you tell within a minute whether they are mutually consistent, and if not, which one is wrong? |
| **Assumptions and limits** | Steady, level, unaccelerated cruise so that lift equals weight and thrust equals drag. ISA atmosphere between −500 m and 20 km, no wind, no temperature deviation. A parabolic drag polar $C_D = C_{D0} + C_L^2/(\pi e\,AR)$ with assumed $C_{D0}$, Oswald factor $e$, and aspect ratio. Plausibility ranges are for subsonic jet airliners; a Concorde, a turboprop, or a drone will be flagged even when its numbers are right. The checker catches gross errors (a factor of 2 or more), not 10 % errors. |
| **Analytical work** | Deriving the ISA pressure law from hydrostatics and the ideal gas; deriving the minimum thrust-to-weight ratio and the $C_L$ at which it occurs; a dimensional check of the lift equation; three Fermi estimates. |
| **Computational deliverable** | `starter/main.py`: `isa`, `lift_coefficient`, `wing_loading`, `mach_number`, `thrust_to_weight_required`, `order_of_magnitude`, and `check_airliner`, which returns one `Flag` per quantity with its plausible range and a reason when it is outside. Tests prove the ISA constants and continuity at the tropopause, the $1/V^2$ scaling of $C_L$, the analytic optimum of $T/W$, that A320-, 777-, and 747-class inputs pass, and that mass-in-pounds, speed-in-km/h, and a tenth-size wing are flagged with reasons. |
| **Practical activity** | **computational experiment**: inject five classic unit mistakes into a correct airliner dataset and record which the checker catches; sweep speed to find the plausible band and compare its minimum-thrust point with the analytic optimum. |
| **Equipment** | None beyond a computer. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | The pytest output, the injected-error table (mistake, flagged quantities, reason), a plot of $T/W$ against speed with the analytic minimum marked, and a paragraph on the one mistake the checker missed and why. |
| **Assessment** | `assessments/phase-00/grading-rubric.md`; competencies 1 (derivations), 2 (Fermi estimates), 3 (code), 4 (judgment about thresholds). Pass: tests green, at least four of five injected mistakes caught and the fifth explained, and problem 6 argued with numbers. |
| **References** | `@nasa-mco-mib` "Mars Climate Orbiter Mishap Investigation Board Phase I Report" (required-free; read the root cause and contributing causes as the case study for this lesson). `@ocw-18-335j` Week 2 "Lecture 5: Condition Numbers" (supplementary-free; graduate level, read only for the idea that a problem can amplify input errors). |
| **Extensions** | ★ Add a `PLAUSIBLE` table for a second vehicle class (regional turboprop, business jet, or small UAV) selected by a keyword argument. ★ Add a Reynolds-number check using the mean chord $S/b$ from a span input. |
| **Time** | unmeasured |

## 1. The engineering problem

The Mars Climate Orbiter was lost in September 1999 because a ground software file produced small-forces data in pound-force seconds while the navigation software expected newton seconds, a factor of 4.45. The investigation board's report lists the root cause and then a list of contributing causes. Read them. Several say the same thing in different words: the trajectory was drifting from the expected one, people noticed, and the discrepancy was never run to ground as *a number that could not be right*. Nobody computed how large each thruster firing's effect on the trajectory ought to be and compared it with how large the navigation team was seeing. The number was off by more than a factor of four for months.

That skill, knowing roughly what size a number ought to be and noticing when it is not, is what engineers mean by judgment. It is not mysterious and it is not only experience. It is four habits that can be taught and practised: estimating to the nearest power of ten from a few remembered values, checking that every formula produces the right units, checking what a formula does in extreme cases, and keeping a mental table of what plausible values look like for the quantities in your field. This lesson teaches all four on one problem small enough to do in a minute by hand.

The problem: a colleague hands you four numbers for an airliner, mass 70 000 kg, wing area 122.6 m², cruise speed 230 m/s, cruise altitude 11 000 m, and asks whether a design study based on them is sound. You will compute what those numbers imply (how hard the wing is working, how fast relative to sound, how much thrust the engines must deliver), compare each implication with the range real airliners occupy, and flag anything outside it with a reason. Then you will do it in code, and find out what a checker like this can and cannot catch. Success means your checker passes three real airliner-class datasets, flags five classic unit mistakes, and you can say in one sentence which kind of mistake no sanity check will ever find.

## 2. Understand and derive

### 2.1 Order-of-magnitude estimation

An order of magnitude is a power of ten. "To order of magnitude" means you aim for the nearest power of ten and accept a factor-of-three error. The method, named after Enrico Fermi, is to break a quantity you do not know into factors you roughly do, estimate each, and multiply. Errors in the factors tend to partially cancel.

**How much does the air in a hangar weigh?** A hangar for one airliner is perhaps 80 m × 80 m × 25 m $= 1.6 \times 10^5$ m³. Air at sea level is about 1.2 kg/m³. So about $2 \times 10^5$ kg: two hundred tonnes, three times the empty mass of the aircraft inside it. Many people guess "a few tonnes". The habit of actually multiplying beats the guess every time.

**How much of an airliner's mass is passengers?** 180 passengers at roughly 80 kg plus 20 kg of baggage each: $180 \times 100 = 1.8 \times 10^4$ kg, about 18 t. Against a 70 t take-off mass that is a quarter. Fuel for a medium flight is about the same again. Those two facts let you check any mass figure you are handed to within a factor of two.

Three reference values make most aerospace estimates possible: sea-level air density 1.2 kg/m³, the speed of sound at sea level 340 m/s (295 m/s at cruise altitude, where the air is cold), and $g \approx 10$ m/s². Add the ISA numbers from `P00.03` at 11 km, pressure about a quarter and density about 0.30 of sea level, and you can estimate any lift, drag, or Mach number to 30 %.

### 2.2 Dimensional reasoning

Units obey algebra. Every term added to another must have the same dimensions, and the two sides of an equation must match. Checking this takes thirty seconds and catches a large fraction of all algebra mistakes. Take the lift equation and write each factor's SI units:

$$ L = \tfrac12\,\rho\,V^2\,S\,C_L : \quad \frac{\text{kg}}{\text{m}^3}\cdot\frac{\text{m}^2}{\text{s}^2}\cdot \text{m}^2 \cdot 1 = \frac{\text{kg}\cdot\text{m}}{\text{s}^2} = \text{N}. \qquad (1) $$

Lift comes out in newtons, as a force must. A proposed formula $L = \tfrac12 \rho V S C_L$ gives kg/s, which is a mass flow rate, not a force; the formula is wrong before you compute anything. In the Mars Climate Orbiter case the quantity was an impulse, force × time; its SI unit is N·s and the file carried lbf·s. Dimensions were consistent (both are force × time), so a dimensional check alone would have passed. That is the limit of dimensional reasoning: it checks the *kind* of quantity, not the *unit system*, and a unit-system error needs a magnitude check instead (§2.4).

Dimensional reasoning also tells you which combinations of variables are meaningful on their own. $C_L$, Mach number $M = V/a$, and Reynolds number are all dimensionless. A dimensionless number has the same value in every unit system, which is exactly why engineers compare aircraft by $C_L$ and $M$ rather than by lift in newtons and speed in m/s: a 747 and a model glider can be compared on the same axes.

### 2.3 Limiting cases

A formula you trust must behave sensibly when one input goes to zero or to infinity. For level flight, lift equals weight, $W = mg$, so

$$ C_L = \frac{W}{\tfrac12 \rho V^2 S} = \frac{W/S}{q}, \qquad q = \tfrac12 \rho V^2. \qquad (2) $$

$W/S$ is the **wing loading** (N/m²) and $q$ is the **dynamic pressure** (Pa). Limiting cases: as $V \to 0$, $C_L \to \infty$: a wing cannot hold an aircraft up at zero speed, correct. As $S \to \infty$, $C_L \to 0$: an infinitely large wing barely needs to work. As altitude rises, $\rho$ falls, so at constant speed $C_L$ rises: this is why cruise altitude is limited for a given speed. Each limit matches physical sense, so equation (2) has passed a check that costs nothing.

Thrust for level flight equals drag, and with the parabolic polar $C_D = C_{D0} + k C_L^2$, $k = 1/(\pi e\,AR)$,

$$ \frac{T}{W} = \frac{D}{L} = \frac{C_D}{C_L} = \frac{C_{D0}}{C_L} + k\,C_L. \qquad (3) $$

Differentiate with respect to $C_L$ and set to zero: $-C_{D0}/C_L^2 + k = 0$, so the thrust-to-weight ratio is smallest at

$$ C_L^* = \sqrt{C_{D0}/k} = \sqrt{\pi e\,AR\,C_{D0}}, \qquad \left(\frac{T}{W}\right)_{\min} = 2\sqrt{C_{D0}\,k}, \qquad \left(\frac{L}{D}\right)_{\max} = \frac{1}{2\sqrt{C_{D0}k}}. \qquad (4) $$

Check the limits of equation (3): $C_L \to 0$ (flying very fast) makes $T/W \to \infty$ because of the $C_{D0}/C_L$ term, the parasite drag of a fast aircraft; $C_L \to \infty$ (flying very slowly) makes $T/W \to \infty$ through $kC_L$, the induced drag of a slow one. There must be a minimum between, and equation (4) says where. For $C_{D0} = 0.02$, $e = 0.8$, $AR = 9$: $k = 0.0442$, $C_L^* = 0.67$, $(T/W)_{\min} = 0.060$, $(L/D)_{\max} = 16.8$. These are the right size for an airliner, which is itself a check on the assumed polar. Notice that $C_L^*$ does not depend on weight or altitude at all; weight and altitude set the *speed* at which that $C_L$ occurs.

![Curve of thrust-to-weight ratio T/W against lift coefficient C_L for C_D0 = 0.02 and k = 0.0442. T/W rises steeply at low C_L (fast flight, parasite drag C_D0/C_L dominates) and rises again at high C_L (slow flight, induced drag k C_L dominates), with a flat minimum (T/W)min = 2√(C_D0 k) = 0.060 at C_L* = 0.67.](figures/thrust-to-weight-vs-cl.svg)
*Figure 1. Equation (3) for C_D0 = 0.02, k = 0.0442. Parasite drag dominates at low C_L (fast), induced drag at high C_L (slow), and the curve is flat near its minimum of 0.060 at C_L\* = 0.67.*

The curve is flat near its minimum, which is why an airliner can cruise 10 % faster than the minimum-drag speed (lower $C_L$) for a few percent more thrust, and why a $C_L$ far from $C_L^*$ in either direction is a warning sign.

### 2.4 "Does this number make sense?" A reference table

The fourth habit is a table in your head. Here is a starting one for subsonic jet airliners in cruise; build your own as you go through the course.

| Quantity | Plausible airliner range | Why |
|---|---|---|
| Mass | 20 t – 600 t | Regional jet to the largest airliner ever certified. |
| Cruise altitude | 6 – 14 km | Engines and wings are optimised for the thin, cold air near the tropopause. |
| Cruise Mach | 0.60 – 0.92 | Below the transonic drag rise; above about 0.85 drag climbs steeply. |
| Wing loading | 3 – 8 kN/m² | Set by the compromise between cruise efficiency and approach speed. |
| Cruise $C_L$ | 0.3 – 0.8 | Near $C_L^*$ from equation (4); a clean wing's maximum is about 1.2–1.5. |
| Cruise $T/W$ | 0.04 – 0.12 | $1/(L/D)$ with $L/D$ between 8 and 25; airliners sit near 15–20. |
| Cruise $L/D$ | 15 – 20 | Sailplanes reach 40–60; fighters are below 10. |
| ISA density at 11 km | 0.36 kg/m³ | About 0.30 of sea level. |
| Speed of sound at 11 km | 295 m/s | Cold air; 340 m/s at sea level. |

A number outside its range is not necessarily wrong. It is a number that needs a reason. The checker you build gives the reason it suspects, which is the start of the conversation, not the end.

### 2.5 The ISA relations

The checker needs density and temperature at altitude. In `P00.03` you built the ISA table; here you write the relations again in the checker so the file stands alone, which is also a chance to see where they come from.

Two laws and one assumption. The hydrostatic equation says pressure falls with height at the rate the air above weighs: $dp/dh = -\rho g$. The ideal gas law says $p = \rho R T$ with $R = 287.05$ J/(kg·K). The ISA assumes the troposphere's temperature falls linearly, $T = T_0 + Lh$ with $T_0 = 288.15$ K and $L = -0.0065$ K/m, up to 11 km, and is constant at 216.65 K from 11 to 20 km. Substituting $\rho = p/(RT)$ into hydrostatics,

$$ \frac{dp}{p} = -\frac{g}{R}\,\frac{dh}{T_0 + Lh} \quad\Rightarrow\quad \frac{p}{p_0} = \left(\frac{T}{T_0}\right)^{-g/(LR)} \quad (0 \le h \le 11\ \text{km}), \qquad (5) $$

with exponent $-g/(LR) = 9.80665/(0.0065 \times 287.05) = 5.2559$. In the isothermal layer $T$ is constant and the same integral gives an exponential,

$$ p = p_{11}\,\exp\!\left[-\frac{g\,(h - 11000)}{R\,T_{11}}\right] \quad (11 \le h \le 20\ \text{km}). \qquad (6) $$

Checks: at $h = 0$ equation (5) gives $p = p_0$. At 11 km, $T/T_0 = 216.65/288.15 = 0.7519$, $0.7519^{5.2559} = 0.2234$, $p = 22\,632$ Pa, the tabulated value. Density follows from $p/(RT)$: 1.225 kg/m³ at sea level, 0.364 at 11 km. Equations (5) and (6) meet at 11 km with the same pressure, so the model is continuous there; your test checks this. Above 20 km the lapse rate changes sign and these formulas are wrong; the function refuses such altitudes rather than extrapolating.

### 2.6 Worked example: the A320-class numbers

Mass 70 000 kg, $S = 122.6$ m², $V = 230$ m/s, $h = 11\,000$ m. Weight $W = 70\,000 \times 9.807 = 6.865 \times 10^5$ N. Wing loading $W/S = 5600$ N/m²: inside 3–8 kN/m². Density 0.364 kg/m³, so $q = 0.5 \times 0.364 \times 230^2 = 9628$ Pa, and $C_L = 5600/9628 = 0.58$: inside 0.3–0.8 and close to $C_L^* = 0.67$. Speed of sound $\sqrt{1.4 \times 287.05 \times 216.65} = 295$ m/s, Mach 0.78: inside. $C_D = 0.02 + 0.0442 \times 0.58^2 = 0.0349$; $T/W = 0.0349/0.58 = 0.060$; $L/D = 16.6$: inside. Every quantity passes. The numbers are mutually consistent, which is what you would expect, since they are real.

Now the colleague's mistake: the mass was read from a US document as 154 000 *pounds* and typed as kilograms. Wing loading becomes 12 300 N/m² (flagged: wing far too small for the mass), $C_L$ becomes 1.28 (flagged: approaching a clean wing's maximum). Two flags, both pointing at the mass–area pair. A good checker does not tell you *which* input is wrong, because it cannot; it tells you which *combinations* are impossible, and you supply the judgment.

The finished checker prints the case like this (`python starter/main.py`):

```
== same aircraft, mass typed in lb ==
FLAGGED: 2 quantity(ies) outside plausible range.
ok   mass               1.54e+05 kg    (plausible 20000–600000)
ok   altitude              11000 m     (plausible 6000–14000)
ok   mach                 0.7795       (plausible 0.6–0.92)
FLAG wing_loading     1.232e+04 N/m²  (plausible 3000–8000)  <- wing far too small for the mass; approach speeds would be unflyable
FLAG lift_coefficient      1.28       (plausible 0.3–0.8)  <- C_L above the cruise range and approaching a clean wing's maximum; too slow, too heavy, or wing too small
ok   thrust_to_weight   0.07224       (plausible 0.04–0.12)
```

Read the two reasons together. "Too heavy" is on both lists; "too slow" is contradicted by the Mach line, which passed. The mass is the suspect.

### 2.7 Sensitivity: how errors grow

Equation (2) says $C_L \propto 1/V^2$. A 10 % error in speed produces a 20 % error in $C_L$; the problem *amplifies* speed errors by a factor of 2. This factor is the relative condition number of the output with respect to that input, the idea behind the lecture on condition numbers in the supplementary reference, and it is worth computing for every quantity you check: it tells you which input to measure most carefully. Mass enters $C_L$ linearly (factor 1), altitude enters through density, which changes about 10 % per kilometre near the tropopause, so a 1 km altitude error is also about a 10 % $C_L$ error.

### 2.8 What a sanity checker cannot do

It cannot detect an error that lands inside the plausible band: a mass 15 % too high is a plausible airliner. It cannot tell which of two inputs is wrong when their combination is implausible. It cannot know that the vehicle is not an airliner; a correct Concorde dataset fails the Mach check. And its ranges embed judgment that goes stale: the 600 t mass limit would have been 400 t before 2005. Say all this in your report. A checker that is presented as a proof of correctness is itself a sanity failure.

## 3. Build the computational model

Open `starter/main.py`. The ISA constants, the `PLAUSIBLE` table with its reasons, the `Flag` and `SanityReport` dataclasses, and `main` are given.

### Step 1: The atmosphere

Implement `isa(altitude_m)` from equations (5) and (6). Refuse altitudes outside −500 m to 20 000 m with `ValueError`.

```python
def isa(altitude_m: float) -> tuple[float, float, float]:
    if not -500.0 <= altitude_m <= H_MAX:
        raise ValueError(f"altitude {altitude_m} m is outside the modelled range")
    exponent = -G0 / (LAPSE * R_AIR)
    if altitude_m <= H_TROPOPAUSE:
        t = T0 + LAPSE * altitude_m
        p = P0 * (t / T0) ** exponent
    else:
        t = T_STRAT
        p11 = P0 * (T_STRAT / T0) ** exponent
        p = p11 * math.exp(-G0 * (altitude_m - H_TROPOPAUSE) / (R_AIR * T_STRAT))
    return t, p, p / (R_AIR * t)
```

`speed_of_sound(T)` is $\sqrt{\gamma R T}$; `dynamic_pressure` is $\tfrac12 \rho V^2$ with the ISA density.

### Step 2: The cruise quantities

`wing_loading` is $mg/S$ (refuse non-positive inputs). `lift_coefficient` divides wing loading by dynamic pressure, equation (2). `mach_number` is $V/a$. `lift_to_drag` implements $C_L/(C_{D0} + kC_L^2)$, and `thrust_to_weight_required` is its reciprocal at the cruise $C_L$, equation (3).

### Step 3: Order of magnitude

`order_of_magnitude(x)` returns the nearest power of ten as an integer: `floor(log10|x| + 0.5)`. So 950 → 3, 20 → 1, 0.002 → −3, and 3 → 0 but 4 → 1 (the boundary is $\sqrt{10} \approx 3.16$). Zero has no order of magnitude; raise.

### Step 4: Flags

`_flag(quantity, value)` looks up the range and reasons in `PLAUSIBLE` and returns a `Flag` with `ok=True` and an empty reason inside the range, or `ok=False` with the low or high reason outside it. `check_airliner` refuses impossible inputs (non-positive, altitude outside the ISA range) with `ValueError`, because those are not implausible values but meaningless ones, then appends flags for mass, altitude, Mach, wing loading, lift coefficient, and thrust-to-weight, in that order.

### Verify

`pytest curriculum/phase-00/lesson-12-engineering-judgment`. Analytic cases: the ISA returns its own defining constants at 0 and 11 km, is isothermal above, and obeys $p = \rho RT$ everywhere; the speed of sound scales as $\sqrt T$. Limiting and scaling cases: $C_L$ falls by exactly 4 when speed doubles and doubles with mass; $L/D \to 0$ both as $C_L \to 0$ and as $C_L \to \infty$. Conservation of consistency: $T/W \times L/D = 1$ by construction. The analytic optimum: sweeping speed, the minimum $T/W$ occurs at $C_L^* = \sqrt{\pi e\,AR\,C_{D0}}$ with value $2\sqrt{C_{D0}k}$, equation (4). Then the judgment tests: three airliner-class datasets pass with no flags; mass in pounds, speed in km/h, and a tenth-size wing are flagged with non-empty reasons.

### Compare

Compare `isa` with the table you built in `P00.03` at 0, 5, 11, and 15 km; they should agree to every printed digit, since both implement the same standard. Then compare your $T/W$ at the A320-class point (0.060) with what you know about real engines: a modern airliner engine pair delivers about 0.25–0.35 of take-off weight at sea level and roughly a quarter of that at cruise altitude, so about 0.06–0.09 of weight is available at cruise. The required 0.060 fits, with margin for climb. The assumed polar is crude, but the answer is the right size, and that is what a sanity check is for.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Two experiments on the A320-class dataset.

**Experiment 1: inject mistakes.** Start from the correct inputs (70 000 kg, 122.6 m², 230 m/s, 11 000 m) and apply, one at a time, five classic unit mistakes: mass in pounds (×2.205), wing area in square feet (×10.76), speed in km/h (×3.6), speed in knots (×1.944), altitude in feet (×3.281). Run `check_airliner` on each. Record, in a table: the mistake, which quantities were flagged, and the reason text. One of the five produces a `ValueError` rather than flags; say which and why that is the correct behaviour. Then design a sixth mistake that the checker does *not* catch (a 15 % error in any input, or two compensating errors), and record that too.

**Experiment 2: the plausible band in speed.** Sweep speed from 120 to 400 m/s at 11 km in 1 m/s steps. For each speed record $C_L$, Mach, $T/W$, and whether the report passes. Plot $T/W$ against speed (save to `tw_vs_speed.png`; use the `Agg` backend) and mark the analytic minimum from equation (4). Report the speed range over which the dataset passes every check, and the speed at minimum $T/W$. Compare that speed with the actual cruise speed: real airliners cruise faster than the minimum-drag speed, and your plot should show why that costs little (the curve is flat near its minimum).

## 5. Validate and improve

Put the hand calculation of §2.6, the checker's output, and the reference table side by side. The $C_L$, wing loading, Mach, and $T/W$ from the code should match the hand values to the number of digits you carried (three). If they differ in the third digit, the likely cause is the density: check that you used 11 000 m, not 10 000 m, and $g = 9.80665$, not 9.81.

Now test the thresholds, which is where the judgment lives. Run the checker on an A380-class dataset (560 000 kg, 845 m², 250 m/s, 11 000 m) and on Concorde-class numbers (185 000 kg, 358 m², 600 m/s, 17 000 m). The first should pass; the second fails on Mach and on $T/W$, and its $C_L$ is low, all of which are *correct descriptions* of a supersonic aircraft and *wrong verdicts* about the data. Decide: do you widen the ranges (fewer false alarms, more missed errors) or keep them and add a vehicle-class argument (the ★ extension)? Write down your choice and its cost. There is no setting that catches every error and raises no false alarm; the cost of each kind is the thing an engineer actually decides.

Finally, estimate the uncertainty of a verdict. The thresholds are round numbers chosen by judgment, so a value within about 10 % of a boundary is neither clearly in nor clearly out. Add a sentence to your report saying which of your flags, if any, were marginal.

## 6. Defend your work

Submit: pytest output, the injected-mistake table with the uncaught sixth mistake explained, `tw_vs_speed.png` with the analytic minimum marked, the A380 and Concorde outputs with your threshold decision, and one sentence naming the kind of error no sanity check will ever find.

### Adversarial self-checks

- What if the mass doubles at fixed wing area, speed, and altitude? Which flags trip, and is the aircraft actually unflyable or just unusual?
- Which assumption fails first for a real aircraft: ISA conditions, the parabolic polar, or steady level flight? What would a 20 K warmer-than-ISA day do to your $C_L$?
- What measurement could falsify the checker? (A real airliner whose certified numbers it flags.)
- Does the verdict survive changing $C_{D0}$ from 0.02 to 0.025? Which flag is most sensitive to the assumed polar?
- Which decision would change if the thresholds were widened by 20 %? Which mistake from Experiment 1 would then slip through?

## Check yourself

1. Give the order of magnitude, as a power of ten, of: an airliner's wing loading in N/m², the ISA density at 11 km in kg/m³, the speed of sound in m/s.
2. Check the dimensions of $q = \tfrac12 \rho V^2$ and name the SI unit it comes out in.
3. If cruise speed is 10 % higher than you were told, by what percentage is the true $C_L$ lower?
4. Why does `check_airliner` raise `ValueError` for an altitude of 30 km instead of flagging it?
5. At what $C_L$ is the required thrust-to-weight ratio smallest, and does that $C_L$ depend on the aircraft's weight?
6. A dataset gives $C_L = 1.3$ in cruise. Name two input errors that could produce it and one physical reason it cannot be right.
7. Why would a dimensional check have passed the Mars Climate Orbiter's small-forces file, and what check would have caught it?
8. A checker flags nothing. List two kinds of error that may still be present.
