# P00.03  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer, and say how many significant figures you trust and why. Altitudes are geopotential unless a problem says "geometric".

## 1. (Warm-up, parameterized) The air where you live

Look up the elevation of the place where you are sitting, z_home in metres above mean sea level (use 500 m if yours is below 100 m, so the numbers are not trivial). Convert it to geopotential altitude with equation (4). Then, by hand with a calculator, compute T, p, ρ and a there using equations (5), (6), (2) and (8) of the lesson. Compare with `isa(h)` from your code. Answer form: a four-row table (hand, code, relative difference) and one sentence on the size of the difference and why it is that size.

Optional: if your phone has a barometer (many do; search its sensor list), read the pressure now and compare. Expect disagreement of up to ±3 %. Write down the two leading causes before you look anything up.

## 2. (Derivation) The isothermal layer and the scale height

Starting from dp/p = −(g₀/(R T)) dh with T constant, derive p(h) = p_b exp(−(h − h_b)/H) and give H in terms of R, T and g₀. State each assumption where you use it. Evaluate H for the 11–20 km layer. Then, using only the tropopause pressure 22 632 Pa and your H, compute the pressure at 20 km and compare with `layer_bases()[1][2]` from your code. Answer form: derivation, H with units, p(20 km) to four significant figures, relative difference from the code.

## 3. (Unfamiliar scenario) A standard atmosphere for Mars

A simplified Mars reference atmosphere, idealized for this exercise, has: surface pressure 610 Pa, surface temperature 210 K, g = 3.71 m/s², a CO₂ atmosphere with R = 188.92 J/(kg K) and γ = 1.3, and a single gradient layer with lapse rate −2.5 K/km from the surface to 40 km. Copy `main.py` to a scratch file, change only the constants and the layer tables, and compute T, p, ρ and a at 10 km. Compute the surface scale height R T₀/g for Mars and for Earth. Answer form: a table of the four properties at 10 km on both planets, both scale heights, and a three-sentence argument for why a parachute that works on Earth needs to be much larger on Mars even though Mars gravity is weaker. Which parts of your Earth code did you *not* have to change, and what does that tell you about how the code is organised?

## 4. (Numerical, uses your code) Density altitude on a hot day

An airport sits at 1 655 m. Take ΔT = 10 + (last two digits of your birth year mod 21) K as today's temperature excess over ISA, keeping pressure at its ISA value. Compute the actual density ρ = p_ISA(h)/(R (T_ISA(h) + ΔT)). Then find the **density altitude**: the ISA altitude whose density equals ρ. Do it numerically with `np.interp` on a fine grid of your `density()` (remember `np.interp` needs increasing x-values), and confirm by evaluating `density()` at the result. Answer form: ρ (kg/m³, 4 s.f.), density altitude in metres with the grid resolution as its uncertainty, the confirmation value, and one sentence on what this means for the aircraft's takeoff distance.

## 5. (Broadcasting) Dynamic pressure and Mach number tables

Let `hs = np.arange(0, 20001, 1000.0)` (21 altitudes) and `V = np.arange(50, 301, 50.0)` (6 speeds, m/s). Without any `for` loop:

(a) Build `q = ½ ρ V²` as a 21 × 6 array using broadcasting and state its shape.
(b) Build the Mach number table `M = V / a(h)` the same way.
(c) List the (altitude, speed) cells where q exceeds 20 kPa and the cells where M exceeds 0.8. Describe the shape of each region in words (which corner of the table, and why).
(d) Verify one cell of each table by hand.

Answer form: the two shapes, the two cell lists (or a Boolean count plus the pattern), and the hand check.

## 6. (Judgment) When the table lies

(a) Give one *measurable* sign that the ISA is unfit for a given day's performance calculation, and name one decision (takeoff, cruise altitude, descent planning, rocket staging) that changes because of it.
(b) A colleague proposes extending your code to 100 km by continuing layer 7 (L = −2.0 K/km) upward. Compute what that extrapolation gives for T and p at 100 km geometric. Argue, in a short paragraph, why it should not be done, considering that the real upper-atmosphere models put the temperature around 100 km near 190–200 K and *rising* with altitude. What makes this extrapolation more dangerous than one that produced an obviously absurd number?
