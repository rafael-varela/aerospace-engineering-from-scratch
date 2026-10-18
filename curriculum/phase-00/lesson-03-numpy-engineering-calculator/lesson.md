# P00.03  NumPy as an Engineering Calculator

> Write the equation once; let the array carry the loop.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the pressure–altitude relation for a gradient layer and for an isothermal layer from the hydrostatic equation and the ideal gas law. 2. Implement and test a vectorized ISA (temperature, pressure, density, speed of sound) for 0–86 km that accepts a scalar or an array of any shape. 3. Use broadcasting to build a two-parameter table (dynamic pressure against altitude and speed) without writing a loop. 4. Decide, from a stated residual, whether a table of atmospheric properties is trustworthy enough for a performance calculation. |
| **Prerequisites** | `P00.01` (working Python with NumPy), `P00.02` (you can commit your work). Skill check: you can rearrange an equation like p = ρRT for ρ, and you can evaluate xʸ and eˣ on a calculator. No calculus is needed; the one integral is shown in full. |
| **Engineering question** | What are the temperature, pressure, density and speed of sound at any altitude from sea level to 86 km, how fast can we compute them for a million points, and how do we know the numbers are right? |
| **Assumptions and limits** | Dry air as an ideal gas with constant composition; hydrostatic equilibrium (no vertical acceleration); gravity constant at g₀ by working in geopotential altitude; the temperature profile is the ISA's fixed piecewise-linear profile. The ISA is a reference, not a forecast: on any real day temperature can differ by ±30 K and sea-level pressure by ±5 %. The model is undefined below 0 m and above 86 km geometric, and it ignores humidity, latitude and season. |
| **Analytical work** | Hydrostatic equation from a force balance on a slab of air; integration in a gradient layer and in an isothermal layer; the geopotential conversion; hand calculation at 5 km and 15 km; unit check of every result. |
| **Computational deliverable** | `starter/main.py`: `temperature`, `pressure`, `density`, `speed_of_sound`, `isa`, `isa_table`, all vectorized. Tests prove the sea-level and 11 km defining values, continuity at every layer boundary, monotonic pressure and density, the ideal gas law, a hydrostatic residual below 10⁻⁵, the √T scaling of the speed of sound, shape preservation and range checking. |
| **Practical activity** | **computational experiment**: a sensitivity study of density to a hot day (T₀ + 15 K) and a cold day (T₀ − 15 K), and a timing comparison of a Python loop against the vectorized code. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Your completed `main.py`, passing tests, the printed ISA table, the hand calculation at 5 km and 15 km with its comparison to the code, the sensitivity table, and a short note on the timing. |
| **Assessment** | `assessments/phase-00/grading-rubric.md`; competencies 1–4. Pass: derivation correct with assumptions named; all public tests pass; the hand calculation agrees with the code to four significant figures; the judgment questions in the problem set answered with a measurable sign. |
| **References** | `@numpy-userguide` "NumPy: the absolute basics for beginners" (required-free); `@numpy-userguide` "NumPy quickstart" (required-free); `@numpy-userguide` "NumPy fundamentals" > "Broadcasting", "Data types", "Copies and views" (required-free); `@numpy-userguide` "NumPy for MATLAB users" (supplementary-free, only if you know MATLAB). |
| **Extensions** | ★ Extend the model with the 1976 upper-atmosphere layers (86–1000 km), where temperature is no longer piecewise linear in altitude and composition changes; ★ add a `pint` unit interface after `P00.05`. |
| **Time** | unmeasured |

## 1. The engineering problem

Every aircraft and rocket calculation starts with "what is the air doing here?" Lift and drag scale with density. The Mach number, which decides whether the flow is subsonic or supersonic, is the speed divided by the local speed of sound, and the speed of sound depends on temperature alone. An altimeter in a cockpit is a pressure gauge with an altitude scale painted on it, and the painting assumes one particular relationship between pressure and height. That relationship is the **International Standard Atmosphere (ISA)**: an agreed, fixed description of how temperature, pressure and density vary with altitude in a dry, calm, average atmosphere. Nobody has ever flown through the ISA; it is a yardstick, and every real day is measured as a deviation from it ("ISA + 15" means a day 15 K warmer than standard at that altitude).

The engineering question of this lesson is: given any altitude between sea level and 86 km, what are the ISA temperature, pressure, density and speed of sound, and how do we know the numbers are right? The deliverable is a Python function you will use in almost every later lesson, from the drag polar in Phase 10 to rocket trajectories in Phase 13. Those lessons will call it with arrays of a million altitudes, so it has to be fast as well as right. That is the second half of the lesson: **NumPy**, the array library that lets you write the equation once and apply it to a million numbers at the speed of compiled code, instead of one number at a time in a Python loop.

Success criteria, written down before we start: the code reproduces the ISA defining values at sea level (288.15 K, 101 325 Pa, 1.225 kg/m³, 340.29 m/s) and at the 11 km tropopause (216.65 K, 22 632 Pa, 0.3639 kg/m³) to four significant figures; temperature, pressure and density are continuous at every layer boundary; the hydrostatic equation is satisfied to one part in 10⁵ everywhere; and the vectorized code evaluates 10⁵ altitudes in under 10 ms. If any of these fail, the atmosphere table is not fit to build on.

## 2. Understand and derive

### 2.1 The hydrostatic equation

Consider a thin horizontal slab of still air with area A, thickness dh and density ρ at height h (Figure 1). Pressure p pushes up on its bottom face; pressure p + dp pushes down on its top face; gravity pulls down on its mass ρ A dh.

![A thin horizontal slab of air of area A, thickness dh and density ρ. Pressure p + dp pushes down on its top face, pressure p pushes up on its bottom face, and its weight ρ g A dh acts downward; height h increases upward.](figures/hydrostatic-slab.svg)
*Figure 1. Forces on a still slab of air. The bottom face is pushed up harder than the top face is pushed down, and that difference, −dp A, is what holds up the weight ρ g A dh.*

**Assumption 1 (hydrostatic equilibrium):** the slab does not accelerate vertically. Then the forces balance:

p A − (p + dp) A − ρ g A dh = 0   ⟹   **dp = −ρ g dh**.   (1)

Units: ρ g dh has units (kg/m³)(m/s²)(m) = kg/(m s²) = N/m² = Pa. Good. The minus sign says pressure falls as you go up, because there is less air above you to hold up.

### 2.2 The ideal gas law

**Assumption 2 (ideal gas, fixed composition):** dry air obeys

**p = ρ R T**,   (2)

with R = 287.05287 J/(kg K) the specific gas constant of air (the universal gas constant divided by the molar mass of air, 28.9644 g/mol). Units: (kg/m³)(J/(kg K))(K) = J/m³ = N m/m³ = N/m² = Pa. Good. The ISA keeps composition fixed up to 86 km, which is realistic because turbulent mixing keeps the lower atmosphere well stirred; above that the gases begin to separate by molar mass and the model stops.

Substitute (2) into (1) to eliminate ρ:

dp/p = −(g/(R T)) dh.   (3)

This is the equation we integrate. The left side is dimensionless; so is the right: (m/s²)/((J/(kg K))(K)) × m = (m/s²)/(m²/s²) × m = 1. Good.

### 2.3 Geopotential altitude: how we keep g constant

Gravity weakens with height: g(z) = g₀ (R_E/(R_E + z))² for a spherical Earth of radius R_E, where z is the **geometric altitude** (tape-measure height above mean sea level). At 86 km that is a 2.7 % reduction, which is not negligible. Rather than carry g(z) through every integral, the ISA defines a **geopotential altitude** h such that g₀ dh = g(z) dz. Integrating with R_E = 6 356 766 m:

**h = R_E z / (R_E + z)**.   (4)

In words: h is the height at which a constant-g₀ world would give the same potential energy as the real world does at z. The ISA's layers are defined in h. The conversion matters a little: at z = 11 000 m, h = 10 981 m, 19 m lower. Our functions take h; `isa_table` does the conversion for you from z.

### 2.4 The temperature profile is *defined*, not derived

The ISA fixes temperature as piecewise linear in h, with a **lapse rate** L (K/m) in each of seven layers:

| layer | base h_b (m) | lapse L (K/m) | base T_b (K) | name |
|---|---|---|---|---|
| 1 | 0 | −0.0065 | 288.15 | troposphere |
| 2 | 11 000 | 0 | 216.65 | tropopause / lower stratosphere |
| 3 | 20 000 | +0.0010 | 216.65 | stratosphere |
| 4 | 32 000 | +0.0028 | 228.65 | stratosphere |
| 5 | 47 000 | 0 | 270.65 | stratopause |
| 6 | 51 000 | −0.0028 | 270.65 | mesosphere |
| 7 | 71 000 | −0.0020 | 214.65 | mesosphere |
| top | 84 852 | | 186.95 | (= 86 000 m geometric) |

The base temperatures are not independent data: each follows from the one below, T_{b+1} = T_b + L_b (h_{b+1} − h_b). Only T₀ = 288.15 K and p₀ = 101 325 Pa are *defined*; your code derives everything else. Inside layer b,

**T(h) = T_b + L_b (h − h_b)**.   (5)

### 2.5 Integrating (3) in a gradient layer (L ≠ 0)

From (5), dh = dT / L. Substitute into (3):

dp/p = −(g₀/(R L)) dT/T.

Both sides are now "d(something)/something", and ∫ dx/x = ln x. Integrate from the layer base (p_b, T_b) to (p, T):

ln(p/p_b) = −(g₀/(R L)) ln(T/T_b)   ⟹   **p = p_b (T/T_b)^(−g₀/(R L))**.   (6)

The exponent is dimensionless, as it must be: (m/s²)/((J/(kg K))(K/m)) = (m/s²)/(m²/(s² m)) = 1. In the troposphere, −g₀/(R L) = 9.80665/(287.05287 × 0.0065) = **5.2559**, a number worth remembering: tropospheric pressure falls as the 5.26th power of the temperature ratio.

### 2.6 Integrating (3) in an isothermal layer (L = 0)

Now T = T_b is constant, so (3) integrates directly:

ln(p/p_b) = −g₀ (h − h_b)/(R T_b)   ⟹   **p = p_b exp(−(h − h_b)/H)**, with **H = R T_b / g₀**.   (7)

H is the **scale height**, the altitude rise that divides pressure by e ≈ 2.718. At T_b = 216.65 K, H = 6 342 m. (You cannot use (6) here: with L = 0 the exponent is infinite and the base is 1; the limit is exactly (7), but a computer evaluating (6) would divide by zero. This is why the code needs two branches.)

### 2.7 Density and speed of sound

With p and T known, (2) gives **ρ = p/(R T)**. The speed of sound in an ideal gas is

**a = √(γ R T)**,   (8)

with γ = 1.4 the ratio of specific heats of air (a fact about diatomic gases from Phase 05; take it as given here). Units: √((J/(kg K))(K)) = √(m²/s²) = m/s. Good. At 288.15 K, a = √(1.4 × 287.05287 × 288.15) = 340.29 m/s. Notice that a depends only on T: in the isothermal layers it is constant, and above 20 km it *rises* with altitude because temperature does.

### 2.8 Two worked examples by hand

**At h = 5 000 m (troposphere, layer 1).** T = 288.15 − 0.0065 × 5000 = 255.65 K. p = 101 325 × (255.65/288.15)^5.2559 = 101 325 × 0.53314 = 54 020 Pa. ρ = 54 020/(287.05287 × 255.65) = 0.7361 kg/m³. a = √(1.4 × 287.05287 × 255.65) = 320.5 m/s. Sanity check: about half the sea-level pressure at 5 km is the rule of thumb every pilot knows, and density is 60 % of sea level, which is why a piston engine loses 40 % of its power there.

**At h = 15 000 m (isothermal, layer 2).** First the layer base: T₁ = 216.65 K and p₁ = 101 325 × (216.65/288.15)^5.2559 = 22 632 Pa (this *is* the tropopause check value; it is derived, not looked up). Then H = 287.05287 × 216.65/9.80665 = 6 341.6 m and p = 22 632 × exp(−4000/6341.6) = 22 632 × 0.5322 = 12 045 Pa. ρ = 12 045/(287.05287 × 216.65) = 0.1937 kg/m³. Sanity check: at 15 km, roughly the ceiling of a business jet, there is 12 % of sea-level pressure left.

### 2.9 Where the model stops being valid

- **Real weather.** On a hot afternoon the surface can be 40 °C, 25 K above ISA; density is then 8 % lower than the table says and a runway that is long enough on paper may not be. The ISA is the reference, never the measurement.
- **Above 86 km geometric.** Composition changes, the mean free path grows, and temperature is no longer piecewise linear. Our code raises `ValueError` there rather than extrapolating silently.
- **Below sea level.** Dead Sea airports exist (about −400 m). The formulas extrapolate smoothly, but we keep the stated domain and raise.
- **Humidity.** Water vapour is lighter than air; humid air at the same p and T is up to 1–2 % less dense. Ignored here.
- **Hydrostatic assumption.** Inside a thunderstorm updraft, air accelerates vertically at several m/s² and (1) fails locally.

## 3. Build the computational model

Open `starter/main.py`. The constants and the layer table are already there as NumPy arrays. Work through the steps; each block runs on its own in a Python session. The NumPy ideas you need are in `@numpy-userguide` "NumPy: the absolute basics for beginners": an **ndarray** is a grid of numbers of one type (here `float64`, the 64-bit floating-point default), and arithmetic on an array applies to every element at once.

### Step 1: Constants and arrays, not scattered numbers

```python
import numpy as np

T0, P0, G0, R_AIR, GAMMA = 288.15, 101325.0, 9.80665, 287.05287, 1.4
H_BASE = np.array([0.0, 11000.0, 20000.0, 32000.0, 47000.0, 51000.0, 71000.0])
LAPSE = np.array([-6.5e-3, 0.0, 1.0e-3, 2.8e-3, 0.0, -2.8e-3, -2.0e-3])
```

Put each physical constant in one named place with its unit in a comment. The Mars Climate Orbiter was lost because one number in one file was in the wrong unit (`P00.05`); named constants are the first defence.

### Step 2: Geopotential conversion

```python
R_EARTH, Z_TOP = 6356766.0, 86000.0

def geometric_to_geopotential(z):
    z = np.asarray(z, dtype=float)       # accept a float, a list or an array
    h = R_EARTH * z / (R_EARTH + z)      # equation (4), element by element
    return float(h) if h.ndim == 0 else h

H_TOP = geometric_to_geopotential(Z_TOP)   # about 84 852 m
```

`np.asarray` turns whatever you were given into an array without copying if it already is one. A 0-dimensional array (`ndim == 0`) is how NumPy represents a single number; we hand a plain `float` back in that case so callers can print it without surprises.

### Step 3: Derive the layer bases with a small loop

```python
def layer_bases():
    n = len(H_BASE)
    t_base, p_base = np.empty(n), np.empty(n)
    t_base[0], p_base[0] = T0, P0
    for i in range(n - 1):
        dh = H_BASE[i + 1] - H_BASE[i]
        t_base[i + 1] = t_base[i] + LAPSE[i] * dh
        if LAPSE[i] == 0.0:
            p_base[i + 1] = p_base[i] * np.exp(-G0 * dh / (R_AIR * t_base[i]))       # (7)
        else:
            p_base[i + 1] = p_base[i] * (t_base[i + 1] / t_base[i]) ** (-G0 / (R_AIR * LAPSE[i]))  # (6)
    return t_base, p_base
```

A loop over seven layers is fine: vectorization matters where the count is large (a million altitudes), not where it is seven. Each base depends on the previous one, which is exactly what a loop expresses. Run it and compare with the table in §2.4: `p_base[1]` should print `22632.04`.

### Step 4: Which layer is each altitude in?

```python
def _prepare(h):
    arr = np.asarray(h, dtype=float)
    scalar = arr.ndim == 0
    arr = np.atleast_1d(arr)
    if np.any(arr < 0.0) or np.any(arr > H_TOP + 1e-6):
        raise ValueError(f"geopotential altitude must lie in [0, {H_TOP:.1f}] m")
    idx = np.searchsorted(H_BASE, arr, side="right") - 1
    return arr, np.clip(idx, 0, len(H_BASE) - 1), scalar
```

`np.searchsorted(H_BASE, arr, side="right")` answers, for every altitude at once, "after which layer base would I insert this value?"; subtracting 1 gives the index of the layer the altitude is in. The result `idx` is an *integer array the same shape as `arr`*. Writing `LAPSE[idx]` then picks the lapse rate for every altitude in one operation: this is **fancy indexing**, and it is the vectorized replacement for "look up the row in a loop". The `1e-6` m tolerance lets `H_TOP` itself through after floating-point rounding. Invalid input raises rather than returning nonsense; a silent `nan` at 90 km would propagate into a trajectory and be noticed, if ever, hours later.

### Step 5: Temperature, vectorized

```python
def temperature(h):
    arr, idx, scalar = _prepare(h)
    t_base, _ = layer_bases()
    t = t_base[idx] + LAPSE[idx] * (arr - H_BASE[idx])      # (5) for every element
    return float(t[0]) if scalar else t
```

One line of arithmetic, no loop, any shape. Try `temperature(np.array([[0, 5000], [11000, 20000]]))` and confirm you get a 2 × 2 array.

### Step 6: Pressure, with both branches and no division by zero

```python
def pressure(h):
    arr, idx, scalar = _prepare(h)
    t_base, p_base = layer_bases()
    lapse = LAPSE[idx]
    t = t_base[idx] + lapse * (arr - H_BASE[idx])
    iso = lapse == 0.0                                  # boolean mask: which points are isothermal
    lapse_safe = np.where(iso, 1.0, lapse)              # never divide by zero, even in the branch we discard
    p_grad = p_base[idx] * (t / t_base[idx]) ** (-G0 / (R_AIR * lapse_safe))      # (6)
    p_iso = p_base[idx] * np.exp(-G0 * (arr - H_BASE[idx]) / (R_AIR * t_base[idx]))  # (7)
    p = np.where(iso, p_iso, p_grad)
    return float(p[0]) if scalar else p
```

This is the vectorized idiom for an `if`: compute **both** branches for every element, then `np.where(mask, a, b)` picks element-wise. It costs twice the arithmetic, which is nothing, but it has a trap: both branches are evaluated *everywhere*, so the gradient formula is evaluated at isothermal points where L = 0. Without `lapse_safe` you would divide by zero, get `inf` and a warning, and `np.where` would discard it, but a warning in your output is a smell that hides real warnings later. Replacing the zero with a harmless 1.0 in the branch you will throw away is the standard fix.

### Step 7: Density, speed of sound, and the convenience wrappers

```python
def density(h):
    return pressure(h) / (R_AIR * temperature(h))          # (2)

def speed_of_sound(t):
    t = np.asarray(t, dtype=float)
    if np.any(t <= 0.0):
        raise ValueError("temperature must be positive (kelvin)")
    a = np.sqrt(GAMMA * R_AIR * t)                          # (8)
    return float(a) if a.ndim == 0 else a

def isa(h):
    t, p = temperature(h), pressure(h)
    return {"T": t, "p": p, "rho": p / (R_AIR * t), "a": speed_of_sound(t)}

def isa_table(z_geometric):
    z = np.atleast_1d(np.asarray(z_geometric, dtype=float))
    h = geometric_to_geopotential(z)
    out = isa(h)
    return np.column_stack([z, h, out["T"], out["p"], out["rho"], out["a"]])
```

`np.column_stack` glues 1-D arrays side by side into an (n, 6) table. `speed_of_sound` checks for a non-positive temperature because `np.sqrt` of a negative number quietly returns `nan`; the physics says kelvin is positive, so we say so in code.

### Step 8: Broadcasting, the idea that replaces nested loops

Dynamic pressure q = ½ ρ V² drives every aerodynamic load. Suppose you want q for 21 altitudes *and* 6 speeds, a 21 × 6 table:

```python
hs = np.arange(0.0, 20001.0, 1000.0)        # shape (21,)
V = np.arange(50.0, 301.0, 50.0)            # shape (6,)
q = 0.5 * density(hs)[:, None] * V[None, :] ** 2    # shape (21, 6)
```

`density(hs)[:, None]` has shape (21, 1); `V[None, :]` has shape (1, 6). NumPy's **broadcasting rule** (`@numpy-userguide` "Broadcasting"): when two arrays are combined, dimensions of size 1 are stretched to match the other array. The product is therefore (21, 6) and `q[i, j]` is the dynamic pressure at `hs[i]` and `V[j]`, with no loop written. Check one cell by hand: at 10 km and 250 m/s, ρ = 0.4127 kg/m³ and q = 0.5 × 0.4127 × 250² = 12 897 Pa. Mach number is one more line: `M = V[None, :] / speed_of_sound(temperature(hs))[:, None]`.

One more NumPy habit from "Copies and views": slicing an array (`table[:, 3]`) gives a *view* that shares memory with the original. Writing into the slice changes the table. Use `.copy()` when you mean a separate array.

### Verify

Run `python -m pytest curriculum/phase-00/lesson-03-numpy-engineering-calculator` from the repo root. The tests check, in order of how much they would hurt if wrong:

1. **Defining values.** T, p, ρ, a at 0 m and at 11 000 m against the values in §2.4 and §2.8. These are the two numbers every ISA implementation in the world agrees on.
2. **Continuity.** T, p and ρ a hair below and a hair above each layer base agree to 10⁻⁹ relative. A discontinuity would mean a layer base was derived from the wrong neighbour.
3. **Monotonicity.** p and ρ strictly decrease on 5000 points from 0 to H_TOP. The atmosphere has no pressure inversions.
4. **The ideal gas law** holds at 1001 points to 10⁻¹² relative. Trivially true by construction, but it catches a wrong constant in one function and not another.
5. **The hydrostatic equation.** The test differentiates your own `pressure` with a central difference, (p(h+1) − p(h−1))/2, and compares with −ρ g₀. The relative residual must be below 10⁻⁵. This is the real check: it asks the code whether it satisfies the equation it was derived from, independently of any tabulated number. A central difference has truncation error ∝ dh², about 10⁻⁸ here, so a 10⁻⁵ tolerance is generous and any formula slip fails it by orders of magnitude.
6. **Scaling laws.** a(4T) = 2 a(T); ln p is linear in h with slope −g₀/(R T) in the isothermal layer.
7. **Interface.** A float in gives a float out; a (3, 4) array in gives a (3, 4) array out; out-of-range input raises `ValueError`; the vectorized result equals a scalar loop to 10⁻¹².

### Compare

There is no ISA function in NumPy or SciPy to compare against, so we compare two ways. First, against the standard's published table: the layer-base pressures that your `layer_bases()` derives (22 632.0, 5 474.9, 868.02, 110.91, 66.939, 3.9564, 0.3734 Pa) match the tabulated ISA values to the fourth significant figure, with differences at the 10⁻⁵ level caused by the standard's use of a slightly different rounding of R. Second, against a plain Python loop:

```python
import timeit
h = np.linspace(0.0, 80000.0, 100_000)
t_vec = timeit.timeit(lambda: pressure(h), number=5) / 5
t_loop = timeit.timeit(lambda: [pressure(float(x)) for x in h[:2000]], number=1) * 50
print(f"vectorized {t_vec*1e3:.1f} ms, loop {t_loop:.1f} s, ratio {t_loop/t_vec:.0f}")
```

On the author's machine: 5 ms against 3.5 s, a factor of 700. The results agree to 10⁻¹² (that is one of the tests). The loop is not wrong; it is just a thousand times too slow to sit inside a trajectory integration that calls it at every time step.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** No measurement is involved here; we are studying the model's own sensitivity. (If your phone has a barometer, the optional check in Problem 1 compares it to the ISA, and you will meet that properly in `P00.14`.)

**Sensitivity of density to a non-standard day.** Airlines and rocket operators care about *density* most, and the thing that moves it day to day is temperature. Build a non-standard atmosphere by shifting the whole temperature profile: T(h) + ΔT with ΔT ∈ {−15, 0, +15} K, keeping pressure at its ISA value (pressure is set by the weight of air above, which does not change much with a surface temperature shift; this is the standard "ISA + ΔT" convention). Then ρ = p/(R (T + ΔT)).

Prediction, written before running: density is inversely proportional to temperature, so +15 K at sea level (288 K) should reduce density by about 15/303 ≈ 5 %, and −15 K should raise it by about 15/273 ≈ 5.5 %.

```python
for dT in (-15.0, 0.0, 15.0):
    h = np.array([0.0, 1500.0, 11000.0])
    rho = pressure(h) / (R_AIR * (temperature(h) + dT))
    print(dT, np.round(rho, 4))
```

Record the table. Then compute the **density altitude**: the ISA altitude at which density equals the hot-day density at your chosen airport elevation. A fine grid and `np.interp` does the inversion (note that `np.interp` needs its x-values increasing, and density *decreases* with h, so reverse both arrays). For an airport at 1 655 m on an ISA + 20 K day, density altitude is about 2 350 m: the aircraft performs as if it were 700 m higher than the map says. Vary ΔT from 0 to 30 K and tabulate.

**Timing experiment.** Repeat the §3 comparison for 10³, 10⁴, 10⁵ and 10⁶ points. Plot nothing yet (`P00.04`); just record the times and the ratio. Notice that the vectorized time is not proportional to n for small n: there is a fixed overhead of a few microseconds per NumPy call, which is why vectorizing a seven-element loop gains nothing.

Raw data format: a plain-text table, one row per (ΔT, h, ρ), and one row per (n, t_loop, t_vec), kept in your lesson notes and committed (`P00.02`).

## 5. Validate and improve

Put the three sources side by side:

| quantity at 11 km | hand (§2.8 method) | code | defining value | difference |
|---|---|---|---|---|
| T (K) | 216.65 | 216.65 | 216.65 | 0 |
| p (Pa) | 22 632 | 22 632.04 | 22 632 | < 0.001 % |
| ρ (kg/m³) | 0.3639 | 0.36392 | 0.3639 | < 0.01 % |
| a (m/s) | 295.07 | 295.07 | 295.07 | 0 |

Every difference must be explained, and here each has a cause you can name:

- **Rounding in the exponent.** Using 5.2559 by hand instead of 5.25588 changes p at 11 km by 0.0003 %. Hand calculations with four figures agree with the code to four figures, as expected, and no better.
- **The gas constant.** Older texts use R = 287.058 or 287.04 J/(kg K). The choice shifts p at the tropopause by 0.002 %, which is invisible to any instrument. If your code and a reference table differ at the fifth figure, look at R before anything else.
- **Geometric against geopotential altitude.** At z = 11 000 m geometric, h = 10 981 m and T = 216.77 K, not 216.65 K. If you compare your table to one indexed by geometric altitude and see a 0.12 K difference at the tropopause, this is why, and it is not an error in either.
- **Discontinuities you did not expect.** If the continuity test fails at 20 km only, the usual cause is using the gradient formula with L = 0 somewhere in `layer_bases()`; the `inf` then poisons every base above.

What is still wrong even when everything agrees: the model matches *itself*. Agreement between hand calculation and code proves the code implements the equations; it says nothing about whether the equations describe today's sky. That is the job of a measurement (`P00.14`, `P00.15`) and a real dataset (`P00.08`, where you read a radiosonde profile and plot it against this table).

Improvements worth making now: cache `layer_bases()` so it is not recomputed on every call (it is cheap, but it is called twice per `isa()`); accept altitude in kilometres through a keyword rather than letting callers guess the unit (or, better, wait for `pint` in `P00.05`). Improvements worth *not* making: extending above 86 km by extrapolating the last layer. The air up there is not the same gas.

## 6. Defend your work

Submit: the completed `starter/main.py` with all tests passing; the printed ISA table every 2 km; your two hand calculations (5 km and 15 km) with the comparison to the code; the ΔT sensitivity and density-altitude tables; the timing table; and a one-page note following `projects/_template/report-template.md` (sections 1, 2, 5, 6 and 7 are enough) stating what the model assumes, what you verified, and where it would mislead a user.

### Adversarial self-checks

- **What if the lapse rate doubled?** The tropospheric exponent halves to 2.63 and pressure at 11 km would be 101 325 × (145.15/288.15)^2.63 ≈ 16.7 kPa. Can your code take a modified `LAPSE` array without any other edit? It should: that is why the layer bases are derived, not typed.
- **What if altitude goes to zero? To H_TOP?** At 0 you must get exactly the defining values, not 288.149999. At H_TOP you must get a number, not a `ValueError` caused by rounding in the geopotential conversion. At H_TOP + 1 m you must get a `ValueError`.
- **Which assumption fails first in the real world?** The fixed temperature profile, every single day, by up to 30 K. Hydrostatic equilibrium and the ideal gas law hold to better than 1 % almost everywhere below 86 km; constant composition holds to 0.1 %.
- **What measurement could falsify the model?** A radiosonde (weather balloon) profile: temperature against pressure. The ISA predicts a tropopause at 11 km at 216.65 K; over the tropics the real tropopause is near 17 km and 195 K; over the poles in winter near 8 km. One sounding is enough to show the ISA is an average, not a law.
- **Does the result survive a change in input type?** Pass a Python list, an integer, a 0-d array, a (3, 4) array, a `float32` array. The answers must be identical to 10⁻⁶ relative (`float32` has only 7 significant figures, so not better).
- **Which decision changes if the uncertainty doubled?** Suppose your pressure were uncertain by 1 % instead of 0.5 %. An altimeter calibrated to it would read 80 m wrong at sea level instead of 40 m. For cruise at 11 km nothing changes; for an approach in cloud with 60 m of terrain clearance, everything does. The acceptable error is set by the use, never by the table.

## Check yourself

1. In which units is the lapse rate L given in this lesson's code, and what is the tropospheric value in K per kilometre?
2. Pressure at 11 km is 22 632 Pa. Without a calculator, roughly what is it at 11 km + one scale height (about 17.3 km)?
3. Why does the speed of sound increase between 20 km and 47 km even though the air is getting thinner?
4. `np.where(mask, a, b)` evaluates `a` and `b` where? Everywhere, or only where the mask selects them?
5. `density(h)[:, None] * V[None, :]**2` produces an array of what shape when `h` has 21 elements and `V` has 6?
6. Geopotential altitude at 86 km geometric is 84 852 m. Is geopotential altitude always smaller than geometric altitude, larger, or does it depend?
7. The hydrostatic test uses a central finite difference with dh = 1 m. Why does the test tolerate a residual of 10⁻⁵ rather than demanding 10⁻¹²?
8. Is the sea-level density 1.225 kg/m³ a *defined* value of the ISA or a *derived* one? Which law determines it, from which two defined values?
