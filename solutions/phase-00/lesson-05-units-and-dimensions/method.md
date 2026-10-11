# P00.05  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. For each conversion write the chain of unity fractions: `ft → m` uses 0.3048 m/ft (exact); `knot → m/s` uses 1852 m/nmi and 3600 s/h; `inHg → Pa` uses the conventional inch of mercury (3386.389 Pa, which `pint` defines); `slug/ft³` uses 1 slug = 1 lbf·s²/ft = 14.5939 kg and (0.3048 m)³. 2. Multiply, cancelling units. 3. Compute the same with `Q_(value, unit).to(target)`. 4. Report the relative difference; it should be at the level of your hand rounding.

## 2.
1. (a) Multiply the two defining constants; count significant figures; note that a product of two exact numbers is exact. 2. (b) Write `F = d(mv)/dt`; assume `m` constant during the burn so `F dt = m dv`; integrate to `I = m Δv`. Quantify the failure: propellant mass expelled during the firing as a fraction of `m`; if it exceeds ~1 % the constant-mass assumption is off by that order and the rocket equation is needed. 3. (c) `lbf·s / lb = (lb · g₀ · s) / lb = g₀ · s`, a velocity; evaluate `g₀ × 1 s`. Explain `I_sp` in seconds as "impulse per unit weight" and its advantage of being system-independent.

## 3.
1. Mass = density × volume. 2. Convert 6.7 lb/gal to kg/L for each gallon definition, or convert the volume to litres each way; multiply. 3. The code: `rho = Q_(6.7, "lb / gallon")`, `V = Q_(2400, "gallon")`; a bare float `V` makes `rho * V` a density, and `.to("kg")` raises. 4. The imperial gallon is 20 % larger than the US gallon; compare 20 % with 2 %.

## 4.
1. `rng = np.random.default_rng(seed)`, `x = rng.uniform(0.2, 0.6, n)`. 2. `dv = accumulated_delta_v_error(x, Q_(m, "kg"))`. 3. Position error `≈ dv × 30 days`, converted to km with `pint`. 4. State that firing size, count and mass are illustrative and that a scalar sum ignores direction and timing; compare the order of magnitude with 170 km.

## 5.
1. `v = Q_(np.array([180, 200, 220, 240, 260.0]), "knot")`. 2. `q = dynamic_pressure(Q_(0.3639, "kg/m**3"), v)`; `.to("psi")` for the second column. 3. Plot `magnitude_in(v, "knot")` against `magnitude_in(q, "Pa")` with `new_figure("True airspeed", "knot", "Dynamic pressure", "Pa")`. 4. Repeat with `"km/h"`; the ratio of the two `q` curves is the square of the knot/(km/h) ratio.

## 6.
1. For each case ask: do the two quantities differ in dimension? If yes, `pint` raises. If no, name the human or structural check (a frame/convention field in the data, a code review, a definition carried with the number). 2. For temperature, consult the pint documentation on offset units: absolute temperatures cannot be added, differences (`delta_degC`) can.

## 7.
1. Read lines; require the first to match `^#\s*impulse\s*\[(.+)\]\s*$`. 2. Build the unit with `ureg(match.group(1))` and the array with `np.loadtxt` or a float loop that raises `ValueError` naming the bad line. 3. Return `Q_(values, unit).to("N*s")`; the `.to` raises `DimensionalityError` for a header like `[lbf]`. 4. Tests with `tmp_path`: write three small files, assert the factor, assert the two exceptions. 5. Paragraph: the format fixes the root cause's mechanism; it does not fix the contributing causes (unfamiliarity, staffing, anomalies not escalated), which would have let a different mismatch through. Sufficient for *this* error, not for the class.
