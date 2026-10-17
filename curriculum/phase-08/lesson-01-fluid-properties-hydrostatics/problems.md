# P08.01  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical. Draw the fluid element or the column you are balancing before you write an equation.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

Shared numbers: g₀ = 9.80665 m/s², R = 287.05287 J/(kg·K) for air, water at 20 °C ρ = 998.2 kg/m³, air at 20 °C ρ = 1.204 kg/m³. Propellant densities are typical values, not specifications: RP-1 (kerosene) 810 kg/m³.

## 1. (Warm-up: a manometer and a propellant tank)

(a) A U-tube manometer filled with water, with air above the water in both legs, reads Δh = 20 + s mm between the two water surfaces. Find the pressure difference. By what fraction does the answer change if you ignore the air's density?
(b) The same pressure difference is read on a manometer inclined at 10° to the horizontal. How long is the reading along the tube? If you can read either scale to ±1 mm, what is the pressure resolution of each?
(c) An RP-1 tank on a small launcher has 1.5 m of liquid above the outlet and an ullage (gas) pressure of 2.5 bar. Find the outlet pressure on the pad, during a burn at a total axial acceleration of 4g (thrust acceleration 3g above gravity's 1g, so `accel` = 3g in `tank_bottom_pressure`), and in a coast with the engines off. Why do upper stages fire small "ullage motors" before relighting the main engine?

Answer form: a pressure and a fraction; a length and two resolutions; three pressures and a sentence.

## 2. (Derivation: the standard atmosphere by hand)

(a) Starting from the force balance on a horizontal slab of fluid at rest, derive dp/dh = −ρg. Say where you used the fact that a fluid at rest carries no shear stress.
(b) With p = ρRT and T = T_b + a(h − h_b), integrate to show p = p_b (T/T_b)^(−g₀/(aR)) in a gradient layer, and p = p_b exp(−g₀(h − h_b)/(RT_b)) in an isothermal layer (a = 0). Show that the first is the a → 0 limit of the second.
(c) Expand the ISA pressure about sea level to second order in h and show that the constant-density result p₀ − ρ₀g₀h is the first-order part. Show that the height at which the constant-density atmosphere would reach zero pressure equals the scale height RT₀/g₀, and evaluate it.
(d) Show that the density in the troposphere varies as (T/T₀)^(−g₀/(aR) − 1) and evaluate the exponent.

Answer form: four derivations, the final formulas boxed, two numbers.

## 3. (Unfamiliar scenario: a helicopter on Mars)

The Martian atmosphere is mostly carbon dioxide (molar mass 44.01 g/mol, so R = 8.3145/0.04401 J/(kg·K)). Near the surface take T = 210 K (treat the lower atmosphere as isothermal), p = 610 Pa (a typical surface value) and g = 3.71 m/s².

(a) Compute R, the density at the surface and the scale height. Compare the scale height with Earth's at 288 K.
(b) Estimate the pressure 10 km above the datum and at the summit of Olympus Mons, about 21.9 km above it.
(c) A helicopter needs a certain lift. At what ISA altitude on Earth is the air as thin as at the Martian surface? Use `isa_density` and a root finder.
(d) The Martian helicopter's rotor tips ran at a Mach number close to the Earth helicopter's. Why does that make its Reynolds number much lower? (One sentence; `P08.05` makes this quantitative.)

Answer form: three numbers and a comparison; two pressures; an altitude; a sentence.

## 4. (Numerical, uses your code: a hot day and density altitude)

The light aeroplane flies at a pressure altitude of 3000 m (the altimeter, set to 1013.25 hPa, reads 3000 m, so the static pressure is the ISA value at 3000 m) on a day that is ISA + ΔT with ΔT = 10 + s/2 K at every height.

(a) Compute the outside air temperature and the density at that pressure. By how many per cent is the density below the ISA value at 3000 m?
(b) Find the **density altitude**: the ISA altitude with the same density. Engine power and wing lift depend on density, so performance charts use it.
(c) Find the true (geopotential) height of the aeroplane: integrate dp/dh = −ρg₀ upward from the standard sea-level pressure with `integrate_hydrostatic(..., T_of_h=lambda h: isa_temperature(h) + ΔT)` and find where p reaches the ISA pressure at 3000 m. Check it against the closed form of Problem 2(b) with T₀ + ΔT at sea level.
(d) Which of the three altitudes matters for terrain clearance, which for take-off distance, and which for air-traffic separation?

Answer form: a temperature, a density and a percentage; an altitude; an altitude by two routes; three short answers.

## 5. (Numerical, uses your code: a stratospheric balloon)

A zero-pressure helium balloon (open at the bottom, so the gas inside is at the outside pressure and, we assume, the outside temperature) must float a payload plus envelope of 50 kg at h = 20 000 + 100s m. Helium: R = 2077.1 J/(kg·K).

(a) Compute the air and helium densities at float altitude, the volume needed and the radius of a spherical envelope.
(b) Check Archimedes' principle for that sphere with `pressure_force_on_sphere` (n = 400) against `displaced_weight` and against ρgV, and show that the surface integration converges at order 2 as n doubles from 25.
(c) The same helium mass is loaded at sea level. What volume does it occupy there, and what is its net lift? Explain the result.
(d) By how much does the air density change from the bottom to the top of the balloon? Is ρgV a good enough model of the buoyancy?

Answer form: two densities, a volume and a radius; three forces and a convergence table; a volume, a lift and an explanation; a ratio and a sentence.

## 6. (Judgment: which simplification fails first?)

A sounding rocket climbs to 47 km. Its flight software estimates altitude from a static-pressure sensor using the ISA.

(a) The ISA layer boundaries are defined in geopotential altitude, but the rocket's tracking radar measures geometric altitude. At a geometric 47 km, how far apart are the two, and what pressure error would follow from confusing them?
(b) Above which altitude does a constant-density model of air make a 1 % error in pressure? A 5 % error?
(c) On a real day the stratosphere can be 10 K warmer than the ISA. Without computing, say whether the pressure-altitude error that follows grows or shrinks with height, and why (think of the scale height).
(d) Rank the three effects for a pressure-based altimeter at 30 km and say which you would correct first.

Answer form: a distance and a percentage; two altitudes; a sentence; a ranking with a reason.
