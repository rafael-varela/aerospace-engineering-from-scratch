# P08.01  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
(a) Δp = (ρ_m − ρ_f) g Δh with `manometer_delta_p`; compare with ρ_m g Δh. (b) Length = Δh/sin 10°; resolution = `manometer_delta_p(0.001, angle_deg=...)` for 90° and 10°. (c) `tank_bottom_pressure(p_u, 810, 1.5, accel)` with accel = 0, 3g, −g. Say what happens to the liquid at the outlet when the effective gravity is zero.

## 2.
(a) Forces on a slab of area A, thickness dh: p A − (p + dp) A − ρ g A dh = 0; no shear on the side faces because the fluid is at rest. (b) Separate variables, dp/p = −g₀ dh/(R T(h)); for a ≠ 0 substitute T; for a = 0 integrate directly; for the limit, write (T/T_b)^(−g₀/(aR)) = exp[−(g₀/(aR)) ln(1 + a(h − h_b)/T_b)] and expand the logarithm. (c) p(h) ≈ p₀ + p′(0) h + ½p″(0)h² with p′(0) = −ρ₀g₀; set p₀ − ρ₀g₀h = 0. (d) ρ = p/(RT) with p from (b).

## 3.
(a) R = R_u/M; ρ = p/(RT); H = RT/g; compare with `scale_height(288.15)`. (b) p = p_s exp(−z/H). (c) `brentq(lambda h: isa_density(h) − ρ_Mars, 0, 47000)`. (d) Re = ρVL/μ: same Mach means similar tip speed, ρ is ~1/80 of Earth's, μ is similar.

## 4.
(a) p = `isa_pressure(3000)`, T = `isa_temperature(3000)` + ΔT, ρ = p/(RT); compare with `isa_density(3000)`. (b) Root of `isa_density(h) − ρ`. (c) `integrate_hydrostatic(4000, 10, T_of_h=...)`, then interpolate the height where p equals the ISA pressure at 3000 m; closed form: invert p₀[(T₀ + ΔT + a h)/(T₀ + ΔT)]^(−g₀/(aR)) = p. (d) Match each altitude to what depends on it: height above ground, density, pressure.

## 5.
(a) ρ_air = `isa_density(h)`, ρ_He = p/(R_He T), V = m/(ρ_air − ρ_He), r = (3V/4π)^(1/3). (b) `pressure_force_on_sphere(isa_pressure, h, r, 400)`, `displaced_weight(isa_density, h, r)`, ρgV; then errors for n = 25, 50, 100, 200 and `observed_order`. (c) Same helium mass, sea-level density: V = m_He/ρ_He,0; lift = (ρ₀ − ρ_He,0) V. (d) ρ(h + r)/ρ(h − r).

## 6.
(a) h = `geopotential_altitude(47000)`; compare `isa_pressure` at both values. (b) Solve p₀ − ρ₀g₀h = (1 − ε) `isa_pressure(h)` for ε = 0.01, 0.05. (c) Scale-height argument. (d) Quantify each at 30 km: the constant-density model, the geopotential confusion, and a +10 K day (integrate with `T_of_h` offset and find the ISA altitude of the resulting pressure).
