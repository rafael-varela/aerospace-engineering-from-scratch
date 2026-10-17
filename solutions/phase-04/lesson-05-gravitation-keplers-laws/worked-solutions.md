# P04.05  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for $s = 15$; numerical results come from `solution.py` and the probe calculations that call it, and are given to the figures shown.

## 1. A circular orbit at 650 km ($s = 15$)

**Assumptions.** Point-mass Earth, μ = 3.986004418 × 10¹⁴ m³/s², R = 6371 km; circular orbit; no drag.

**Derivation and result.** r = 7021 km.
- g(h) = μ/r² = **8.086 m/s²** (82 % of the surface value).
- v = √(μ/r) = **7534.8 m/s**; T = 2πr/v = **97.58 min**.
- ε = −μ/2r = **−28.386 MJ/kg**.
- v_esc = √(2μ/r) = **10 655.8 m/s**; one burn along the velocity needs Δv = v_esc − v = (√2 − 1)v = **3121.0 m/s**.

**Check.** v²/r = 7534.8²/7.021 × 10⁶ = 8.086 m/s² = g(h). Kinetic energy per unit mass v²/2 = 28.386 MJ/kg = μ/2r, half the magnitude of the potential −μ/r = −56.77 MJ/kg, so ε = −28.386 MJ/kg as stated. Units: √(m³ s⁻² / m) = m/s.

**Common mistakes.** Using the altitude instead of the radius in μ/r; using g = 9.81 m/s² at altitude; reporting T in seconds when minutes were asked; forgetting that escape from orbit needs only the *difference* v_esc − v_c.

## 2. Vis-viva, eccentricity, the third law, and the leapfrog

(a) At the apsides v ⊥ r, so h = r_p v_p = r_a v_a and v_a = v_p r_p/r_a. Energy conservation: v_p²/2 − μ/r_p = v_p² r_p²/(2r_a²) − μ/r_a. Collect: (v_p²/2)(1 − r_p²/r_a²) = μ(r_a − r_p)/(r_p r_a). Divide by (r_a − r_p)(r_a + r_p)/r_a²: **v_p² = 2μr_a/(r_p(r_p + r_a))**. Then ε = v_p²/2 − μ/r_p = μr_a/(r_p(r_p + r_a)) − μ(r_p + r_a)/(r_p(r_p + r_a)) = −μ/(r_p + r_a) = **−μ/2a** since r_p + r_a = 2a. Because ε = v²/2 − μ/r at every point, **v² = μ(2/r − 1/a)**.

(b) h² = r_p² v_p² = 2μ r_p r_a/(r_p + r_a) = 2μ a²(1 − e)(1 + e)/(2a) = **μa(1 − e²)**. With a = −μ/2ε: 1 − e² = h²/(μa) = −2εh²/μ², so **e = √(1 + 2εh²/μ²)**. e = 1 is ε = 0: the parabolic escape trajectory.

(c) The areal velocity is constant, h/2, so one period sweeps the whole ellipse: πab = (h/2)T. With b = a√(1 − e²) and h = √(μa(1 − e²)): T = 2πa²√(1 − e²)/√(μa(1 − e²)) = **2π√(a³/μ)**, independent of e.

(d) Kick: v′ = v + (h/2) f r, so r × v′ = r × v + (h/2) f (r × r) = r × v. Drift: r′ = r + hv, so r′ × v = r × v + h(v × v) = r × v. Every substep preserves r × v exactly, so the composition does: the leapfrog **conserves angular momentum exactly** for any central force, limited only by rounding. RK4 updates v with a weighted sum of accelerations evaluated at four different intermediate positions; that sum is not parallel to the r that the update ends at, so r × v changes at the level of the local truncation error, O(h⁵) per step.

**Check.** The code's ten-orbit runs in Problem 4 show |Δh/h| ≤ 2 × 10⁻¹⁴ at every step size from 250 to 4000 per orbit, no trend with h: rounding, not truncation.

**Common mistakes.** In (a), equating speeds rather than angular momenta at the apsides; in (c), using the circle's area πa²; in (d), arguing from "symplectic" alone (symplecticity bounds the energy error; exact angular momentum comes from the central-force structure).

## 3. Mars

**Assumptions.** Point-mass Mars; areostationary means circular, equatorial, period one sidereal day.

(a) r = (μT²/4π²)^{1/3} = **20 427.7 km**, altitude **17 038.2 km**, v = √(μ/r) = **1448.0 m/s**.

(b) g = μ/R² = **3.728 m/s²** (38 % of the Earth's 9.82); v_esc = √(2μ/R) = **5027.0 m/s** (45 % of the Earth's 11 186 m/s from the same formula).

(c) T = 2π√(a³/μ) = **7.657 h**, much shorter than the 24.62 h sidereal day. Phobos goes round faster than the surface turns, so it overtakes the observer: it **rises in the west and sets in the east**, about twice a day.

**Check.** Ratio r_areo/R_Mars = 6.03 against r_GEO/R_Earth = 6.62: similar because both planets have similar rotation periods and Mars' lower mean density is partly offset by its smaller size (r/R ∝ (ρ̄T²)^{1/3}).

**Common mistakes.** Using the solar day (the sol) instead of the sidereal day; mixing km and m inside μ.

## 4. Step size for an energy tolerance ($s = 15$: e = 0.125, a = 8000 km)

(a) Perigee altitude **629 km**, apogee altitude **2629 km**; v_p = **8003.8 m/s**, v_a = **6225.2 m/s**.

(b)

| steps/orbit | position error, 1 orbit (m) | max \|ΔE/E\|, 10 orbits | max \|Δh/h\|, 10 orbits |
|---|---|---|---|
| 250 | 16 030 | 8.769 × 10⁻⁵ | 4.9 × 10⁻¹⁵ |
| 500 | 4007 | 2.191 × 10⁻⁵ | 3.5 × 10⁻¹⁵ |
| 1000 | 1002 | 5.475 × 10⁻⁶ | 1.2 × 10⁻¹⁴ |
| 2000 | 250.5 | 1.369 × 10⁻⁶ | 8.9 × 10⁻¹⁵ |
| 4000 | 62.6 | 3.422 × 10⁻⁷ | 1.9 × 10⁻¹⁴ |

**Observed order 2.00.**

(c) N = 1000 × √(5.475 × 10⁻⁶ / 10⁻⁶) = **2340 steps per orbit**; the check run at 2340 gives **9.999 × 10⁻⁷**. The h² law is exact to four figures here because the leapfrog's energy error is a smooth function of h²: there is no pre-asymptotic regime at these steps.

(d) Simulated period 7121.097 s against Kepler 7121.082 s: **+0.016 s**, so the simulated satellite is **late**: the leapfrog's phase error makes the numerical orbit slightly slow, and the one-orbit position error is mostly along-track.

**Common mistakes.** Measuring the energy error over one orbit and calling it the ten-orbit value without checking (here they are equal, but you must show it); using N ∝ error ratio instead of its square root; comparing positions at mismatched times (interpolate the leapfrog onto the oracle's times or vice versa).

## 5. Transfer to the GPS radius ($s = 15$: parking orbit 315 km)

(a) r₁ = 6686 km, r₂ = 26 560 km, a_t = 16 623 km. **Δv₁ = 2038.7 m/s, Δv₂ = 1417.1 m/s, total 3455.8 m/s; time of flight 2.962 h; e_t = 0.5978.**

(b) After the first burn the state is (r₁, 0) with speed v_c(r₁) + Δv₁ along +y; `orbital_elements` returns r_a = **26 560.000 km** (to 10⁻⁹ m relative). Independent of `hohmann_transfer`: it uses only ε and h of the post-burn state.

(c) A plane change of Δi at speed v costs 2v sin(Δi/2), proportional to the speed. The speed is lowest at the transfer's apoapsis (2457 m/s there by vis-viva, against 3874 m/s in the final circular orbit), so the cheapest plan combines the plane change with Δv₂ at apoapsis: one burn whose vector is the sum of the circularisation and the rotation, smaller than the two done separately.

**Common mistakes.** Using altitude 26 560 km (it is a radius); adding the plane change as a separate burn in the parking orbit, where the speed is highest.

## 6. Judgment

(a) **Not adequate.** J₂ regresses the ISS orbit plane by 5.0° per day (about 30 km of ground-track shift per orbit at the equator crossing), and drag changes the along-track position by kilometres per day; both must be added, and the drag coefficient and atmospheric density are the dominant uncertainty, not the integrator.

(b) **Adequate.** The Hohmann closed form gives 3856.7 m/s; J₂ and drag change it negligibly for a feasibility study. The integrator is not even needed.

(c) **Not adequate.** Over ten years lunar and solar perturbations tilt the orbit and solar radiation pressure changes its eccentricity; these are what station-keeping corrects, so they are the whole point. Add third bodies and radiation pressure; the leapfrog's bounded energy error is an advantage over ten years, but its fixed step must resolve the lunar period too.

(d) **The physics is adequate (for one return, ignoring planets); the numerics are not.** At e = 0.99 a fixed step that resolves perihelion wastes nearly all its steps far out; use a time transformation (dt ∝ r^{3/2}) or an adaptive solver with tight tolerance, and include Jupiter, which perturbs comet periods by months.

(e) **Adequate, and the answer is "the same".** The acceleration is μ/r², independent of the satellite's mass (and the satellite's own mass changes μ = G(M + m) by 10⁻²¹).

(f) **Adequate numerically, and in fact the right choice.** The leapfrog's energy error is bounded over 10⁵ orbits while RK4's drifts linearly; only the energy distribution matters, so the phase error is irrelevant. The physics (J₂, drag) still needs adding for a real cloud.

## Check-yourself answers

1. g(h) = μ/(R + h)². At 400 km, 88.5 % of the surface value (8.694 against 9.820 m/s²).
2. v_c = √(μ/r), v_esc = √(2μ/r); ratio √2.
3. ε = −μ/2a. ε = +1 MJ/kg is positive: unbound (it escapes on a hyperbola).
4. v² = μ(2/r − 1/a). At r = a it gives the circular speed √(μ/a); as a → ∞ it gives the escape speed √(2μ/r).
5. Conservation of angular momentum (a central force exerts no torque). Each kick changes v along r and each drift moves r along v; neither changes r × v, so the leapfrog keeps h exact and the areal velocity constant.
6. Neither: the period depends only on a (Kepler's third law), so they are equal.
7. Thrust raises the energy ε = −μ/2a, so a grows; in the higher orbit the circular speed √(μ/r) is lower. The added energy goes into potential energy, and the satellite gives back kinetic energy as it climbs.
8. A monotone growth means the integrator is not symplectic as implemented: typically the half-kicks are wrong (a full kick at the start, or the acceleration evaluated at the old position at the end), which turns the leapfrog into a first-order, non-symplectic method.
