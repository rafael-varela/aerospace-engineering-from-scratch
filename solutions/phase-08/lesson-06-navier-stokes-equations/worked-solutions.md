# P08.06  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe scripts that call it, given to the figures shown. Fluid properties are typical values, not data-sheet values.

## 1. A turbopump journal bearing (s = 15: N = 9000 rpm)

**Assumptions.** Concentric shaft, steady laminar film, no pressure gradient round the film, constant viscosity, the narrow gap (h/R = 1/600) treated as plane Couette flow.

(a) Ω = 2π × 9000/60 = **942.48 rad/s**; U = ΩR = **14.137 m/s**. The shear stress is uniform across a Couette film: τ = μU/h = 0.03 × 14.137/25 × 10⁻⁶ = **16 965 Pa**. Torque = τ (2πRL) R = **1.199 N·m**; power = torque × Ω = **1130 W**.

(b) Re = ρUh/μ = 870 × 14.137 × 25 × 10⁻⁶/0.03 = **10.25**; Ta = Re² h/R = **0.175**, four orders of magnitude below the Taylor-vortex threshold of about 1700. The laminar Couette solution is valid; the flow is far from any instability.

(c) The 1.1 kW is dissipated as heat in a film only 25 µm thick. It is carried away by the oil flowing axially through the bearing and by conduction into the shaft and housing. If the oil is not replaced fast enough it warms, its viscosity falls (steeply: oil viscosity is a strong function of temperature), the torque falls with it, and so does the film's ability to carry load. The constant-viscosity assumption is the first to fail in a high-speed bearing.

**Independent check.** Dimensional: τ has units Pa·s × (m/s)/m = Pa ✓. The torque from `wall_shear_stress(25e-6, 14.137, 0, 0.03)` gives the same stress at both walls, as Couette flow must.

**Common mistakes.** Using the diameter for R; forgetting the second factor R (stress → force → torque).

## 2. Derivations

(a) **Assumptions:** incompressible Newtonian fluid with constant ρ and μ; steady; fully developed (∂/∂x of every velocity component is zero); no body force along x. Continuity: ∂u/∂x + ∂v/∂y = 0 → ∂v/∂y = 0, so v is constant across the gap, and v = 0 at the walls makes **v = 0**. y-momentum: every velocity term vanishes, leaving **∂p/∂y = 0** (ρg_y would add a hydrostatic part that does not drive the flow). x-momentum: ρ(u u_x + v u_y) = 0, so −∂p/∂x + μ u″ = 0. The left term depends only on x, the right only on y: both are constant, dp/dx = −G. Integrating μu″ = −G twice, u = −Gy²/(2μ) + C₁y + C₂; u(0) = 0 ⇒ C₂ = 0; u(h) = U ⇒ C₁ = U/h + Gh/(2μ):

**u = Uy/h + Gy(h − y)/(2μ),  q = ∫₀ʰ u dy = Uh/2 + Gh³/(12μ),  τ(0) = μU/h + Gh/2,  τ(h) = μU/h − Gh/2.**

τ(h) − τ(0) = −Gh: a slab of length dx feels the net pressure force G h dx forward and the wall stresses τ(h) dx on its top face and −τ(0) dx on its bottom face; their sum is zero ✓.

(b) For w(r) with no radial or swirl velocity, continuity and the convective terms vanish identically and the axial equation is μ(1/r)(r w′)′ = −G. Integrating, r w′ = −Gr²/(2μ) + C₁, so w′ = −Gr/(2μ) + C₁/r; a finite shear on the axis needs **C₁ = 0**. Then w = G(R² − r²)/(4μ) with w(R) = 0, and **Q = ∫₀ᴿ 2πr w dr = πR⁴G/(8μ)**. With V = Q/(πR²) = GR²/(8μ) and τ_w = μ|w′(R)| = GR/2: f = 8τ_w/(ρV²) = 8(GR/2)/(ρV · GR²/(8μ)) = 32μ/(ρVR) = **64μ/(ρVD) = 64/Re** ✓.

(c) With η = y/(2√(νt)): η_y = 1/(2√(νt)), η_t = −η/(2t). Then u_t = UF′(η)(−η/(2t)) and u_yy = UF″/(4νt). Substituting into u_t = νu_yy and multiplying by 4t/U: −2ηF′ = F″, i.e. **F″ + 2ηF′ = 0**. So (ln F′)′ = −2η, F′ = C e^{−η²}, F = 1 − (2/√π)∫₀^η e^{−s²} ds = **erfc η** with F(0) = 1 and F(∞) = 0 (which fixes C = −2/√π). The problem has no length scale: U, ν and t are all it contains, and the only length they make is √(νt). So u/U can depend only on y/√(νt). The initial condition u(y, 0) = 0 and the far condition both map to η → ∞, which is why one ODE condition covers both.

**Common mistakes.** Forgetting that the separation-of-variables argument (x-only = y-only) is what makes dp/dx constant; dropping the 1/r term in the pipe.

## 3. A de-icing fluid film (s = 15: θ = 17.5°)

**Assumptions.** Steady, fully developed (an infinitely long panel), laminar, Newtonian, constant properties, a flat free surface at atmospheric pressure, no air shear.

(a) x down the slope, y normal to the panel. Body force per unit mass: gx = g sin θ, gy = −g cos θ. No slip at y = 0; at the free surface **μ du/dy = 0 at y = H** (no shear) and p = p_a.

(b) y-momentum: ∂p/∂y = −ρg cos θ ⇒ **p = p_a + ρg cos θ (H − y)**. x-momentum: μu″ = −ρg sin θ (the pressure does not vary along x). With u(0) = 0 and u′(H) = 0: **u = (ρg sin θ/μ)(Hy − y²/2)**, q = ρg sin θ H³/(3μ), so **H = (3μq/(ρg sin θ))^{1/3}**.

(c) H = (3 × 0.01 × 10⁻⁵/(1050 × 9.80665 × sin 17.5°))^{1/3} = **0.459 mm**; surface speed u(H) = ρg sin θ H²/(2μ) = **0.0327 m/s** (1.5 times the mean speed 0.0218 m/s, as for any half-parabola); Re = ρq/μ = **1.05**: firmly laminar, and slow enough that the film surface stays smooth.

(d) `ns_residuals(u, 0, rho*g*cos(theta)*(H - y), rho, mu, gx=g*sin(theta), gy=-g*cos(theta))` returns **(0, 0, 0)**, and `diff(u, y).subs(y, H)` is 0 ✓.

**Independent check.** The film is half of a plane Poiseuille flow of height 2H driven by the "pressure gradient" ρg sin θ: the free surface is its symmetry plane. Its flow rate per unit width is half of G(2H)³/(12μ) with G = ρg sin θ, i.e. ρg sin θ H³/(3μ) ✓.

**Common mistakes.** Imposing u = 0 at the free surface; forgetting the cube root's sensitivity: doubling the flow rate thickens the film by only 26 %.

## 4. Your own actuator (s = 15: h = 7 µm)

**Assumptions.** As in the lesson: concentric piston, plane slot of width πD, typical fluid properties, rigid walls, constant viscosity.

(a)

| piston speed | leakage (L/min) | Re_gap | shear force on piston (N) | t₉₉ (µs) |
|---|---|---|---|---|
| 0 | 0.01964 | 0.298 | 5.690 | 1.644 |
| +0.1 m/s (with the leak) | 0.02129 | 0.323 | 5.531 | 1.644 |
| −0.1 m/s (against it) | 0.01799 | 0.273 | 5.850 | 1.644 |

The pressure-driven part, G(h/2) × πDL = **5.690 N**, is the same at every speed; only the Couette part, ∓μU/h × πDL = ∓**0.159 N**, changes. Either way the force is a few newtons against an actuator force of order 10 kN: negligible.

(b) q = 0 when U = −Gh²/(6μ) = **−1.19 m/s**, a piston moving against the leak at more than ten times typical actuator rates. Not reachable; the leak always wins.

(c) h = 6 µm: factor **0.630**; h = 8 µm: factor **1.493**. A ±14 % tolerance on h is a −37 %/+49 % spread in leakage (the h³ law). A leakage specification must come with a tight clearance tolerance, or be verified by test.

(d) `crank_nicolson` across the 7 µm gap from rest to t₉₉ = 1.644 µs gives **0.98998** of the steady centreline speed, against 0.99000 from the full series: the one-term formula for t₉₉ is accurate because the n = 3 term has decayed by e^{−8π²νt/h²} more than the first.

**Common mistakes.** Using the diameter instead of the circumference for the slot width; reporting leakage per unit width.

## 5. Stokes's second problem

(a) With u = Re{U e^{iωt − (1+i)ky}}: u_t = iω u and u_yy = (1 + i)²k² u = 2ik² u. So u_t = νu_yy requires iω = 2iνk², **k = √(ω/(2ν))** ✓ (or `ns_residuals` returns zeros). The amplitude e^{−ky} falls to 1 % at **y = ln(100)/k = 4.605/k**; for ν = 1, ω = 2π, k = 1.7725 and that depth is 2.598.

(b) Errors after three periods 1.580 × 10⁻³, 3.995 × 10⁻⁴, 1.001 × 10⁻⁴ for n = 40, 80, 160: observed orders **1.98 and 2.00**.

(c) From rest, the maximum difference from the periodic solution is 6.2 × 10⁻³ after one period, 1.6 × 10⁻³ after two, **7.6 × 10⁻⁴ after three** and 3.4 × 10⁻⁴ after five: three periods. The transient is the start-up wave diffusing away; it decays slowly (algebraically), not exponentially.

(d) k = √(πf/ν) with ν = 14 mm²/s: 1/k = **668 µm at 10 Hz** and **211 µm at 100 Hz**, 20–70 times the clearance. Inside a 10 µm gap the oscillating shear layer is far thicker than the gap, so the gap flow follows the wall motion instantaneously: quasi-steady Couette flow, as assumed.

**Common mistakes.** Imposing u = 0 at a far boundary that is too close (the e^{−6} = 2.5 × 10⁻³ residual at y = 6/k shows up as a floor in the error); expecting second-order convergence from an impulsive start, where the initial data are not smooth.

## 6. When the exact solution stops describing the actuator

(a) ΔT = Δp/(ρc_p) = 20.7 × 10⁶/(845 × 1900) = **12.9 K**. The leaking fluid warms along the land, its viscosity falls, and the leakage, proportional to 1/μ, rises. A real calculation integrates the temperature along the gap (a thermal Couette–Poiseuille problem).

(b) **2.5**: eccentricity alone can multiply the leakage by up to 2.5, more than any other single effect here. Real pistons are pushed to one side by side loads, so the concentric figure is a lower bound.

(c) With Re_gap = 0.87 and h = 10 µm, 0.06 Re h = **0.52 µm**, against a 12 mm land: the flow is fully developed after one twenty-thousandth of the land. Justified.

(d) Ranked by effect on leakage: eccentricity (up to ×2.5), clearance tolerance (×0.63 to ×1.49 for ±1 µm on 7 µm; ×0.73 to ×1.33 for ±1 µm on 10 µm), viscous heating (a 12.9 K rise, whose effect depends on the fluid's viscosity–temperature curve), entrance length (negligible). Measure the leakage on the assembled actuator at operating temperature first: it captures all four at once. The model's value is in telling you which variable to control: h, cubed.

**Common mistakes.** Treating the leakage as known to three figures because the formula is exact; an exact solution of an idealised problem is not an exact prediction.

## Check-yourself answers

1. ρ(∂u/∂t + u·∇u) = −∇p + μ∇²u + ρg: unsteady and convective acceleration (mass × acceleration per volume), pressure gradient, viscous force, body force; every term in N/m³ (Pa/m).
2. A Newtonian fluid has a stress that is linear in the rate of strain, τ = −pI + 2μS, with one constant viscosity. Incompressibility (∇·u = 0) removes the bulk-viscosity term and turns ∇·(2μS) into μ∇²u for constant μ.
3. The central difference's truncation error is (Δy²/12)u'''' + …, and u'''' = 0 for a parabola; the scheme is exact for polynomials up to degree three, so the discrete solution equals the exact one at the nodes on any grid.
4. Slot leakage per unit width ∝ h³ (Gh³/(12μ)); pipe flow ∝ R⁴ (πR⁴G/(8μ)): halving the gap divides the leakage by 8, halving the radius divides the flow by 16.
5. t ~ h²/ν (the start-up of plane Poiseuille flow reaches 99 % at 0.47 h²/ν). If the forcing changes slowly compared with h²/ν, the flow is quasi-steady; for a 10 µm hydraulic gap h²/ν is 7 µs.
6. From the momentum balance of a slab of fluid: the net pressure force G h per unit length must be carried by the difference between the shear stresses on the two walls.
7. The problem has no length scale of its own, so u/U can depend on y and t only through y/√(νt). Its 99 % thickness is 3.64√(νt).
8. Any two of: the piston is not concentric (eccentricity raises leakage, up to ×2.5); the fluid warms as it leaks (lower μ, more leakage); the bore deforms under 20 MPa (a larger gap, more leakage); the clearance varies along the land and with tolerance (h³).
