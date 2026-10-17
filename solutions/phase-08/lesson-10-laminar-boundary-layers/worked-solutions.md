# P08.10  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe script that calls it, given to the figures shown. The Blasius constants used are the code's: f″(0) = 0.332057, δ99 = 4.910, δ* = 1.72079, θ = 0.66411, H = 2.5911 (in units of x/√Re_x).

## 1. A wind-tunnel fin (s = 15: U = 27.5 m/s)

**Assumptions.** Flat plate at zero incidence, zero pressure gradient, laminar from the leading edge, steady, incompressible (M ≈ 0.08), sea-level ISA air.

(a) Re_L = 27.5 × 0.20/1.4607 × 10⁻⁵ = **3.765 × 10⁵**. With x/√Re_L = 0.20/613.6 = 3.259 × 10⁻⁴ m: δ99 = **1.600 mm**, δ* = **0.561 mm**, θ = **0.2165 mm**.

(b) c_f = 0.664/√Re_x: at mid-chord (Re_x = 1.883 × 10⁵) **c_f = 1.531 × 10⁻³, τ_w = 0.709 Pa**; at the trailing edge **c_f = 1.082 × 10⁻³, τ_w = 0.501 Pa** (q = 463.2 Pa). The ratio is **√2 = 1.414**: c_f ∝ x^(−1/2) because the wall gradient U/δ falls as δ ∝ √x grows.

(c) C_F = 1.328/√Re_L = **2.165 × 10⁻³**; D = q C_F × (2 × 0.20 × 0.6 m²) = **0.241 N**.

(d) Re_L = 3.8 × 10⁵ is below the range where flat-plate transition is usually seen in a quiet tunnel, so laminar flow over the whole chord is plausible but not guaranteed: tunnel turbulence, a blunt leading edge or roughness can trip it earlier. `P08.12` computes where small disturbances start to grow and where the e^N method puts transition.

**Check.** δ*/θ = 0.561/0.2165 = 2.59 = H. Integrating τ_w along the chord (the test `test_integrated_drag_equals_integral_of_local_shear` does this) gives the same drag to 10⁻⁶.

**Common mistakes.** Forgetting the second side; using the trailing-edge c_f as the average (it is half of C_F).

## 2. Derivations

(a) With ψ = √(νUx) f(η), η = y√(U/(νx)): η_y = √(U/(νx)), η_x = −η/(2x).
u = ψ_y = √(νUx) f′ η_y = **U f′(η)**.
v = −ψ_x = −[½√(νU/x) f + √(νUx) f′ (−η/(2x))] = **½√(νU/x)(η f′ − f)**.
u_x = U f″ η_x = −U η f″/(2x); u_y = U f″ √(U/(νx)); u_yy = U f‴ U/(νx).
Then u u_x + v u_y = −U²ηf′f″/(2x) + ½(U²/x)(ηf′ − f)f″ = −U² f f″/(2x), and ν u_yy = U² f‴/x. Equating and dividing by U²/x: **f‴ + ½ f f″ = 0**. Assumptions: the boundary-layer approximations (thin layer, ∂²/∂x² ≪ ∂²/∂y², pressure impressed by the outer flow), zero pressure gradient (U constant, so no U dU/dx term), steady, incompressible, Newtonian. Conditions: no slip u = 0 ⇒ f′(0) = 0; impermeable wall v = 0 ⇒ f(0) = 0; u → U outside ⇒ f′(∞) = 1.

(b) If f = aF(aη) then f′ = a²F′, f″ = a³F″, f‴ = a⁴F‴ and f f″ = a⁴FF″; so f‴ + ½ff″ = a⁴(F‴ + ½FF″) = 0. The wall conditions survive (f(0) = aF(0) = 0, f′(0) = 0), f″(0) = a³F″(0) = a³ and f′(∞) = a²F′(∞). Choosing a = [F′(∞)]^(−1/2) gives f′(∞) = 1 and **f″(0) = a³ = [F′(∞)]^(−3/2)**. The Blasius problem has a scaling symmetry, so the "unknown" initial condition can be found from one integration.

(c) Integrate u u_x + v u_y = ν u_yy from y = 0 to h (above the layer). Integrating v u_y by parts with v(0) = 0 and v_y = −u_x: ∫v u_y dy = v(h)U + ∫u u_x dy, and v(h) = −∫₀ʰ u_x dy. So ∫(2u u_x − U u_x) dy = −ν u_y(0). Since U is constant, 2u u_x − U u_x = −∂/∂x[u(U − u)], giving
**d/dx ∫₀^∞ u(U − u) dy = ν (∂u/∂y)₀ = τ_w/ρ**, i.e. **dθ/dx = c_f/2**.
With θ = C√(νx/U) and c_f = 2f″(0)/√Re_x: C/(2√Re_x) · ... precisely dθ/dx = (C/2)√(ν/(Ux)) = C/(2√Re_x), and c_f/2 = f″(0)/√Re_x, so **C = θ√Re_x/x = 2f″(0)**. The code's θ integral, 0.66411, equals 2 × 0.332057 = 0.664115 to the trapezoid rule's accuracy (the test asserts 10⁻⁴).

**Common mistakes.** Dropping the η f′ term in v; treating U as varying (it does not on a flat plate).

## 3. Why tunnel test sections diverge (s = 15: U = 18 m/s)

(a) Re_x = 18 × 2/1.4607 × 10⁻⁵ = 2.465 × 10⁶ (laminar is an assumption here; see (d)). x/√Re_x = 1.274 × 10⁻³ m: **δ99 = 6.26 mm, δ* = 2.19 mm**.

(b) The core area is (0.5 − 2 × 0.00219)² against 0.25 m²: **U_exit/U_in = 1.0178, a 1.78 % rise**. The linear approximation 4δ*/W = 1.75 % is close because δ* ≪ W. A 1.8 % velocity change is a 3.6 % change in dynamic pressure along the model: too much for force testing.

(c) Each wall must move outward as fast as δ* grows: dδ*/dx = δ*/(2x) = **5.48 × 10⁻⁴ rad = 0.031°**. The outward velocity at the edge is v/U = 0.8604/√Re_x = **5.48 × 10⁻⁴**, the same number: the flow just outside the layer is deflected outward by exactly the slope of δ*. That is what "displacement thickness" means.

(d) Larger: a turbulent layer's δ* grows roughly as x^0.8 rather than x^0.5 and is several times thicker at the same Re_x, so real tunnels diverge more, and their exact divergence is set by calibration.

**Common mistakes.** Using δ99 instead of δ* for the blockage; counting only two walls.

## 4. Shooting accuracy (s = 15: η_max = 9)

(a) Reference `blasius_bvp(9, tol=1e-10)`: 0.3320573908. Errors: n = 90: **8.45 × 10⁻⁹**; 180: **5.56 × 10⁻¹⁰**; 360: **3.55 × 10⁻¹¹**; 720: **2.18 × 10⁻¹²**; observed orders **3.93, 3.97, 4.03**: RK4 is fourth order. The reference must use the same η_max because the truncation error of the domain (part b) is about 5 × 10⁻⁸ at η_max = 9, larger than the RK4 error you are trying to measure.

(b) f″(0) − 0.332057336: η_max = 4: **+2.07 × 10⁻²**; 5: **+4.10 × 10⁻³**; 6: **+5.09 × 10⁻⁴**; 7: **+3.90 × 10⁻⁵**; 8: **+1.86 × 10⁻⁶**; 9: **+5.5 × 10⁻⁸**; 10: **+1.2 × 10⁻⁹**; 12: 2 × 10⁻¹⁰. **η_max ≥ 9** gives six correct digits. Forcing f′ = 1 at a finite η_max squeezes the layer: the profile must reach the free stream sooner than it naturally would, so its wall slope is steeper. The error falls roughly like exp(−(η − 1.72)²/4), the decay of f″ outside the layer.

(c) One integration with F″(0) = 1 to ζ = 10 (2000 steps) gives F′(10) = 2.0854092 and **f″(0) = 0.3320573362**, correct to about ten digits, because in the scaled variable the layer reaches only to aζ ≈ 6.9 and F″ has decayed much further at ζ = 10. The cost is one integration against eight for the secant shooting started from (0.2, 0.5).

**Common mistakes.** Measuring the RK4 order against the published value (the domain error dominates and the order looks like zero); reading "more digits than the domain allows".

## 5. The PDE forgets its starting profile (s = 15: U = 10 m/s)

(a) c_f√Re_x at x = 1 m: nx = 200: **0.65673** (error 7.38 × 10⁻³); 400: **0.66047** (3.65 × 10⁻³); 800: **0.66222** (1.90 × 10⁻³); 1600: **0.66306** (1.05 × 10⁻³). Observed orders **1.02, 0.94, 0.85**: first order, so the x-step (backward Euler with lagged coefficients) dominates; the y-error (second order, 100 cells across δ99 at x = L) is smaller, and its presence is why the order drifts below 1 on the finest grids. At nx = 1600 the error for ny = 75, 150, 300, 600 is 2.39, 1.50, 1.05 and 0.88 × 10⁻³: going from 300 to 600 changes it by only 1.7 × 10⁻⁴, so at ny = 300 the y-contribution is about 2 × 10⁻⁴ and the x-step owns the rest.

(b) Richardson with p = 1: 2 × 0.66306 − 0.66222 = **0.66391**, within 0.03 % of 2f″(0) = 0.66411.

(c) Largest |u/U − f′(η)|: **1.46 × 10⁻³ at x = 0.25 m, 1.13 × 10⁻³ at x = 1 m**. A uniform start is nothing like the Blasius profile, but after a short distance the layer depends only on the local Re_x: it has no memory of its initial condition, which is why one similarity solution describes every station. (The same is true for the impulsive start of Stokes's first problem in `P08.06`, which is the time analogue.)

**Common mistakes.** Measuring c_f with a one-sided first-order wall derivative (that adds an O(Δy) error larger than the effect you are measuring); refining nx and ny together and calling the order "1.5".

## 6. Where Blasius stops describing the wing

(a) δ99/x = 4.91/√Re_x = 0.1 ⇒ **Re_x = 2.41 × 10³**, x = Re_x ν/U = **0.75 mm** from the leading edge at 60 m/s and 3000 m. Upstream of that the "thin layer" approximation is poor and the leading-edge shape matters; the region is too short to change the drag measurably.

(b) (1) **Pressure gradient** (an aerofoil accelerates the flow to a suction peak and then decelerates it): accelerating flow thins the layer and raises c_f near the nose; decelerating flow thickens it and lowers c_f, ending in laminar separation (`P08.11`). (2) **Transition** to turbulence somewhere along the chord (`P08.12`): turbulent friction is several times larger, so the drag rises well above 29.7 N (`P08.14` estimates by how much). (3) **Surface roughness, insects and rain**: they trip transition early and raise drag (`P08.12`, `P08.14`). Also: the wing's thickness and three-dimensionality (sweep, tip flow), which the flat plate ignores.

(c) A row of flush hot-film sensors or a surface temperature camera (infrared thermography) along the chord: a laminar layer gives a low, smoothly falling heat-transfer or shear signal ∝ x^(−1/2); transition shows as a sudden rise of several times. On a tunnel model, an oil-flow or china-clay pattern shows the same line. If the laminar model held, the signal would fall smoothly to the trailing edge with no jump.

**Common mistakes.** Treating 0.75 mm as important for drag; claiming that the pressure gradient changes nothing because "Blasius is exact" (it is exact only for a flat plate).

## Check-yourself answers

1. δ99: the height at which u = 0.99U (4.91 x/√Re_x). δ*: the distance the wall would have to move to carry the same mass flux with uniform flow, ∫(1 − u/U)dy (1.7208 x/√Re_x). θ: the same for the momentum flux, ∫(u/U)(1 − u/U)dy (0.6641 x/√Re_x). H = δ*/θ = 2.591.
2. δ ∝ √(νx/U): it grows as √x and shrinks as 1/√U. Quadrupling the speed halves δ at the trailing edge.
3. c_f = 0.664/√Re_x falls because the layer thickens, so the velocity gradient at the wall weakens. Averaging x^(−1/2) over 0..L gives 2 × the value at L: C_F = 1.328/√Re_L.
4. f″(0) is the dimensionless wall shear; the problem fixes f, f′ at the wall and f′ at infinity, so the missing wall condition f″(0) must be chosen (shot) until the far condition is met.
5. Imposing f′ = 1 at a finite η_max forces the profile to reach the free stream early, steepening it at the wall; the error decays rapidly with η_max (2 × 10⁻² at 4, 2 × 10⁻⁶ at 8).
6. Steady, incompressible, Newtonian, two-dimensional, laminar flow with a layer thin compared with x (Re_x ≫ 1), pressure set by the outer flow (here constant). It fails near the leading edge (Re_x of order 10³ and below), at separation, and where the flow becomes turbulent.
7. The layer slows the fluid, so continuity pushes fluid outward; at the edge v = 0.8604 U/√Re_x, which equals U dδ*/dx.
8. Mass (continuity) and x-momentum in differential form (Prandtl's simplification of Navier–Stokes), with the Newtonian laminar stress τ = μ ∂u/∂y; the momentum integral (`dθ/dx = c_f/2`) is the same law in integral form across the layer.
