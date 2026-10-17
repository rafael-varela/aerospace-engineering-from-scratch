# P08.09  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe script that calls it, given to the figures shown.

## 1. A Lamb–Oseen vortex (s = 15: Γ = 25 m²/s)

**Assumptions.** An axisymmetric Lamb–Oseen vortex at one instant, r_c² = 0.25 m².

(a) r_max = 1.1209 × 0.5 = **0.5605 m**; u_max = (25/(2π × 0.5605))(1 − e^{−1.2564}) = **5.078 m/s**.

(b) Γ(r) = Γ(1 − e^{−r²/r_c²}): within 0.5 m **15.803 m²/s (63.21 %)**, within 1 m **24.542 m²/s (98.17 %)**, within 2 m **25.000 m²/s (100.00 %)**. `circulation_polygon` on 2000-point circles gives 15.8030, 24.5421, 25.0000 ✓. Within one core radius the vortex holds 1 − 1/e of its circulation.

(c) ω(0) = Γ/(π r_c²) = **31.83 s⁻¹**. Near the axis 1 − e^{−r²/r_c²} ≈ r²/r_c², so u_θ ≈ Γr/(2π r_c²) = (ω(0)/2) r: solid-body rotation at **15.92 rad/s**, half the vorticity, as for any solid-body rotation (whose vorticity is twice its angular velocity).

(d) Off-centre disk: line integral **4.386838**, area integral **4.386845** m²/s: Stokes's theorem to 1.6 × 10⁻⁶, the midpoint rule's O(Δs²) error on 2000 chords.

**Common mistakes.** Taking r_c as the radius of peak swirl; forgetting that circulation counts all the vorticity inside, so it saturates outside the core.

## 2. Derivations

(a) For a material loop, d/dt ∮ u·dl = ∮ (Du/Dt)·dl + ∮ u·d(dl/dt), and dl/dt = du for a line element moving with the fluid. Euler with a barotropic fluid (∇p/ρ = ∇P with P(p) = ∫dp/ρ) and a conservative body force g = −∇Φ: Du/Dt = −∇(P + Φ). So dΓ/dt = −∮ d(P + Φ) + ∮ d(|u|²/2) = **0**, since the integral of an exact differential round a closed loop vanishes. The three conditions (inviscid, barotropic, conservative force) are each needed: a baroclinic fluid (density varying on pressure surfaces, as in a sea breeze) or a non-conservative force (Coriolis in a rotating frame is one) generates circulation. A viscous fluid adds **ν∮∇²u·dl = −ν∮(∇×ω)·dl**, which is zero only where the vorticity field is uniform; a vortex core with a Gaussian ω has ∇×ω ≠ 0 wherever ω varies, which is where the loop of Problem 4 cuts the core.

(b) With ω = ω(r, t) and u = u_θ(r)e_θ, (u·∇)ω = (u_θ/r)∂ω/∂θ = 0: the convective term vanishes, and the vorticity equation is ω_t = ν(ω_rr + ω_r/r). For ω = Γ/(πR) e^{−r²/R} with R = r_c0² + 4νt: ω_t = ω(−1/R + r²/R²) Ṙ = 4νω(r²/R² − 1/R); ω_r = −2rω/R; ω_rr = ω(4r²/R² − 2/R); ν(ω_rr + ω_r/r) = νω(4r²/R² − 4/R) = **ω_t** ✓. Then Γ(r) = ∫₀^r 2πr′ω dr′ = **Γ(1 − e^{−r²/R})** and **u_θ = Γ(r)/(2πr)**.

(c) Put the filament on the x axis from A = (x_A, 0, 0) to B = (x_B, 0, 0), the field point at (0, d, 0). The Biot–Savart segment formula gives the speed (Γ/(4πd))(cos θ_A − cos θ_B), with θ_A, θ_B the angles between the filament and the lines to its ends. Infinite line: θ_A → 0, θ_B → π, **Γ/(2πd)**. Semi-infinite, the point abeam its end: θ_A = π/2, θ_B → π, **Γ/(4πd)**. The horseshoe at the centre of its bound vortex: the bound vortex induces nothing on its own line, and each trailing leg is semi-infinite with the point abeam its end at distance b′/2, giving 2 × Γ/(4π b′/2) = **−Γ/(πb′)** (downward, by the right-hand rule for legs trailing behind a lifting wing). Far downstream, each leg looks infinite: **−2Γ/(πb′)**, twice the value at the wing.

**Common mistakes.** Forgetting the term from the moving loop in (a) (it integrates to zero, but it is there); dropping ω_r/r in polar coordinates.

## 3. Crossing an airliner's wake

**Assumptions.** Elliptic loading (b′ = πb/4), straight line vortices, sea-level air, laminar diffusion only in (c).

(a) b′ = π × 34/4 = 26.70 m; Γ₀ = 70 000 × 9.80665/(1.225 × 75 × 26.70) = **279.8 m²/s**. Descent rate Γ₀/(2πb′) = **1.67 m/s**: the pair sinks about 100 m a minute.

(b) At 5 m from the axis the vortex induces Γ₀/(2π × 5) = **8.91 m/s**, up at one tip and down at the other: a change in angle of attack of **±8.4°** at 60 m/s. The light aeroplane's ailerons cannot oppose a rolling moment that large: the wake rolls it faster than its controls can stop. This is why light aircraft are kept clear of airliners' wakes.

(c) r_c² = 4νt with ν = 1.46 × 10⁻⁵ m²/s: after 60 s the peak swirl is at **6.6 cm** and its speed is **480 m/s**; after 120 s, at **9.4 cm** and **339 m/s**. Halving the circulation within 5 m needs e^{−25/(4νt)} = ½, t = 25/(4ν ln 2) = **6.2 × 10⁵ s, 7.1 days**. The numbers do not make physical sense: a 480 m/s swirl in a 7 cm core is supersonic, and no wake lasts a week. Laminar diffusion from a zero core is the wrong model for a wake.

(d) A real wake vortex rolls up from the whole span's shed vorticity into a core of order metres, is turbulent from the start (an eddy viscosity many orders of magnitude above ν), and is destroyed in minutes by atmospheric turbulence, by the mutual-induction (Crow) instability that links the pair into rings, and by the ground. The Lamb–Oseen solution is exact, but its assumptions are not met: laminar, axisymmetric, isolated, and starting from a point. It stays useful as the profile *shape* fitted to a measured core.

**Common mistakes.** Using b instead of b′; treating the induced angle as small.

## 4. Kelvin's theorem for your own loop (s = 15: centre (0.8, 0))

(a) u_θ(0.8) = 0.09404 m/s (Γ = 1, r_c = 1), turnover **53.45** (time units of the problem).

(b) 4000 points:

| t (turnovers) | ν | line integral | area integral | length |
|---|---|---|---|---|
| 0 | 0 | 0.17733267 | 0.17733274 | 3.770 |
| 1 | 0 | 0.17733252 | 0.17733274 | 8.708 |
| 2 | 0 | 0.17733207 | 0.17733274 | 16.126 |
| 4 | 0 | 0.17733028 | 0.17733274 | 31.471 |
| 0 | 0.005 | 0.17733267 | 0.17733274 | 3.770 |
| 1 | 0.005 | 0.12030067 | 0.12030077 | 5.807 |
| 2 | 0.005 | 0.08940770 | 0.08940783 | 7.305 |
| 4 | 0.005 | 0.05866224 | 0.05866238 | 8.703 |

Inviscid: constant (to 1.4 × 10⁻⁵ relative after four turnovers, while the loop stretches eightfold). Viscous: falling by a factor of three, and the line and area integrals agree to better than 2 × 10⁻⁶ relative at every time: Stokes's theorem holds on the deformed material loop, and the circulation changes because the vorticity inside it diffuses out.

(c) Inviscid change after four turnovers: −3.829 × 10⁻⁵, −9.568 × 10⁻⁶, −2.392 × 10⁻⁶ for 1000, 2000, 4000 points: ratios 4.00 and 4.00, **order 2**. The drift is the midpoint rule's O(Δs²) error on chords that lengthen as the loop stretches; it vanishes as the points are refined, so it is not a change of the true circulation.

(d) Viscosity spreads the vortex: the swirl profile flattens, the differential rotation between neighbouring circles falls, and the loop is sheared less. The inviscid vortex keeps its sharp profile and keeps shearing the loop at the same rate.

**Common mistakes.** Evaluating the velocity with the initial core size at later times; reading the inviscid drift as physics.

## 5. Downwash at your tailplane (s = 15: l_t = 5.0 m)

(a) b′ = 7.854 m, Γ₀ = 7355.0/(0.90913 × 60 × 7.854) = **17.168 m²/s**. At (5.0, 0, 0): w = **−1.5805 m/s**, ε = **1.509°**; at (5.0, 0, 0.5): w = **−1.5547 m/s**, ε = **1.484°**.

(b) In the wing plane the bound vortex contributes **−0.3375 m/s** and the trailing legs **−1.2430 m/s**. The far-wake value is −2Γ₀/(πb′) = **−1.3916 m/s** (1.329°). The tail sees more downwash than the far wake because, at only 5 m behind, the bound vortex still adds its share while the legs have not reached their full (infinite-line) value.

(c) At z = 0.5 m: −1.5547, −1.6385 and −1.9607 m/s at y = 0, 1 and 2 m. Moving outboard brings the point closer to the trailing leg on that side, whose induced velocity grows as 1/distance; the horseshoe concentrates all the trailing vorticity at ±b′/2, which exaggerates the growth compared with a real, distributed wake.

(d) Downwash is proportional to Γ, so a 30 % rise raises ε from 1.48° to about 1.93°, reducing the tail's angle of attack by **0.45°**. Lowering flaps changes the trim, and the pilot (or the trim system) must correct the tail's lift.

**Common mistakes.** Forgetting the bound vortex's contribution at the tail; using the 2-D strip circulation (13.48 m²/s) in a horseshoe model, which needs the centre-span value.

## 6. How fast does lift build up?

**Assumptions.** Plane flow; the bound vortex at the wing and a point starting vortex of the opposite sign at distance d behind; Γ = πcV sin α = 13.484 m²/s (c = 1.5 m, V = 60 m/s, α = 2.733°).

(a) w = Γ/(2πd), so w/(Vα) ≈ c/(2d): **50 %, 10 %, 5 %, 1 %** of α at d = 1, 5, 10, 50 chords (w = 1.43, 0.286, 0.143, 0.0286 m/s).

(b) **5 chords**: 7.5 m of travel, an eighth of a second at 60 m/s.

(c) The model gets the mechanism right: Kelvin's theorem forces an equal and opposite starting vortex, which induces a downwash at the wing that reduces its effective angle of attack until the wing has moved several chords away. It gets the scale right too: the influence fades over a few chords. It is not accurate in detail, because the real wake is a continuous sheet shed while Γ grows, and its nearest part, the part with the largest effect, is the most recently shed. For a gust-response calculation I would use the Wagner (or, for gusts, Küssner) indicial function and check the time scale against c/V, not the point-vortex fraction.

**Common mistakes.** Using a 3-D horseshoe circulation in a 2-D model; forgetting the small-angle form when dividing by α.

## Check-yourself answers

1. ω = ∇×u (in plane flow ω = ∂v/∂x − ∂u/∂y), twice the local angular velocity of a fluid particle. Γ = ∮ u·dl round a closed curve. Stokes's theorem: Γ = ∬ ω·n dA over any surface bounded by the curve.
2. The circulation round a material loop is constant in time if the fluid is inviscid, barotropic (p a function of ρ only) and acted on by conservative body forces.
3. Before the wing moves, a large loop round it has zero circulation. Kelvin keeps it zero as the loop moves with the air. The wing's bound circulation must therefore be balanced by an equal and opposite vortex shed from the trailing edge into the wake and left behind.
4. Zero everywhere except at the centre (the flow is irrotational outside); the circulation round any loop that encloses it is Γ, and round any loop that does not is zero.
5. r_c² = r_c0² + 4νt: the core grows as √(νt). The total circulation stays Γ; far outside the core the flow is still the point vortex Γ/(2πr).
6. dV = (Γ/4π) dl × r/|r|³ integrated along the filament; for a straight segment, (Γ/(4πd))(cos θ_A − cos θ_B) normal to the plane of filament and point; for an infinite line, Γ/(2πd).
7. At the wing each trailing leg is semi-infinite with the point abeam its end, Γ/(4π · b′/2) each; far behind they look infinite, Γ/(2π · b′/2) each: twice as much.
8. The polygon's midpoint rule has an O(Δs²) error, and the loop's segments stretch as it is sheared. Halving Δs (doubling the points) divided the drift by 4.00, so the drift is discretisation error that tends to zero, while the exact circulation (the area integral) is constant.
