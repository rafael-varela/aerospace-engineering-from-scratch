# P10.09  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and a probe script that calls it, given to the figures shown.

## 1. The downwash at the tail

**Assumptions.** One horseshoe with legs πb/4 = 7.854 m apart, flat wake, the tail on the bound vortex's plane.

(a) Γ = W/(ρVb′) = 7357.5/(1.225 × 55 × 7.854) = **13.90 m²/s** (the same as the elliptic Γ0, since πb/4 × Γ0 is the elliptic load's lift integral). Centre downwash Γ/(πb′) = **0.564 m/s** (0.587°).

(b) At x = 4.5 m: **1.311 m/s (1.37°)**; far downstream: **1.127 m/s (1.17°)**.

(c) Just behind the wing the bound vortex itself still adds downwash (behind a lifting vortex the flow moves down), and the legs have nearly their full strength; far downstream the bound vortex's effect has died away. The far-field angle of the horseshoe (1.17°) is below the elliptic wing's 2α_i = 1.45°: the horseshoe keeps the lift and the centroid of the shed vorticity, but concentrating the sheet into two lines changes the velocity on the centreline, to 2Γ0/(πb′) = 8Γ0/(π²b), which is 8/π² = 0.81 of the elliptic sheet's Γ0/b.

**Common mistakes.** Using b instead of πb/4 for the legs; forgetting the bound vortex at the tail.

## 2. Derivations

(a) Put a point of the segment at A + t r0, t ∈ [0, 1]; the vector from it to P is r = r1 − t r0 and dl = r0 dt. Then dl × r = r0 × r1 dt, a fixed direction, so V = (Γ/4π)(r0 × r1)∫₀¹ dt/|r1 − t r0|³. Measure distance along the line by σ = t|r0| from the foot of the perpendicular; with h the distance of P from the line, |r|² = h² + σ² and ∫ dσ/(h² + σ²)^{3/2} = σ/(h²√(h² + σ²)), which between the ends gives (cos θ₁ − cos θ₂)/h², where θ₁, θ₂ are the angles between r0 and r1, r2. Now |r0 × r1| = h|r0|, r0 × r1 = r1 × r2 (because r2 = r1 − r0), |r1 × r2| = h|r0|, and |r0| cos θ_k = r0·r_k/|r_k|. Substituting, **V = (Γ/4π)(r1 × r2)/|r1 × r2|² · r0·(r1/|r1| − r2/|r2|)**. Symmetric case: cos θ₁ − cos θ₂ = 2L/√(L² + h²), giving **V = (Γ/2πh) L/√(L² + h²)**.

(b) With B = A + Ld: r0 = Ld, r2 = r1 − Ld, r1 × r2 = Ld × r1, r2/|r2| → −d. Then V → (Γ/4π)(Ld × r1)/(L²|d × r1|²) · L(d·r1/|r1| + 1) = **(Γ/4π)(d × r1)/|d × r1|² (1 + d·r1/|r1|)**. In the plane through A, d·r1 = 0 and |d × r1| = h: V = Γ/(4πh), half the infinite line.

(c) Γ = Γ0 sin φ, y = −(b/2) cos φ, so dΓ/dy = (dΓ/dφ)/(dy/dφ) = Γ0 cos φ/((b/2) sin φ). Then (3) of the lesson is w(θ) = (1/4π)∫₀^π (dΓ/dφ) dφ/(y(θ) − y(φ)) with y(θ) − y(φ) = (b/2)(cos φ − cos θ): w = (Γ0/(2πb))∫₀^π cos φ dφ/(cos φ − cos θ) = (Γ0/(2πb)) × π = **Γ0/(2b)**, independent of θ.

(d) CD/CL = CD0/CL + kCL with k = 1/(πeAR); zero derivative at CD0 = kCL²: **CL* = √(CD0/k) = √(πeAR CD0)**, and CL/CD = CL*/(2CD0) = **½√(πeAR/CD0)**. CD/CL^{3/2} = CD0 CL^{−3/2} + kCL^{1/2}; derivative zero at 3CD0 = kCL²: **CL_mp = √3 CL***.

**Common mistakes.** Losing the sign of the right-hand rule; integrating (3) without the principal value.

## 3. A wake encounter

**Assumptions.** Elliptic load, fully rolled-up point vortices of strength Γ0 at ±πb/8, infinite lines, no decay.

(a) Γ0 = 4W/(ρVπb) = 4 × 686 700/(1.225 × 70 × π × 60) = **170 m²/s**; separation πb/4 = **47.1 m**.

(b) Midway each vortex is 23.6 m away and both push down: w = 2 × 170/(2π × 23.6) = **2.30 m/s**.

(c) At 5 m from one core: v = 170/(2π × 5) = **5.41 m/s**, about 10 % of the light aircraft's 55 m/s, and ten times the downwash its own wing makes (0.56 m/s at its centre). Across the light aircraft's 10 m span the vortex field reverses (up on one side, down on the other), giving a rolling moment the ailerons may not be able to hold. Hence the regulated separations behind heavy aircraft, and the advice to stay above an airliner's flight path, since wakes sink.

**Common mistakes.** Using the vortex strength of a single horseshoe per wing *half* (it is Γ0 for the whole tip vortex); forgetting that the airliner's wake sinks and decays (which helps).

## 4. A non-elliptic load (s = 15: a₃ = −0.15)

(a) e = **0.93725, 0.93689, 0.93680** on 20, 40, 80 strips, against **1/(1 + 3 × 0.0225) = 0.93677**.

(b) Near/far = **1.0000000** (to 10⁻¹⁵): for a flat wake the legs at the bound line are exactly half-infinite, so w = w_T/2 identically and (4) and (5) of the lesson are the same sum.

(c) Γ = sin φ − 0.15 sin 3φ is, near the tips (small φ), φ − 0.45φ = 0.55φ: **less** load near the tips than elliptic, more at mid-semispan (it is a step towards the bell of `P10.12`).

**Common mistakes.** Using a₃ > 0 and a₃ < 0 interchangeably: e depends on a₃², the shape does not.

## 5. Your own polar (s = 15: CD0 = 0.040, e = 0.825)

(a) (L/D)max = **11.62** at CL* = **0.929**, V_md = **32.8 m/s**; CL_mp = **1.610**, V_mp = **24.9 m/s**.

(b) At 55 m/s: D = **1002 N**, induced share **11.2 %**, P = **55.1 kW**.

(c) V_mp = 24.9 m/s is *below* the clean stall speed (26.4 m/s at CL_max 1.44). The aircraft cannot fly there clean: its practical best-endurance speed is a little above the stall, where power is close to the minimum anyway (the power curve is flat near its minimum: at V_mp, P = 18.2 kW), and a pilot would fly it with a margin, or with flaps.

**Check.** V_md ∝ 1/√CL* and CL* = √(πeAR CD0) grew from 0.793 to 0.929 (both CD0 and e are larger than the lesson's), so V_md falls by √(0.793/0.929) = 0.924, from 35.5 to 32.8 m/s ✓.

**Common mistakes.** Mixing the wing-alone e with the aircraft e; quoting a minimum-power speed below the stall without comment.

## 6. When is a span efficiency wrong?

(a) For a planar wing the elliptic load minimises induced drag for given lift and span, so e ≤ 1; e = 1.031 means the discrete drag sum is wrong. With control points at the y-midpoints of a cosine layout, the tip strip's control point sits very close to the tip vortex (the strip is narrow), the sampled downwash there is huge (relative error 8.2) and the sum weights it badly. The semicircle layout (control points at the φ-midpoints) has the right relation between control points and vortices.

(b) Changed a lot: the downwash and sidewash at a tail or a following aircraft (the vortices move inboard and down). Hardly changed: the induced drag, which depends on the shed circulation, not on where the vortices go (Trefftz).

(c) "The lesson's induced drag is exact only for the elliptic load on a planar wing; the trainer's rectangular wing has e ≈ 0.93–0.97 (`P10.10`, `P10.11`), and the fuselage, trim and the lift-dependent part of the profile drag lower the aircraft's Oswald factor further. Use a computed wing e times an estimated aircraft factor, or measure the polar in a glide test."

## Check-yourself answers

1. A vortex filament has constant strength along its length, and it cannot end in the fluid (it closes, goes to infinity or ends on a surface).
2. Γ/(2πh) for the infinite line; Γ/(4πh) in the plane through the start of a semi-infinite line.
3. Because Γ0 sin θ is the single term whose Glauert integral is a constant: the trailing sheet's strength varies along the span in exactly the way that makes its induced velocity uniform.
4. CDi = CL²/(πeAR). D_i = L²/(½ρV²πeb²): with L = W and V fixed, only b (and e) appear; the area cancels.
5. At the lifting line each trailing leg extends only downstream (half of an infinite line); far downstream it extends both ways.
6. Equal at best L/D; induced drag three times parasite at minimum power.
7. A planar wing cannot have e > 1: the computation is in error (layout, sampling or a sign).
8. Its two concentrated trailing vortices have finite circulation in zero area: the kinetic energy per unit length of wake, and so the Trefftz drag, diverges logarithmically at the cores.
