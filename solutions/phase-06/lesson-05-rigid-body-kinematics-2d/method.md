# P06.05  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. `grashof(*links)` for each, or by hand. 2. `position_circle_intersection(radians(10 s), L, ±1)`; `transmission_angle`. 3. `velocity_analysis(θ₂, θ₃, θ₄, 10, L)`; |v_A| = ω₂r₂, |v_B| = |ω₄|r₄; compute v_A + ω₃ × (B − A) and compare with ω₄ × (B − O₄).

## 2.
(a) Real and imaginary parts of (3); differentiate each by θ₃ and θ₄; det = r₃r₄(sin θ₃ cos θ₄ − cos θ₃ sin θ₄). (b) d/dt(r e^{iθ}) = r iω e^{iθ}; second derivative r(iα − ω²)e^{iθ}; move the known terms (θ₂ terms and the ω² terms) to the right. (c) Velocity of a point as part of a body is perpendicular to its line to that body's instant centre relative to the ground; equal velocities need the same direction. Then use I₂₄ on the ground line. (d) Show θ₃ = 0 and θ₄ = θ₂ solve (3); substitute into (5).

## 3.
(a) Grashof with crank shortest. (b) Triangle O₂O₄B with sides r₁, r₄, r₃ ± r₂; angle at O₄; θ₄ = π − that angle (B above the ground line). Check with `sweep_positions` over 0–360°. (c) Direction of O₂→B (extended) or B→O₂ (folded) gives θ₂; difference of the two crank angles; ratio of the larger to the smaller. (d) ω₂ = 2π × 5; `velocity_analysis`, `acceleration_analysis` along the sweep; `transmission_angle`.

## 4.
1. `grashof`, `position_circle_intersection(π/2, L, 1)`. 2. `sweep_positions` over π/2 ± 45°; gearing from `velocity_analysis` with ω₂ = 1; `np.polyfit` of degree 1; `transmission_angle`. 3. θ₂(t) = θ₂,₀ + ω₂(t − t₀) + ½α₂(t − t₀)²; three-point histories; `central_rates`; `observed_order`.

## 5.
(a) For each horn: travel, gearing range (max/min), minimum μ, and T = 0.3 × max gearing. (b) Free play = 1° × gearing at neutral. (c) A change-point linkage can fold into its crossed branch. (d) Combine.
