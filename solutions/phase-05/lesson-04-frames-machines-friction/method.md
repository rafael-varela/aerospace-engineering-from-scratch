# P05.04  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. μ = 0.4 + s/100. 2. All wheels braked: along the slope W sin α ≤ μ W cos α, so α_max = arctan μ. 3. Mains only: friction ≤ μ × 0.9 W cos α, so α_max = arctan(0.9 μ). 4. W cancels; the load share on the mains changes with slope and CG height.

## 2.
(a) Element dθ: normal force N = T dθ (to first order); tangential balance dT = μN; integrate dT/T = μ dθ. (b) Contact point moves until the resultant is at the friction angle to the normal; perpendicular distance from the centre r sin φ_f; sin(arctan μ) = μ/√(1 + μ²). Pulley: moments about the centre (T_t − T_s)R = r_f |T_t + T_s| (vector sum, law of cosines); divide by T_s; square; quadratic in k; root > 1; at θ = π the resultant is T_t + T_s, so k − 1 = ρ(k + 1). (c) δW = F L_g δφ + H δδ = 0 with δδ = (dδ/dφ) δφ.

## 3.
1. Centre pin: two link forces at ±θ to the line and P; perpendicular balance 2 T sin θ = P. 2. Slider: line component T cos θ = Q. 3. Divide. 4. Evaluate at three angles. 5. Virtual work P δy = Q δx with numerical derivatives. 6. Over-centre: the load drives the toggle against a stop instead of folding it.

## 4.
1. dδ/dφ ≈ −(a_s/b₁)(b₂/h). 2. `Linkage(b2=...)`; `gearing`; `solve_kinematics` at ±20°; `linkage_forces(lk, phi)["F"]`. 3. `brentq` on δ(φ = 20°; b₂) + 25°; evaluate the same quantities at the root. 4. Asymmetry: the horn and second arm swing through arcs, so equal stick angles give unequal elevator angles.

## 5.
1. Single ratios with `pulley_ratio` and `capstan_ratio`. 2. Build the run lists; `cable_run_ratio`; F_drive = F_ideal k, F_back = F_ideal/k; band = difference. 3. Replace the fourth pulley with a fairlead of 90° at μ = 0.3; recompute.

## 6.
1. Pulleys: compare the plain and ball-bearing runs of Problem 5. 2. Pin: the friction moment ∝ r, so doubling the bolt diameter at that one pivot roughly doubles its share of the band; use the lesson's 12.7 mm result as an upper bound (all pins). 3. Discuss breakout, centring, scaling with q, wear, cold and contamination. 4. Propose a stick-force-versus-position test with loads applied both ways.
