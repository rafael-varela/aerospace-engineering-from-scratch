# P02.06  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for $s = 15$; numerical results are from `solution.py` and given to the figures shown. All simulated records are synthetic.

## 1. Your own swing test ($s = 15$: m = 4.5 kg, d = 0.30 m, I_cg = 0.25 kg m²)

**Assumptions.** Rigid body, frictionless pivot, no damping, g = 9.80665 m/s².

**Derivation and result.** (a) I_pivot = 0.25 + 4.5 × 0.09 = **0.655 kg m²**; L_eff = 0.655/(4.5 × 0.30) = **0.4852 m**; ω₀ = √(9.80665/0.4852) = **4.4958 rad/s**; T₀ = 2π/ω₀ = **1.3976 s**. (b) At 15° (0.2618 rad): T_exact = **1.4036 s** (T/T₀ = 1.00430), inferred I_cg = **0.25565 kg m²**, relative error **+2.26 %**. (c) Lever ratio 0.655/0.25 = **2.62**; estimate 2.62 × 0.2618²/8 = 2.24 %, within 1 % of the exact 2.26 % (the θ₀⁴ term). (d) `max_amplitude_for_inertia_error(0.01, case)` gives 0.1745 rad = **10.0°**.

**Check.** Dimensions: m g d T² has units kg m s⁻² m s² = kg m² ✓. Limiting case: I_cg → 0 makes L_eff → d, the point-mass pendulum. The lever ratio 2.62 is smaller than the lesson's 3.45 because this airframe has more inertia relative to m d², so it tolerates a larger amplitude (10.0° against 8.7°).

**Common mistakes.** Using the amplitude in degrees in θ₀²/8; forgetting to subtract m d² when inferring I_cg (then the error is only 0.86 %, the I_pivot error, and the lever amplification is missed).

## 2. The inverted equilibrium and the period series

**Derivation.** (a) With **y** = (θ, θ̇): **f** = (θ̇, −2ζω₀θ̇ − ω₀² sin θ). Equilibria: θ̇ = 0 and sin θ = 0, i.e. θ = kπ. J = [[∂f₁/∂θ, ∂f₁/∂θ̇], [∂f₂/∂θ, ∂f₂/∂θ̇]] = [[0, 1], [−ω₀² cos θ, −2ζω₀]].

(b) At θ = π, J = [[0, 1], [ω₀², −2ζω₀]]. det(λI − J) = λ(λ + 2ζω₀) − ω₀² = λ² + 2ζω₀λ − ω₀² = 0, so λ = −ζω₀ ± √(ζ²ω₀² + ω₀²) = ω₀(−ζ ± √(1 + ζ²)). Since √(1 + ζ²) > ζ for every ζ ≥ 0, λ₊ > 0: the product of the roots is −ω₀² < 0, which is the saddle condition. A tilt u grows as e^{λ₊t}, doubling in t₂ = ln 2/λ₊. UAV (ω₀ = 4.4607 s⁻¹): ζ = 0 gives λ₊ = 4.4607, **t₂ = 0.1554 s**; ζ = 0.01 gives λ₊ = 4.4163, **t₂ = 0.1570 s**; ζ = 0.2 gives λ₊ = 4.4607 × (−0.2 + 1.0198) = 3.6569, **t₂ = 0.1895 s**. Damping slows the fall but can never stop it: no viscous damping makes a saddle stable.

(c) Energy: ½θ̇² + ω₀²(1 − cos θ) = ω₀²(1 − cos θ₀), so θ̇ = ω₀√(2(cos θ − cos θ₀)). Quarter period: T/4 = ∫₀^θ₀ dθ/θ̇. With cos θ = 1 − 2 sin²(θ/2): cos θ − cos θ₀ = 2[sin²(θ₀/2) − sin²(θ/2)] = 2[k² − sin²(θ/2)], k = sin(θ₀/2). Substitute sin(θ/2) = k sin φ (φ runs 0 → π/2 as θ runs 0 → θ₀): ½ cos(θ/2) dθ = k cos φ dφ, so dθ = 2k cos φ dφ/√(1 − k² sin² φ), and √(2(cos θ − cos θ₀)) = 2k cos φ. The integrand becomes dφ/√(1 − k² sin² φ), hence T = (4/ω₀)∫₀^{π/2} dφ/√(1 − k² sin² φ) = (4/ω₀)K(k).

(d) K = (π/2)(1 + k²/4 + …) and k² = sin²(θ₀/2) = θ₀²/4 + O(θ₀⁴). So T = (2π/ω₀)(1 + θ₀²/16 + O(θ₀⁴)) = T₀(1 + θ₀²/16 + …). The 1/16 comes from ¼ × ¼: a quarter from K's series and a quarter from k² ≈ θ₀²/4.

**Check.** At θ₀ = 90°, k² = ½ and the series gives 1 + 1/8 + 9/256 = 1.160 against the exact 1.18034: the series converges slowly at large k, as expected for a series in k². Dimensions: ω₀ in s⁻¹, so 4K/ω₀ is a time ✓.

**Common mistakes.** Writing the Jacobian at π with −ω₀² (forgetting cos π = −1); dropping the 2 in dθ = 2k cos φ dφ/cos(θ/2) and getting T = 2K/ω₀.

## 3. Gravity-gradient libration ($s = 15$: k = 0.65, n = 1.1300 × 10⁻³ rad/s)

**Assumptions.** Circular orbit, rigid satellite, pitch motion decoupled from roll and yaw, no aerodynamic or magnetic torques.

**Derivation.** (a) sin θ cos θ = ½ sin 2θ, so the equilibria in [0, π] are where sin 2θ = 0: θ = 0, π/2, π. With **y** = (θ, θ̇), J = [[0, 1], [−3n²k cos 2θ, 0]], tr J = 0, det J = 3n²k cos 2θ. At θ = 0 and π: det = +3n²k = 2.49 × 10⁻⁶ s⁻², **center** (eigenvalues ±1.578 × 10⁻³ i). At θ = π/2: det = −3n²k, **saddle** (eigenvalues ±1.578 × 10⁻³ s⁻¹). θ = 0 (and the physically identical θ = π) is the orientation with body z along the local vertical; since k > 0 means I_x > I_z, z is the axis of smaller inertia, so the stable attitude has the *long* axis pointing at the Earth. The undamped equation has an energy integral ½θ̇² + (3n²k/2) sin² θ, so the linear centers are true centers.

(b) Linearised at 0: θ̈ + 3n²k θ = 0, ω_lib = n√(3k) = 1.1300 × 10⁻³ × √1.95 = **1.578 × 10⁻³ rad/s**, T_lib = 2π/ω_lib = 3982 s = **66.4 min**; ratio T_lib/T_orbit = 1/√(3k) = **0.716**. Only the inertia ratio k enters because the gravity-gradient torque and the inertia both scale with the mass distribution in the same way; a 1 m and a 10 m satellite of the same shape librate at the same rate.

(c) With φ = 2θ: φ̈ = 2θ̈ = −3n²k sin 2θ = −3n²k sin φ: a pendulum with ω₀² = 3n²k. A 30° libration is a 60° pendulum swing: `period_exact(np.radians(60), omega_lib)` over the small-amplitude period gives **1.0732**, so the period at 30° is 7.3 % longer, 71.2 min.

(d) k < 0 means I_z > I_x: the axis along the vertical at θ = 0 is now the larger-inertia one. det J changes sign at every equilibrium: 0 and π become saddles and π/2 becomes the center. The satellite swings round until its smaller-inertia axis points at the Earth, which is the general rule.

**Check.** Dimensions: 3n²k has units s⁻² ✓. Limiting case k = 1 (a thin rod, I_z = 0): ω_lib = √3 n, the textbook value for a dumbbell satellite. The saddle's eigenvalues are ±ω_lib, just as the pendulum's saddle has ±ω₀.

**Common mistakes.** Linearising sin θ cos θ as θ cos θ and missing the cos 2θ in the Jacobian; forgetting that φ = 2θ doubles the amplitude as well as the angle when using `period_exact`.

## 4. Identifying damping from a simulated swing ($s = 15$: ζ = 0.015)

**Assumptions.** Synthetic record from `integrate_rk4` at h = 1 ms, 12 T₀ long, released from rest at 20°; linear viscous damping.

**Result.** (a) Peaks (the release counts as the first): 20.00°, 18.18°, 16.53°, 15.03°, 13.67°, 12.43° at t = 0, 1.419, 2.835, 4.251, 5.665, 7.078 s. Taking the six *detected* peaks after the release (18.18° … 11.31°, the seventh maximum being 11.31° at t = 8.491 s): δ = (1/5) ln(18.18/11.31) = **0.0950**; ζ = δ/√(4π² + δ²) = **0.0151**, against the 0.015 set (the 0.8 % excess is the amplitude dependence of the period, which the linear decrement formula does not know about).

(b) `period_from_trajectory` gives **1.41264 s**. T₀ = 1.40858 s, the exact undamped period at 20° is 1.41938 s, and the damped linear period T₀/√(1 − ζ²) = 1.40874 s. The measured value lies between T₀ and T(20°) because the amplitude decays from 20° to 7° during the record, so the mean period corresponds to a mean amplitude of about 12.5°; the damping correction to the period is only ζ²/2 = 1.1 × 10⁻⁴, a hundred times smaller than the amplitude correction (θ₀²/16 = 7.6 × 10⁻³ at 20°).

(c) `inertia_from_period(1.41264, 4.2, 0.35)` = 0.21419 kg m², error **+2.00 %**, almost all of it amplitude: with the 10° release the error is 0.56 %.

(d) Two-release procedure: T₁ = T₀(1 + θ₁²/16), T₂ = T₀(1 + θ₂²/16) with θ₁ = 20°, θ₂ = 10° gives T₀ = (T₁θ₂² − T₂θ₁²)/(θ₂² − θ₁²). With T₁ = 1.41264 s and T₂ = 1.40971 s: T₀ = 1.40874 s, and the inferred I_cg is 0.21017 kg m², error **+0.08 %**. The residual is exactly the damping correction, ζ²/2 = 1.1 × 10⁻⁴ in T times the lever factor 6.9 = 0.078 %. (Using the record's rms amplitude instead of the nominal release gives 1.40873 s, the same answer: the correction is insensitive to how the amplitude is defined because only the *difference* of the two corrections matters.)

**Check.** The decrement predicts the sixth peak from the first: 18.18° × e^{−5 × 0.0950} = 11.31° ✓. Energy at the last peak is (7.05/20)² = 12 % of the release energy, consistent with e^{−2ζω₀ t} = e^{−2 × 0.015 × 4.46 × 15.5} = 0.125 at t = 11 T₀.

**Common mistakes.** Taking the first detected sample peak as peak 1 when the release itself is the true first maximum; attributing the period shift to damping (it is O(ζ²), negligible) instead of to amplitude.

## 5. What limits a real swing test

**Derivation.** (a) I_cg = m g d T²/(4π²) − m d². ∂I_cg/∂T = 2 m g d T/(4π²) = 2 I_pivot/T = 2 × 0.7245/1.4086 = **1.0287 kg m²/s**; ±2 ms gives ±0.00206 kg m², **±0.98 %** of I_cg. ∂I_cg/∂d = m g T²/(4π²) − 2 m d = I_pivot/d − 2 m d = 2.070 − 2.940 = **−0.870 kg m²/m**; ±2 mm gives ∓0.00174 kg m², **±0.83 %**.

(b) Amplitude bias: 5°: **+0.33 %**; 10°: **+1.32 %**; 15°: **+2.97 %**.

(c) Ranked at the nominal 10°: amplitude (+1.3 %, and up to +3.0 % at 15°), timing (±1.0 %), tape (±0.8 %). The amplitude error is the largest *and* it is a bias, not a random error: averaging repeated swings does not reduce it. It is also the cheapest to remove: read the release angle on a protractor and apply the correction T₀ = T/(1 + θ₀²/16), or release at 5°. Money, if spent, goes on a photogate or an IMU to time the period (10⁻⁴ s rather than 2 × 10⁻³) because after the amplitude correction the timing error is the largest remaining one. The tape error is harder to reduce because d is the distance to a centre of mass that itself is located by a separate balance test.

(d) Pivot compliance: a knife edge resting on a flexible beam adds a second period and shifts the first. The revealing measurement is the period at two suspension distances d₁, d₂: the rigid model predicts T₁²d₁ − T₂²d₂ = 4π²(d₁² − d₂²)/g exactly (I_cg cancels), and a compliant pivot breaks the identity. Quadratic damping is revealed by a log decrement that falls with amplitude.

**Check.** The three relative errors add in quadrature to about 1.8 % at 10° release, which is why swing tests quote I to about 2 %. Dimensions of ∂I/∂d: kg m²/m = kg m ✓ (and I_pivot/d has the same).

**Common mistakes.** Treating the amplitude bias as random and dividing it by √N; forgetting the −2 m d term in ∂I/∂d (its sign makes the two contributions to the sensitivity partly cancel).

## 6. The Van der Pol limit cycle

**Derivation.** (a) With **y** = (x, ẋ): **f** = (ẋ, μ(1 − x²)ẋ − x), J(0) = [[0, 1], [−1, μ]], tr J = μ, det J = 1, discriminant μ² − 4. μ = 0.1: **unstable spiral** (λ = 0.05 ± 0.9987i); μ = 1: **unstable spiral** (0.5 ± 0.866i); μ = 3: **unstable node** (λ = 0.382, 2.618). In every case the origin repels.

(b) `limit_cycle_amplitude`: μ = 0.1: **2.0001** from (0.5, 0) and **2.0001** from (3, 0); μ = 1: **2.0086** and **2.0086**; μ = 3: **2.0233** and **2.0233**. The same orbit is reached from inside and outside, and its amplitude moves away from 2 only at O(μ²).

(c) With x = A cos t, ẋ = −A sin t: W = ∮ μ(1 − x²)ẋ² dt = μA² ∫₀^{2π} (1 − A² cos² t) sin² t dt = μA²[π − A² (π/4)] = **μπA²(1 − A²/4)**. W = 0 only at A = 2 (and the trivial A = 0); W > 0 for A < 2 (the "negative damping" near the origin pumps energy in and the amplitude grows), W < 0 for A > 2 (the x² term removes energy and the amplitude shrinks). Both sides move toward A = 2, so the cycle attracts from inside and outside: it is a stable limit cycle.

(d) A linearisation keeps only the terms of first order in the displacement from one point and is therefore valid in a small neighbourhood of that point; the amplitude of the cycle is fixed by an energy balance along a whole closed orbit on which the cubic term x²ẋ is as large as the linear ones. Nothing evaluated at a single point can contain that information.

**Check.** The integral ∫₀^{2π} sin² t cos² t dt = ∫ ¼ sin² 2t dt = π/4 ✓. For μ = 0.1 the measured amplitude 2.0001 confirms the leading-order 2 to 5 × 10⁻⁵, consistent with a correction of order μ² = 0.01 times a small coefficient.

**Common mistakes.** Calling the origin "stable" because the cycle is stable; integrating the work with x = A sin t and ẋ = A cos t and dropping the sign (the answer is the same, but only if the signs are kept consistently).

## Check-yourself answers

1. An equilibrium is a state **y*** with **f**(**y***) = **0**, where the system can remain forever; the Jacobian is the matrix of first partial derivatives J_ij = ∂f_i/∂y_j at **y***, so that **u̇** ≈ J**u** for a small displacement **u**. For the pendulum J(θ) = [[0, 1], [−ω₀² cos θ, −2ζω₀]].
2. [[0, 1], [−9, −0.6]]: tr = −0.6, det = 9, disc = 0.36 − 36 < 0: stable spiral. [[0, 1], [9, −0.6]]: det = −9 < 0: saddle.
3. A center has pure imaginary eigenvalues, so the neglected quadratic terms can turn the closed linear orbits into an inward or outward spiral; the linearisation is silent (Hartman–Grobman needs no eigenvalue on the imaginary axis). For the undamped pendulum the conserved energy E = ½θ̇² + ω₀²(1 − cos θ) forces the orbits to be closed level curves, so the center is real.
4. E_sep = 2ω₀² (the energy of the inverted pendulum at rest). A state with E > 2ω₀² goes over the top and rotates in one direction forever (undamped).
5. 1 + (π/2)²/16 = 1.154 against the exact 1.180: the series in θ₀² is poor at 90° because k² = ½ is not small and K(k) has a logarithmic singularity at k = 1.
6. θ₀ = 0.2618 rad: error ≈ 3.45 × 0.2618²/8 = 2.96 % (exact 2.97 %).
7. Fourth order: 10⁻⁶/16 = 6 × 10⁻⁸ at h/2 and 1.6 × 10⁻⁵ at 2h.
8. Wing flutter settling at a fixed amplitude (or a Van der Pol circuit, or a hunting autopilot): a limit cycle. The linearisation at the equilibrium shows an unstable spiral (or node), i.e. that the equilibrium repels, but nothing about where the trajectory ends up.
