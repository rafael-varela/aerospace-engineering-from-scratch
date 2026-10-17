# P04.04  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for **s = 15**; numerical results are from `solution.py` (and the probe scripts that call it) and are given to the figures shown. All "simulated" values are outputs of the model, not measurements.

## 1. A reaction wheel (s = 15: m = 3.5 kg)

**Assumptions.** Uniform solid disk; rigid; balanced; the spacecraft and wheel start at rest.

(a) I₃ = ½ × 3.5 × 0.1² = **0.0175 kg·m²**; I₁ = 3.5(3 × 0.01 + 0.0004)/12 = **0.008867 kg·m²**.

(b) ω = 6000 × 2π/60 = 628.3 rad/s; h = I₃ω = **10.996 N·m·s**; T = ½I₃ω² = **3454 J**. The motor's electrical energy becomes the wheel's kinetic energy; when the wheel is slowed, the motor either returns it to the bus (regenerative braking) or dissipates it as heat.

(c) `parallel_axis(I, 3.5, [0.5, 0, 0])` = **diag(0.008867, 0.8839, 0.8925) kg·m²**: m|d|² = 0.875 kg·m² is added about axes 2 and 3 and nothing about axis 1, which passes through the wheel's centre. There are **no products of inertia**, because the offset lies along a principal axis of the wheel, so ddᵀ is diagonal.

(d) t = h/τ = 10.996/0.1 = **110.0 s**. The total angular momentum about axis 3 stays zero, so the body turns the other way at 10.996/135.4 = **0.0812 rad/s = 4.65 °/s**.

**Check.** Perpendicular-axis theorem for a thin disk: I₁ + I₂ = I₃ when ℓ → 0; here 2 × 0.008867 = 0.01773 is slightly more than 0.0175 because of the 2 cm thickness (+ mℓ²/6).

**Common mistakes.** Using rpm directly in ½Iω²; forgetting that the spacecraft counter-rotates.

## 2. Euler's equations and their consequences

(a) In the inertial frame dL/dt = τ. By (4), dL/dt|i = dL/dt|b + ω × L, and in body axes I is constant, so dL/dt|b = Iω̇. Hence **Iω̇ + ω × (Iω) = τ**. In principal axes: **I₁ω̇₁ − (I₂ − I₃)ω₂ω₃ = τ₁, I₂ω̇₂ − (I₃ − I₁)ω₃ω₁ = τ₂, I₃ω̇₃ − (I₁ − I₂)ω₁ω₂ = τ₃.**

(b) With τ = 0, dot with ω: ω · Iω̇ = −ω · (ω × Iω) = 0, and since I is symmetric ω · Iω̇ = d(½ω·Iω)/dt, so **T is constant**. Dot with Iω: (Iω) · Iω̇ = −(Iω) · (ω × Iω) = 0, i.e. d(½|Iω|²)/dt = 0: **|L| is constant**. In the inertial frame, L_i = RL and dL_i/dt = R(dL/dt|b + ω × L) = Rτ = 0: **L_i is constant**, the stronger statement (direction too).

(c) Let ω = (ε₁, ω₂ + ε₂, ε₃) with ε small. Keeping first-order terms: I₁ε̇₁ = (I₂ − I₃)ω₂ε₃, I₂ε̇₂ = 0, I₃ε̇₃ = (I₁ − I₂)ω₂ε₁. Differentiating the first and substituting the third: **ε̈₁ = (I₂ − I₃)(I₁ − I₂)ω₂²ε₁/(I₁I₃)**. If I₂ lies between I₁ and I₃, both brackets have the same sign, the coefficient is positive and ε grows as e^{σt} with **σ = |ω₂|√((I₁ − I₂)(I₂ − I₃)/(I₁I₃))**. The same calculation about axis 1 gives ε̈₂ = −(I₁ − I₃)(I₁ − I₂)ω₁²ε₂/(I₂I₃); if I₁ is the largest or the smallest, the product is positive and the motion is an oscillation at **√((I₁ − I₂)(I₁ − I₃)/(I₂I₃))|ω₁|**. For the box at 0.5 rad/s: 0.2402 rad/s about axis 1, 0.2049 rad/s about axis 3, and σ = 0.1641 s⁻¹ about axis 2.

(d) Body axes with axis 1 horizontal: the vertical is ẑ = (0, sin θ, cos θ), ω = (0, Ω sin θ, ω₃), L = (0, I₁Ω sin θ, I₃ω₃). Steady precession at Ω about the vertical means dL/dt = Ω ẑ × L, whose only non-zero component is the first: Ω[sin θ · I₃ω₃ − cos θ · I₁Ω sin θ]. The gravity torque (l e₃) × (m g_b), g_b = −g ẑ, is (mgl sin θ, 0, 0). Equating: Ω sin θ(I₃ω₃ − I₁Ω cos θ) = mgl sin θ, i.e. **I₁ cos θ Ω² − I₃ω₃Ω + mgl = 0**. Real roots need **(I₃ω₃)² ≥ 4I₁mgl cos θ**. The slow root written as Ω = 2mgl/(I₃ω₃ + √((I₃ω₃)² − 4I₁mgl cos θ)) tends to 2mgl/(2I₃ω₃) = **mgl/(I₃ω₃)** as ω₃ → ∞, which is (12).

(e) With Euler angles (φ precession, θ tilt), for a symmetric top E = ½I₁(θ̇² + φ̇² sin²θ) + ½I₃ω₃² + mgl cos θ and L_z = I₁φ̇ sin²θ + L₃ cos θ, with ω₃ and L₃ = I₃ω₃ constant. At release θ̇ = φ̇ = 0: L_z = L₃ cos θ₀, E = ½I₃ω₃² + mgl cos θ₀. At another turning point θ̇ = 0: φ̇ = L₃(cos θ₀ − cos θ)/(I₁ sin²θ), and energy gives ½I₁φ̇² sin²θ = mgl(cos θ₀ − cos θ). Substituting: L₃²(cos θ₀ − cos θ)²/(2I₁ sin²θ) = mgl(cos θ₀ − cos θ); dividing by the non-zero factor (cos θ₀ − cos θ), **cos θ₀ − cos θ = 2I₁mgl sin²θ/L₃²**. For small θ − θ₀ the left side is ≈ sin θ₀ (θ − θ₀), so θ − θ₀ ≈ 2I₁mgl sin θ₀/L₃².

**Checks.** (c): a sphere (all moments equal) gives σ = 0 and no wobble; (d): at θ = 90° the quadratic is linear and the root is exactly (12); (e): L₃ → ∞ gives no nod.

**Common mistakes.** Using I_cm instead of the moment about the pivot in the top; dropping the I₁Ω cos θ term and then being puzzled by the 3 % difference; taking the larger root of (13).

## 3. A pencil-shaped satellite (s = 15: m = 13 kg)

**Assumptions.** Uniform cylinder; no external torque; the antennas dissipate energy but are light enough not to change the inertia appreciably.

(a) I₃ = ½ × 13 × 0.08² = **0.0416 kg·m²**; I₁ = I₂ = 13(3 × 0.0064 + 4)/12 = **4.354 kg·m²**. ω₀ = 750 rpm = 78.54 rad/s; L = I₃ω₀ = **3.267 N·m·s**; T₀ = ½I₃ω₀² = **128.3 J**.

(b) **Yes**, in the rigid model: the long axis has the smallest moment, an extreme axis, and small disturbances oscillate.

(c) At fixed |L| = L, the energy of spin about a principal axis i is L²/(2Iᵢ); for a general state T ≥ L²/(2I_max) (write T = Σ Lᵢ²/(2Iᵢ) with ΣLᵢ² = L²). The minimum is spin about the axis of largest moment, a transverse axis: ω_f = L/I₁ = **0.7504 rad/s = 7.17 rpm**, T_f = **1.226 J**. The antennas dissipate **127.1 J**, 99 % of the initial energy. The satellite ends up tumbling end over end at a slow rate, with the same angular momentum.

(d) **A passively spin-stabilised body must spin about its axis of maximum moment of inertia** (or carry an active damper that keeps it from dissipating its way off a minor-axis spin).

**Check.** T_f/T₀ = I₃/I₁ = 0.00955, as L²/(2I) predicts.

**Common mistakes.** Conserving energy instead of angular momentum; concluding from (b) that the satellite is safe.

## 4. Slewing a spacecraft (s = 15: 35° in 60 s)

**Assumptions.** Rigid body; ideal wheels; no external torque; the axis-3 wheel's inertia is lumped with the body.

(a) `bang_bang_slew(135.42, 35°, 60)`: τ = **91.91 mN·m**, ω_max = **1.167 °/s**, h_max = **2.757 N·m·s**.

(b) Final angle **35.0000°** (2 atan2(q_z, q_w)), final body rate **9 × 10⁻¹² rad/s**; Iω + h is zero at every sample to 10⁻¹² N·m·s.

(c) With h₁ = 2 N·m·s: max |ω₁| = **1.31 × 10⁻³ rad/s**, max |ω₂| = **5.17 × 10⁻³ rad/s**, final tilt of axis 3 **9.55°**, final body rate **5.8 × 10⁻³ rad/s**. With h₁ = 10 N·m·s: **2.77 × 10⁻³** and **1.28 × 10⁻² rad/s**, final tilt **21.1°** (largest 21.6°), final body rate **2.7 × 10⁻² rad/s**. The total inertial angular momentum stays at its initial (h₁, 0, 0) to 10⁻¹⁴. Explanation: in (10) the term ω × h with ω along axis 3 and h along axis 1 is a gyroscopic torque along axis 2, proportional to h₁ω₃; it tilts axis 3 and leaves the body tumbling after the wheel torque stops. Physically, the bias momentum is fixed in space, so the body cannot simply rotate about axis 3 and carry the bias with it. A momentum-biased spacecraft must slew with all three axes controlled, or about the bias axis.

(d) Momentum limit: h_max = 0.02 × 628.3 = 12.57 N·m·s; T ≥ 2IΘ/h_max = 2 × 135.4 × (π/2)/12.57 = **33.9 s**. Torque limit: T ≥ √(4IΘ/τ_max) = √(4 × 135.4 × 1.571/0.1) = **92.2 s**. The **torque binds**: at 92.2 s the peak wheel momentum is only 4.61 N·m·s, well under 12.57.

**Check.** (a): τ(T/2)²/I = 0.0919 × 900/135.4 = 0.611 rad = 35°.

**Common mistakes.** Applying the wheel torque to the body with the same sign; letting the integrator step over the torque reversal (hence `max_step`); assuming the slew in (c) stays about axis 3.

## 5. Your own top (s = 15: 55 rev/s)

**Assumptions.** Frictionless pivot; rigid disk; I₁ = 1.8263 × 10⁻⁴ kg·m² about the pivot, I₃ = 4.5 × 10⁻⁵ kg·m², mgl = 0.03923 N·m.

(a) ω₃ = 345.58 rad/s. Fast-top **2.5225 rad/s**; steady root **2.5907 rad/s** (period **2.43 s**); simulated **2.5907 rad/s** (equal to the root to 10⁻⁹).

(b) | tilt | 10° | 30° | 60° | 80° | 89° |
|---|---|---|---|---|---|
| minimum spin (rev/s) | 18.8 | 17.6 | 13.4 | 7.9 | 2.5 |

At 90° the term I₁Ω² cos θ vanishes, (13) becomes linear and has the root mgl/(I₃ω₃) for any spin.

(c) Tilt swing **1.8386°**; turning point from (15) **1.8386°**; fast-top estimate **1.697°** (8 % low); nutation frequency from successive tilt maxima **80.7 rad/s** against L₃/I₁ = **85.1 rad/s**; mean precession **2.593 rad/s**, close to the steady root.

(d) **17.6 rev/s** (b). Below it no cone of steady precession exists at 30°: the nods grow, the axis drops, and the top falls over or its rim touches the ground.

**Check.** The released top's mean precession and the steady root agree to 0.1 %; the swing equals (15) to the printed digits.

**Common mistakes.** Using I₁ about the centre of mass; reading the nutation frequency from a coarse FFT (a 1 s record has bins 6.3 rad/s wide); starting the "steady" case with Ω = 0.

## 6. Judgment

(a) The torque adds angular momentum at 10⁻⁵ N·m·s per second, and the wheels must absorb it to keep the telescope still; a 50 N·m·s wheel is full after 5 × 10⁶ s ≈ **58 days**. Then an external torque (magnetic torquers or thrusters) must "dump" the momentum.

(b) Angular momentum about the vertical is conserved (no external torque about it, friction aside). Kinetic energy is not: T = L²/(2I) rises as I falls. The extra energy is the work her arms do pulling inward against the centrifugal tendency.

(c) Its rotors spin in opposite pairs. Speeding up one pair and slowing the other changes the net rotor angular momentum, and the body takes the reaction torque about the vertical, while total thrust is held constant.

(d) Any of: the torque expressed in the wrong frame; a sign error in ω × Iω; R computed from an unnormalised quaternion; a loose solver tolerance or a fixed-step integrator with too large a step.

(e) Gravity's torque is what makes the axis precess; without gravity the axis would stay fixed. When the top is released without precession, the centre of mass drops slightly (the nod) and gravity's work goes into the kinetic energy of the precession and the nodding.

(f) Energy dissipation at fixed L drives any body towards its major axis; a minor-axis spinner (long axis) therefore needs an active damper that removes nutation without letting the spin migrate.
