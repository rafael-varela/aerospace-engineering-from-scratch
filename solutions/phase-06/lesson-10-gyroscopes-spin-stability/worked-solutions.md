# P06.10  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for **s = 15**; numerical results are from `solution.py` (and the probe scripts that call it) and are given to the figures shown. All "simulated" values are outputs of the model, not measurements.

## 1. Nutation of an oblate satellite (s = 15: ω₃ = 1.25 rad/s)

**Assumptions.** Rigid, axisymmetric, torque-free.

(a) H = √((1200 × 1.25)² + (800 × 0.02)²) = √(1500² + 16²) = **1500.085 N·m·s**; λ = (1200 − 800) × 1.25/800 = **0.625 rad/s**; Ω_n = H/I_t = **1.875107 rad/s**; θ = atan(16/1500) = **0.6111°**.

(b) Body-frame period 2π/λ = **10.05 s**; inertial coning period 2π/Ω_n = **3.351 s**. The ground sees the axis cone about the fixed H at the inertial rate: a 3.35 s wobble. The 10.05 s motion is what an observer riding on the satellite sees the rate vector do.

(c) Simulated coning rate over 40 s: **1.875107 rad/s**, equal to H/I_t to 10⁻¹² relative.

**Check.** Ω_n − λ = 1.875 − 0.625 = 1.250 ≈ ω₃ (for small θ the inertial coning rate is the spin plus the body nutation rate: the body turns at ω₃ while the transverse rate turns back at λ). 

**Common mistakes.** Quoting λ as the wobble seen from the ground; using degrees in tan θ.

## 2. Energy, dissipation and the damper

(a) The component of H along the symmetry axis is I₃ω₃ = H cos θ and across it I_tω_t = H sin θ. So T = ½I₃ω₃² + ½I_tω_t² = **(H²/2)(cos²θ/I₃ + sin²θ/I_t)**. For any rigid body, in principal axes T = Σ Hₖ²/(2Iₖ) ≥ Σ Hₖ²/(2I_max) = H²/(2I_max), with equality only when every Hₖ with Iₖ < I_max is zero: spin about the axis of largest moment.

(b) dT/dt = (H²/2)(2 sin θ cos θ)(1/I_t − 1/I₃) θ̇ = −P, so **θ̇ = −P/(H² sin θ cos θ (1/I_t − 1/I₃))**. Prolate (I₃ < I_t): 1/I_t − 1/I₃ < 0, so θ̇ > 0 for P > 0 (nutation grows). Oblate: θ̇ < 0 (it decays).

(c) Adding the two equations of (8): d(Iω)/dt + d(Jω_d)/dt + ω × (Iω + Jω_d) = 0, i.e. dH/dt|body + ω × H = 0, which by the transport theorem is dH/dt|inertial = 0. Dot the first with ω: d(½ω·Iω)/dt = c ω·(ω_d − ω); the second with ω_d: d(½J|ω_d|²)/dt = −c ω_d·(ω_d − ω) (the term ω_d·(ω × ω_d) is zero). Adding: **−c|ω_d − ω|²**.

(d) With ω = (z, Ω), ω_d = (z_d, Ω) to first order (ω₃ and ω_d3 are constant to first order), (ω × Iω) in complex form is −i(I₃ − I_t)Ωz and (ω × ω_d) is iΩ(z_d − z) (lesson 2.5); substituting gives (9). For c = 0: **i(I₃ − I_t)Ω/I_t**, the body's free nutation (1); and **−iΩ**, the sphere at rest in space seen from the body turning at Ω.

(e) For large c the sphere is dragged along: z_d → z. Add I_t × (first row) and J × (second): (I_t + J)ż = i(I₃ − I_t)Ωz − iJΩ(z_d − z) + … → i(I₃ − I_t)Ωz. So the eigenvalue tends to **i(I₃ − I_t)Ω/(I_t + J)**, purely imaginary: body and sphere move as one rigid body of transverse inertia I_t + J (and axial I₃ + J), z_d − z → 0, and the dissipation c|z_d − z|² falls to zero. A probe gives −1.73077i s⁻¹ (real part 5 × 10⁻⁷ s⁻¹) at c = 10⁷ for the lesson's prolate spinner, against −1.73077i from the formula.

**Common mistakes.** Writing the sphere's angular momentum as Jω_d only after forgetting that it, too, needs the ω × term in body axes.

## 3. A spinning upper stage (s = 15: I₃ = 175, I_t = 525 kg·m²)

**Assumptions.** Rigid; constant thrust and offset; no propellant motion in (a)–(c).

(a) τ = 60 000 × 0.001 = **60 N·m**; ω₃ = 2π rad/s, H = 175 × 2π = **1099.6 N·m·s**; (5): 2τ/(I₃ω₃²) = 2 × 60/(175 × 39.48) = **0.01737 rad = 0.9952°**.

(b) Simulated over 80 s: largest wander of H **0.9950°**; largest tilt of the spin axis from its start **1.297°**. The axis also nutates around the wandering H: the body-fixed torque excites a transverse rate, so the axis cones about H with a small half-angle of its own, which adds to the wander.

(c) From rest about a transverse axis, ½(τ/I_t)t² = ½ × (60/525) × 80² = 365.7 rad ≈ **20 954°**: about 58 full turns. Spin turns an unacceptable drift into a 1° wobble.

(d) e^{80/300} = **1.31**: during the burn a small nutation grows by about a third, which is acceptable if the initial nutation is small. After the burn, if the stage stays attached to the payload and spinning, the growth continues without limit: the payload must separate (or the spin be removed) within a few e-folding times, or the stack will head for a flat spin.

**Check.** The simulated H wander matches (5) to 0.02 %, so the averaging formula is reliable for this body.

**Common mistakes.** Using I_t in (5) (it is the axial H that resists); forgetting that the axis error is H's wander plus nutation.

## 4. A gyroscope (s = 15: τ = 2.5 N·m, θ = 35°)

(a) Ω = τ/(H sin θ) = 2.5/(1200 × sin 35°) = **3.63218 × 10⁻³ rad/s**; period **1729.9 s** (28.8 min).

(b) Fitted: **3.63218 × 10⁻³ rad/s** (equal to (4) to 5 × 10⁻¹¹ relative); as in the lesson, |H| and the 35° tilt stay constant to the solver tolerance because the torque is perpendicular to H and to z.

(c) τ/H = **2.0833 × 10⁻³ rad/s**, 43 % low. H's tip moves at τ along a horizontal circle of radius H sin θ, not H, so the angular rate is τ/(H sin θ); τ/H is right only at θ = 90°.

**Common mistakes.** Using the fast-top formula τ/H at a tilt; thinking the torque tilts H downwards (it moves it sideways).

## 5. Tuning a damper (s = 15: c = 60 N·m·s)

(a) Growth rates (s⁻¹): **0.00976 (c = 10), 0.02472 (30), 0.03324 (60), 0.03202 (100), 0.02199 (200), 0.01581 (300)**. The maximum is **0.03368 s⁻¹ at c ≈ 71 N·m·s**: weak damping dissipates little, strong damping locks the sphere to the body.

(b) Fitted early growth rate **0.03315 s⁻¹** (against 0.03324 s⁻¹ from (9), 0.3 % low); the axis passes **45° at 214.6 s** and **89° at 327.6 s**; final energy **1256.5392 J** against H²/(2(I_t + J)) = **1256.5392 J** (from 2970.0 J); H conserved to 8.5 × 10⁻¹⁰.

(c) ln(0.2/10⁻³)/0.03324 = **159.4 s**. The simulation reaches 45° at 214.6 s, later, because 0.2 rad/s is still a small nutation (atan(1500 × 0.2/(600 × 3)) = 9.5°): reaching 45° needs a few more e-folds, and the nonlinear growth slows as the motion leaves the linear range.

(d) For an oblate spinner the damper should maximise the *decay* rate: for I_t = 800, I₃ = 1200 kg·m² at 1 rad/s with J = 60, the optimum is **c ≈ 83 N·m·s**, giving −0.01765 s⁻¹ (a 57 s time constant). Too little or too much viscosity both damp more slowly, for the same reason as in (a).

**Check.** The final rate equals H/(I_t + J) and the axis is at 90° from H: the flat spin of 2.4.

**Common mistakes.** Choosing the largest c "for the most damping"; fitting the growth rate over the nonlinear part of the run.

## 6. Judgment

(a) The rigid-body equations are right for a rigid body, and they correctly say long-axis spin of a rigid cylinder is stable. Explorer 1 was not rigid: its flexible antennas dissipated energy, and with dissipation the major-axis rule applies (16.07 L27).

(b) An oblate body spins about its axis of largest moment, the one stable spin with dissipation, so any energy loss (fuel slosh, flexing) damps its nutation, and a nutation damper helps.

(c) When the spin lasts only a small fraction of the nutation growth time (or an active controller holds it): Problem 3 shows a growth by e^{80/300} = 1.31 during an 80 s burn, acceptable; the stack must then separate or despin.

(d) By (7) a damper on a prolate spinner adds dissipation, which makes θ grow: it destabilises it faster (Experiment 5 of the lesson shows the oblate case decays with the same damper that drives the prolate one to a flat spin).

(e) The direction of the drift (7), the major-axis rule and the end state H²/(2I_max) depend only on energy and angular momentum, not on the mechanism. The rate (the growth rate, the time to a flat spin) depends on the sink's power, so on the mechanism.

(f) The cylinder is uniform, which Explorer 1 was not, and the damper (J = 0.02 kg·m², c = 0.01 N·m·s) is an illustrative stand-in for the antennas, chosen to make the simulation short; the real dissipation rate was never modelled here, so 42 s says nothing about the real satellite's timescale.

## Check-yourself answers

1. Body-frame λ = (I₃ − I_t)ω₃/I_t; inertial H/I_t. For an oblate body the inertial rate is faster (Problem 1: 1.875 against 0.625 rad/s).
2. tan θ = I_tω_t/(I₃ω₃).
3. Towards +y: H changes by τ dt, so the spin axis turns towards the torque.
4. The torque turns with the body, so its effect on H averages out each revolution; the wander is about 2τ/(I₃ω₃²), falling as 1/spin².
5. With dissipation only spin about the axis of largest moment is stable. Proof: T = Σ Hₖ²/(2Iₖ) ≥ H²/(2I_max) at fixed H, and dissipation lowers T while H is conserved.
6. A rigid body spinning about an extreme axis is stable (no energy is lost, so it stays on its polhode); in reality the antennas dissipated energy, and the only state of least energy at fixed H is spin about the major axis.
7. Conserved: the total angular momentum Iω + Jω_d. Dissipated: kinetic energy, at c|ω_d − ω|². Final energy: H²/(2(I_t + J)) for a prolate body (spin about the major axis of body plus sphere).
8. Dissipation lowers T at fixed H; for an oblate body T is lowest at θ = 0 (so nutation decays), for a prolate body at θ = 90° (so it grows).
