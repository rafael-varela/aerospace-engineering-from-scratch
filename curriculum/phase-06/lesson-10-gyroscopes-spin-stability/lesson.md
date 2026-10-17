# P06.10  Gyroscopes & Spin Stability

> A spinning body resists being turned and turns sideways when pushed. It keeps its spin axis only if that axis has the largest moment of inertia, because anything that wastes energy will find the axis that needs least.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the nutation of an axisymmetric spinner: the body-frame rate (I₃ − I_t)ω₃/I_t, the inertial coning rate H/I_t and the half-cone angle, and verify all three against simulation. 2. Derive and verify gyroscopic precession: an inertially fixed torque tilts H towards itself at τ/H, and a torque kept perpendicular to H and to a fixed axis turns H at τ/(H sin θ). 3. Show that a torque fixed in a fast-spinning body only makes H wander by about 2τ/(I₃ω₃²), and use it to explain spin stabilisation of an upper stage. 4. Prove the major-axis rule from energy at fixed angular momentum, derive the energy-sink rate dθ/dt, and explain Explorer 1's flat spin. 5. Implement a body with a Kane damper, derive its linear nutation growth or decay rate as the eigenvalue of a 2 × 2 complex matrix, and verify it, the conservation of H and the predicted flat-spin energy against the nonlinear simulation. |
| **Prerequisites** | `P06.09` (Euler's equations, conservation, the simulator), `P06.08` (quaternions), `P06.07` (coning as a 3-1-3 motion), `P04.04` (the heavy top, the gyroscope estimate, the rigid-body stability of extreme axes), `P01.16` (eigenvalues, complex eigenvalues of a real system), `P02.06` (linearisation). Skill check: you can write a pair of real linear ODEs ẋ = −ay, ẏ = ax as one complex equation ż = iaz with z = x + iy. |
| **Engineering question** | A spin-stabilised satellite or upper stage must hold its axis while a thruster or a misaligned motor pushes on it, and it must not tumble on its own. How large is the pointing error a 1 N·m body-fixed torque causes at 1, 2 and 4 rad/s? Which spin axis is safe, and how long does a body with an internal energy sink take to go from a long-axis spin to a flat spin, as Explorer 1 did in 1958? |
| **Assumptions and limits** | Rigid body plus, where stated, one internal damper modelled as a sphere of moment of inertia J in a viscous fluid (a Kane damper) with a linear viscous torque c(ω_d − ω); axisymmetric bodies for the closed forms (I₁ = I₂ = I_t); no external torque except the stated ones; small nutation for the linear rates. Explorer 1's inertia is modelled as a uniform cylinder of its approximate size and mass, which it was not; its whip antennas, the real energy sink, are not modelled (the damper stands in for any energy sink, which is the point of the energy-sink argument). |
| **Analytical work** | Nutation of a symmetric body in the body and in space; the half-cone angle; precession under an inertially fixed torque and under a torque that follows H; spin averaging of a body-fixed torque; T at fixed H as a function of the nutation angle and the major-axis rule; the energy-sink rate; the Kane damper's equations, its conserved H and its dissipation c\|ω_d − ω\|²; the linearised complex 2 × 2 system and its eigenvalues; the flat-spin end state H²/(2(I_t + J)). |
| **Computational deliverable** | `starter/main.py`: `nutation`, `axis_coning_rate`, `precession_rate`, `gyroscope_torque`, `tilt_under_fixed_torque`, `spin_averaged_wander`, `energy_at_fixed_momentum`, `energy_sink_theta_rate`, `spin_stability`, `flat_spin`, `kane_linear_matrix`, `kane_nutation_mode`, `simulate_with_damper`, `damper_energy` (the `P06.08`/`P06.09` helpers, `simulate_rigid`, `inertial_momentum`, `spin_axis_inertial`, `fit_growth_rate`, `explorer_like` and `main` are given). Tests prove: the nutation closed forms against simulation for an oblate and a prolate body (coning rate and body rate to 10⁻⁷, cone angle to 10⁻⁸); gyroscopic precession to 10⁻⁷ with \|H\| and the tilt held; the tilt under a fixed torque; the spin-averaged wander to 2 % and its 1/ω² scaling; the energy at fixed H and the stability rule; the damper's undamped and weak-damping limits; prolate destabilised and oblate stabilised for four damping coefficients; H conserved and energy never rising; the simulated growth rate equal to the linear eigenvalue to 1 %; the flat spin's energy and rate to 10⁻⁶; the oblate nutation decaying at its predicted rate. |
| **Practical activity** | **computational experiment**: spin averaging against spin rate; a damper sweep of the nutation growth rate of an Explorer-1-like spinner; a prolate spinner and an Explorer-like spinner simulated from long-axis spin to flat spin; the energy-sink rate checked along the way; the same damper on an oblate spinner. |
| **Equipment** | None. A computer with the environment from `P00.01`. (A toy gyroscope or a spinning bicycle wheel held by its axle shows precession; nothing is measured with it.) |
| **Safety** | None beyond normal computer use. If you try a bicycle wheel, hold it by the axle, keep fingers clear of the spokes and keep it away from your face. |
| **Evidence** | Completed `main.py`, passing tests, `spin_stability.png`, the derivations of the nutation rates, the spin-averaged wander, the major-axis rule, the energy-sink rate and the damper's linear system, the experiment tables, and a half-page note answering the engineering question. |
| **Assessment** | `assessments/phase-06/grading-rubric.md`; competencies 1–4. Pass: the major-axis rule and the energy-sink rate are derived (not quoted); all public tests pass; the simulated nutation growth is compared with the linear eigenvalue and the flat-spin energy with H²/(2(I_t + J)), with numbers; the note states which axis to spin about and why. |
| **References** | `@mit-ocw-16-07-dynamics` L27 "3D rigid body dynamics: kinetic energy, instability, equations of motion" (its Explorer 1 passage: launched 31 January 1958, spin-stabilised about its minimum-inertia axis, flexible antennas dissipated energy, and it went into spin about its maximum-inertia axis); L28 "Stability of a Gyrostat"; L29 "Description of Free Motions of a Rotating Body Using Euler Angles"; L30 "3D rigid body dynamics: tops and gyroscopes" (Steady Precession: Gyroscopic Motion; Extreme Tops: Gyroscopes) (required-free). `@libretexts-cline-variational-principles` 13.20 Torque-free rotation of an inertially-symmetric rigid rotor; 13.22 Stability of torque-free rotation of an asymmetric body; 13.23 Symmetric rigid rotor subject to torque about a fixed point (required-free). `@mit-ocw-16-61-aerospace-dynamics` lecture 13 "Axisymmetric Rotations"; lecture 14 "Gyroscopes"; lecture 15 "Spacecraft Attitude Dynamics" (supplementary-free). `@openstax-university-physics-vol1` 11.4 Precession of a Gyroscope (required-free). Context: `@nasa-basics-of-space-flight` Chapter 11 "Onboard Systems", Attitude and Articulation Control Subsystems (AACS), the "Spin" and "3-Axis" paragraphs (supplementary-free, qualitative). |
| **Extensions** | ★ Replace the Kane damper by a mass on a spring and dashpot moving along a body axis (a closer model of a flexible antenna) and compare its growth rate with the energy-sink estimate. ★ Model a dual-spin spacecraft (a despun platform on a spinning rotor, the gyrostat of 16.07 L28) and show that a damper on the despun part can stabilise spin about the minor axis of the whole. ★ Fire a body-fixed thruster in short pulses once per revolution and turn the spin axis by a chosen angle (a precession manoeuvre); compare the number of pulses with τΔt/H per pulse. |
| **Time** | unmeasured |

## 1. The engineering problem

The cheapest way to point a spacecraft is to spin it. A spinning body has angular momentum H along its spin axis, and in space nothing changes H except external torque, so the axis stays where it was put; small disturbing torques change it only slowly and, as this lesson shows, mostly cancel themselves out over each revolution. Early satellites, many communications satellites, interplanetary probes and solid-motor upper stages have all been spin-stabilised (Basics of Space Flight, AACS, contrasts spin and three-axis stabilisation). The same physics makes a gyroscope precess instead of falling, which `P04.04` met for a heavy top.

There is a trap. The first American satellite, Explorer 1, launched on 31 January 1958, was a long thin cylinder spun about its long axis, the axis of *minimum* moment of inertia. Rigid-body theory says that is stable (`P04.04`, `P06.09`). But its flexible whip antennas flexed and dissipated a little energy, and the satellite ended up spinning end over end about an axis of *maximum* moment of inertia, a flat spin (16.07 L27 tells the story). Nothing external had acted. The cause is a theorem about energy at fixed angular momentum, the **major-axis rule**, and it has governed spacecraft design since.

This lesson builds the pieces with closed forms and checks each against simulation: the nutation (coning) of a symmetric spinner; precession under a torque; spin averaging of a torque fixed in the body, the reason a spinning upper stage keeps its axis while its motor's thrust is slightly misaligned; and the passage from minor-axis spin to flat spin when the body carries an energy sink, modelled as a **Kane damper**, a heavy sphere in a viscous fluid inside a spherical cavity. Success criteria, written first: every closed form matched by simulation to the solver tolerance or to the stated linearisation error; the total angular momentum conserved to 10⁻⁸ with the damper; the end state's energy equal to H²/(2(I_t + J)).

## 2. Understand and derive

### 2.1 Nutation of a symmetric spinner

Take an axisymmetric body, I₁ = I₂ = I_t (transverse) and I₃ (spin axis), with no torque. Euler's third equation gives ω̇₃ = 0. The first two, with z = ω₁ + iω₂ (Cline 13.20):

I_t ż = i(I₃ − I_t)ω₃ z  ⇒  **z(t) = ω_t e^{iλt},  λ = (I₃ − I_t)ω₃/I_t.**   (1)

The transverse rate keeps its size ω_t and turns in the body at λ, with the spin for an oblate body (I₃ > I_t), against it for a prolate one. In space the picture is different. H = (I_tω₁, I_tω₂, I₃ω₃) has fixed size H = √((I₃ω₃)² + (I_tω_t)²) and fixed direction; the symmetry axis e₃ makes a constant angle θ with it,

**tan θ = I_t ω_t / (I₃ ω₃),**   (2)

and turns about H on a cone (the **space cone**) at the **inertial nutation rate**

**Ω_n = H / I_t.**   (3)

To see (3), write H = I_tω + (I₃ − I_t)ω₃e₃, so ω = H/I_t − ((I₃ − I_t)ω₃/I_t) e₃. The second term spins the body about its own axis and does not move e₃; the first turns everything, including e₃, about the fixed H at H/I_t. This is the coning of `P06.07` with ψ̇ = H/I_t, β = θ and φ̇ = −λ (16.07 L29 describes the same free motion with Euler angles). For the lesson's oblate satellite (I_t = 800, I₃ = 1200 kg·m², ω₃ = 1 rad/s, ω_t = 0.05 rad/s): H = 1200.6665 N·m·s, λ = 0.5000 rad/s, Ω_n = **1.500833 rad/s**, θ = **1.9092°**. The simulation's axis cones at 1.500833 rad/s at a constant 1.9092°.

### 2.2 Precession under a torque

dH/dt = τ in the inertial frame. A body spinning fast has H nearly along its spin axis, so the axis goes where H goes: **towards the torque, not towards the push that makes it**.

- **Inertially fixed torque perpendicular to H₀.** H(t) = H₀ + τt: the tilt grows as atan(τt/H), in the direction of τ. A 5 N·m torque on the 1200 N·m·s satellite tilts it 9.5° in 40 s (a test checks the whole history).
- **Torque kept perpendicular to H and to a fixed axis ẑ**, like gravity's torque on a top: |H| and the tilt θ from ẑ cannot change, and H turns about ẑ at

**Ω = τ / (H sin θ).**   (4)

At θ = 90° this is OpenStax 11.4's Ω = τ/L; 16.07 L30 derives it as steady gyroscopic motion, and `P04.04` found the same formula as the fast-top limit. With τ = 2 N·m, H = 1200 N·m·s and θ = 30°: Ω = **0.003333 rad/s**, simulated **0.003333 rad/s** with the tilt held at 30.000000° and |H| to 9 × 10⁻¹² N·m·s.

### 2.3 Spin averaging

A torque fixed in the body, such as a thrust misaligned by a few millimetres from the centre of mass of a spinning upper stage during its motor burn, does not stay fixed in space: it turns with the body at about the spin rate ω₃. Integrate dH/dt = τ with τ rotating at ω₃ in the plane perpendicular to H: the tip of H runs round a circle of radius τ/ω₃ (N·m·s) that passes through its starting point. The direction of H wanders by at most

**Δθ_max ≈ 2τ / (H ω₃) = 2τ / (I₃ ω₃²),**   (5)

and returns every revolution. A non-spinning body under the same torque would turn without limit. The spin converts a secular drift into a bounded wobble, and doubling the spin quarters it. For a 1 N·m body-fixed torque on the 1200 kg·m² (axial) satellite: 0.09549°, 0.02387° and 0.00597° at 1, 2 and 4 rad/s, by (5) and by simulation alike; the same torque would turn a non-spinning 800 kg·m² body through 129° in 60 s.

### 2.4 The major-axis rule

At fixed angular momentum, the kinetic energy of a symmetric body whose axis is at θ from H is (Problem 2)

**T(θ) = (H²/2) (sin²θ / I_t + cos²θ / I₃).**   (6)

It runs from H²/(2I₃) at θ = 0 (pure spin) to H²/(2I_t) at θ = 90° (a **flat spin** about a transverse axis). For any body, T ≥ H²/(2I_max): of all motions with the same H, spin about the axis of largest moment has the least energy. Internal dissipation (flexing antennas, sloshing fuel, a damper) cannot change H, which only external torques can, but it removes energy. So it drives T down to H²/(2I_max). Hence **the major-axis rule: with any energy dissipation, spin about the axis of maximum moment of inertia is the only stable spin** (16.07 L27–L28; Cline 13.22). For an oblate body (I₃ > I_t) that is the spin axis: dissipation damps nutation. For a prolate body (I₃ < I_t), like Explorer 1, dissipation makes nutation grow until the body tumbles end over end.

**The energy-sink rate.** If dissipation removes power P, then dT/dt = −P with T from (6):

**dθ/dt = −P / (H² sin θ cos θ (1/I_t − 1/I₃)),**   (7)

positive for a prolate body. Equation (7) is a statement about energy and momentum only; P can come from any mechanism. To use it one needs P, which requires a model of the sink.

**Explorer 1, illustrated.** As a uniform cylinder of 13.97 kg, 2.03 m long and 0.152 m across (approximate size and mass; the real satellite was not uniform): I_t = 4.8176 kg·m², I₃ = 0.04035 kg·m². Spun at 750 rpm (78.540 rad/s), H = 3.1687 N·m·s and T₀ = 124.44 J. The flat spin at the same H turns at H/I_t = 0.6577 rad/s (6.28 rpm) with T = 1.0421 J: dissipating 99.2 % of the energy is all it takes. The rigid model calls the long-axis spin `stable`; with dissipation it is `unstable`.

### 2.5 A concrete energy sink: the Kane damper

Put a sphere of moment of inertia J (about any axis) in a spherical cavity filled with viscous fluid at the body's centre. Let ω_d be the sphere's absolute angular velocity in body axes. The fluid exerts a torque c(ω_d − ω) on the body and −c(ω_d − ω) on the sphere. Writing dH/dt = τ for each in body axes (the sphere's H is Jω_d because it is spherical):

**I ω̇ + ω × (I ω) = c(ω_d − ω),   J ω̇_d + ω × (J ω_d) = −c(ω_d − ω).**   (8)

Adding them, the total H = Iω + Jω_d obeys dH/dt + ω × H = 0, so it is conserved in inertial axes; and the energy T = ½ω·Iω + ½J|ω_d|² falls at exactly c|ω_d − ω|² ≥ 0 (Problem 2). When nothing moves relative to the body, ω_d = ω and nothing dissipates: the only motions without loss are rigid rotations of body and sphere together, which by 2.4 end at spin about the major axis of the combined inertia, with **T_final = H²/(2(I_t + J))** for a prolate body.

**Linear growth rate.** Linearise (8) about spin Ω (body and sphere co-rotating) with small transverse rates z = ω₁ + iω₂ and z_d = ω_d1 + iω_d2:

I_t ż = i(I₃ − I_t)Ω z + c(z_d − z),   J ż_d = −c(z_d − z) − iJΩ(z_d − z).   (9)

This is ẋ = Ax with a complex 2 × 2 matrix A. With c = 0 its eigenvalues are i(I₃ − I_t)Ω/I_t, the body's nutation (1), and −iΩ, the sphere staying still in space while the body turns around it. With c > 0 the nutation eigenvalue moves off the imaginary axis: its real part is the nutation growth rate (positive) or decay rate (negative). For weak damping it is proportional to c (a test); for very strong damping the sphere is locked to the body and dissipation stops again, so there is an optimum c.

### 2.6 Where the model fails

- **The sink.** A real antenna or fuel tank is not a Kane damper; its P depends on frequencies and modes. The direction of the drift (7) does not; its speed does.
- **Uniform cylinder.** Explorer 1's mass was not uniform, so its real inertia ratio differed; the conclusion (prolate, so unstable with dissipation) does not depend on the details.
- **No external torque.** Gravity gradient, magnetic and aerodynamic torques act on a real satellite at low altitude over hours; they change H, which this lesson's closed forms hold fixed.
- **Small nutation in (1)–(5) and (9).** The nonlinear simulation is not limited to it; (7) is exact for a symmetric rigid body, given P.

## 3. Build the computational model

Open `starter/main.py`. The quaternion helpers, `simulate_rigid` (the `P06.09` simulator), `inertial_momentum`, `spin_axis_inertial`, `fit_growth_rate`, `explorer_like` and `main` are given.

### Step 1: closed forms

`nutation` returns H, λ, H/I_t and θ from (1)–(3). `axis_coning_rate` measures the coning rate from a simulation: the unwrapped azimuth of e₃ in the plane perpendicular to H, fitted by a straight line. `precession_rate` is (4), `tilt_under_fixed_torque` atan(τt/H), `spin_averaged_wander` (5), `energy_at_fixed_momentum` (6), `energy_sink_theta_rate` (7), `flat_spin` H/I_max and H²/(2I_max), `spin_stability` the rule of 2.4 for a rigid or a dissipative body.

### Step 2: the gyroscope torque

`gyroscope_torque(I, tau0)` returns a torque function of (t, q, ω): it computes H_i = Cᵀ(Iω), takes τ_i = τ₀ (ẑ × H_i)/|ẑ × H_i| and returns Cτ_i in body axes.

### Step 3: the damper

```python
def simulate_with_damper(I, J, c, w0, wd0, t_end, ...):
    def rhs(t, y):
        q, w, wd = y[:4], y[4:7], y[7:]
        rel = wd - w
        return np.concatenate([0.5 * quat_multiply(q, [0, *w]),
                               Iinv @ (c * rel - cross3(w, I @ w)),
                               -c * rel / J - cross3(w, wd)])
```

`damper_energy` is the total kinetic energy; `kane_linear_matrix` is A of (9) and `kane_nutation_mode` picks the eigenvalue nearest the undamped nutation i(I₃ − I_t)Ω/I_t.

Running `main()` prints the nutation check, the gyroscope, the spin-averaging table, the Explorer-like closed forms and damper sweep, the prolate spinner simulated to its flat spin with the energy-sink check, and the oblate damper, and writes `spin_stability.png`: the oblate spinner's constant cone angle; the prolate spinner's transverse rate against the linear prediction; and its nutation angle and energy on the way to the flat spin.

### Verify

Run `python -m pytest curriculum/phase-06/lesson-10-gyroscopes-spin-stability`. The tests:

1. **Nutation**: the closed forms, including the prolate sign; simulated oblate and prolate spinners with H fixed to 10⁻⁹, the cone angle constant at (2), the coning rate equal to (3) and the body rate to (1) to 10⁻⁷.
2. **Precession**: the gyroscope at (4) to 10⁻⁷ with |H| and the tilt fixed; τ/L at 90°; an inertially fixed torque tilting H along itself as atan(τt/H).
3. **Spin averaging**: the largest wander within 2 % of (5) at two spins and a ratio of one quarter.
4. **Energy and stability**: (6) at its ends and its monotonic behaviour for prolate and oblate bodies; the signs of (7); the flat spin; the stability rule for a box, the Explorer-like cylinder and the oblate satellite.
5. **Damper, linear**: the undamped eigenvalues; growth proportional to c for weak damping; prolate growing and oblate decaying for four values of c.
6. **Damper, nonlinear**: H conserved to 10⁻⁸ and energy never rising; no loss with c = 0; the prolate spinner growing at the eigenvalue's real part to 1 % and ending in a flat spin with T = H²/(2(I_t + J)) and rate H/(I_t + J) to 10⁻⁶; the oblate spinner's nutation decaying at its predicted rate to 3 %.
7. `main` writes the figure.

### Compare

The independent checks are the closed forms against the simulation and the energy-sink prediction against the measured motion. For the prolate spinner (I_t = 1500, I₃ = 600 kg·m², 3 rad/s, J = 60 kg·m², c = 100 N·m·s) the linear theory gives a nutation mode of 0.03202 − 1.7546i s⁻¹; the simulation's transverse rate grows at **0.03194 s⁻¹** while it is small. Over 600 s the total angular momentum, direction and size, drifts by **5.9 × 10⁻¹⁰**, the energy never rises, the axis passes 89° from H at 341.6 s, and the run ends with the axis at 89.9999° from H, the body turning at **1.269231 rad/s** against H/(I_t + J) = 1.269231 rad/s, and **T = 1256.5392 J** against H²/(2(I_t + J)) = 1256.5392 J (from T₀ = 2970.0 J). The energy-sink rate (7), with the measured dissipation averaged over 10 s and the co-rotating sphere counted in the inertias (I_t + J, I₃ + J), predicts dθ/dt within 0.4–2.1 % of the simulated rate at nutation angles of 5°, 20°, 45° and 70° (table in part 5). The usual mistakes show at once: the damper torque applied with the same sign to body and sphere breaks H conservation; leaving out ω × Jω_d in the sphere's equation does too; and the wrong eigenvalue (the sphere's) gives a "growth rate" with the opposite sign.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** All results are outputs of the model, not measurements. Write each prediction down before running it.

1. **Spin averaging.** A 1 N·m torque fixed along body axis 1 of the oblate satellite, spins of 1, 2 and 4 rad/s, 60 s. Predict the largest wander of H from (5). *Finding:* 0.09549°, 0.02387°, 0.00597°, equal to (5) to the digits shown. A spinning upper stage at 4 rad/s turns a large, unknown misalignment torque into a pointing error of a few thousandths of a degree.
2. **Damper sweep for the Explorer-like spinner.** At 750 rpm, compute the nutation mode for (J, c) = (0.004, 10⁻⁴), (0.004, 3 × 10⁻³), (0.02, 10⁻²) (kg·m², N·m·s). Predict which grows fastest. *Finding:* growth rates 2.45 × 10⁻³, 3.08 × 10⁻², 1.26 × 10⁻¹ s⁻¹ (e-folding in 408 s, 32.5 s and 7.9 s, i.e. 5095, 406 and 99 revolutions); nutation frequency −77.88, −77.85, −77.78 rad/s, close to (I₃ − I_t)Ω/I_t = −77.88 rad/s. Even a strong sink needs about a hundred revolutions per e-fold, because the long thin body's nutation is fast and the damper is poorly coupled to it. The real antennas' rate is not computed here.
3. **Explorer-like spinner to the flat spin.** Simulate (J, c) = (0.02, 0.01) from a transverse rate of 0.05 rad/s for 120 s (rtol 10⁻⁹; about 1.6 s of computer time). Predict the end state. *Finding:* the total H is 4.7456 N·m·s (the sphere, spinning with the body at 78.5 rad/s, carries 1.57 N·m·s of it) and is conserved to 9.7 × 10⁻⁹; energy never rises; the axis passes 10°, 45°, 80° and 89° from H at 9.1, 21.9, 31.3 and 42.2 s; the fitted early growth rate is 0.1246 s⁻¹ against 0.1263 s⁻¹ from (9); the run ends in a flat spin at 0.98099 rad/s with T = 2.32771 J against H²/(2(I_t + J)) = 2.32771 J, from 186.13 J.
4. **The prolate spinner and the energy sink.** The case of the Compare paragraph. Predict dθ/dt from (7) at 5°, 20°, 45° and 70° using the dissipated power measured in the simulation. *Finding:* see the table in part 5.
5. **The same damper on an oblate spinner.** I_t = 800, I₃ = 1200 kg·m², 1 rad/s, J = 60, c = 100. Predict the sign of the nutation rate. *Finding:* −0.01732 + 0.4792i s⁻¹: the nutation decays with a 58 s time constant, and a simulation from ω_t = 0.05 rad/s ends after 300 s below 1 % of it. This is why oblate spinners carry nutation dampers, and why prolate ones must not.

Raw data format: CSV tables `spin_radps, wander_sim_deg, wander_closed_deg` (1); `J, c, growth_per_s, freq_radps, efold_s, efold_revs` (2); `t_s, nutation_deg, T_J, H_drift` (3, 4); `nutation_deg, P_W, dtheta_pred, dtheta_sim` (4); `case, eigenvalue` (5); plus `spin_stability.png`.

## 5. Validate and improve

| quantity | closed form | simulation | difference | cause |
|---|---|---|---|---|
| oblate coning rate (rad/s) | H/I_t = 1.500833 | 1.500833 | < 10⁻⁷ relative | — |
| oblate half-cone angle (°) | 1.9092 (2) | 1.9092 (constant) | < 10⁻⁸ rad | — |
| gyroscope precession (rad/s) | τ/(H sin θ) = 0.003333 | 0.003333 | < 10⁻⁷ relative | — |
| spin-averaged wander, 1 rad/s (°) | 2τ/(I₃ω₃²) = 0.09549 | 0.09549 | < 2 % (tested) | nutation, small-angle algebra |
| prolate growth rate (s⁻¹) | 0.03202 (9) | 0.03194 | 0.3 % | finite amplitude, fit window |
| flat-spin energy (J) | H²/(2(I_t + J)) = 1256.5392 | 1256.5392 | < 10⁻⁷ relative | — |
| dθ/dt at 5° (rad/s) | (7): 0.00278 | 0.00279 | −0.4 % | averaging over 10 s |
| dθ/dt at 20° | 0.01031 | 0.01039 | −0.8 % | sphere not exactly co-rotating |
| dθ/dt at 45° | 0.01644 | 0.01654 | −0.7 % | same |
| dθ/dt at 70° | 0.01121 | 0.01098 | +2.1 % | same, faster change |
| oblate damper (s⁻¹) | −0.01732 (decays) | decays below 1 % in 300 s; fit within 3 % | — | — |

Every closed form is reproduced; the energy-sink rate, which treats the sphere as part of the body, is right to about 2 %, and better where the nutation changes slowly. The engineering answer: **a 1 N·m body-fixed torque wanders the axis by 0.095° at 1 rad/s, 0.024° at 2 rad/s and 0.006° at 4 rad/s; spin about the axis of maximum moment of inertia, because any energy sink drives a minor-axis spinner to a flat spin, as it drove Explorer 1. With the illustrative damper of Experiment 3, the Explorer-like cylinder goes from long-axis spin to 89° of nutation in 42 s; with a weak sink the same path takes thousands of revolutions per e-fold, but it is the same path.**

What to improve: the sink model (★ a spring–mass damper or a flexible-antenna model would give a P with the right frequency dependence); and for a minor-axis design that must exist (a long launcher stage), an active nutation controller or a dual-spin configuration (16.07 L28, "Stability of a Gyrostat").

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (1)–(7), the conservation and dissipation properties of (8), and the linear system (9); `spin_stability.png`; the five experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) answering the engineering question.

### Adversarial self-checks

- **What if the spin doubles?** The spin-averaged wander falls fourfold (5); the inertial nutation rate doubles; for the damper, the growth rate per revolution stays the same if c scales with the spin, so the path to the flat spin takes the same number of revolutions in half the time.
- **What if I₃ = I_t?** The body is a "sphere" for nutation: λ = 0, and T at fixed H does not depend on θ, so dissipation has nothing to drive; (7) divides by zero, correctly signalling that the rule gives no preference.
- **Which assumption fails first in the real world?** For an upper stage, rigidity (its burning propellant moves and sloshes); for a long-lived satellite, the absence of external torque.
- **What measurement could falsify the model?** Telemetry of a spinner's nutation angle against time: with a known dissipation mechanism, (7) predicts dθ/dt from the measured P; growth on an oblate body would point at an external torque or a wrong inertia.
- **Does the result survive a change of solver?** The conserved H and the monotonic energy are the test; tightening rtol from 10⁻⁹ to 10⁻¹⁰ changes the Explorer-like end state by less than its printed digits. The closed forms need no solver.
- **Which decision changes if the uncertainty doubles?** If the moments were uncertain enough that I₃ and I_t could swap order (a nearly spherical body), the stable axis would be unknown; the design would need a clear margin between the moments, or active control.

## Check yourself

1. Write the body-frame and inertial nutation rates of an axisymmetric spinner. Which is faster for an oblate body?
2. What is the half-cone angle between the spin axis and H, in terms of the rates?
3. A torque acts on a fast spinner along +y while H is along +z. Which way does the spin axis move?
4. Why does a torque fixed in a spinning body cause only a bounded pointing error? How does it scale with spin?
5. State the major-axis rule and prove it in two lines from T and H.
6. Explorer 1 was spun about its long axis. Why was that stable in rigid-body theory and unstable in reality?
7. What is conserved and what is dissipated in a body with a Kane damper? What is the final kinetic energy?
8. Why does a damper stabilise an oblate spinner and destabilise a prolate one?
