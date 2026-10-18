# P10.20  Ground Effect & Interference

> The ground is a streamline, so put a mirror-image wing under it. One reflection explains why aircraft float in the flare, why a tail feels the ground first, and why wing-in-ground-effect craft exist.

| | |
|---|---|
| **Learning outcomes** | 1. Explain with the method of images why a wing near the ground makes more lift at a given angle and less induced drag at a given lift, and derive the closed-form downwash reduction of a single horseshoe vortex. 2. Implement and test an image system in a vortex-lattice model: the influence matrix with mirrored horseshoes, the Trefftz-plane induced drag with an image wake, and the downwash a wing induces at its tail, with and without the ground. 3. Verify the model three ways: the ground-plane normal velocity vanishes, free-air results return as h → ∞, and a single horseshoe reproduces its closed form; compare the induced-drag ratio with the classical empirical estimate. 4. Quantify ground effect for the reference trainer in the flare, a large airliner and a wing-in-ground-effect craft, and the wing–tail downwash interference, and judge where the flat, steady image model fails. |
| **Prerequisites** | `P10.09` (horseshoe vortex, Biot–Savart, Trefftz-plane induced drag), `P10.11` (the vortex lattice method), `P10.10` (lifting line, span efficiency), `P01.11` (vector fields, circulation), `P03.08` (direct linear solvers), `P03.11` (observed order). Skill check: you can write the velocity of a 2-D point vortex at a point and add two of them. |
| **Engineering question** | When the reference trainer's wing is 1.0 m above the runway in the flare (h/b = 0.1), by how much do its lift-curve slope and its induced drag change, how much does the downwash at its tail change, and what does that do to the pitching moment? The answers from a 4 × 40 vortex lattice with images: CL_α rises 16.6 %; at the same lift the induced drag falls 46.9 % (K ratio 0.531; the classical empirical estimate 0.515), from 314 N to 167 N at 30 m/s (4 × 40 lattice; about 323 N to 173 N when converged); the downwash gradient at the tail falls from 0.318 to 0.124, so at 12° the tail sees 2.3° less downwash and the aircraft pitches nose-down by ΔCm ≈ −0.095. |
| **Assumptions and limits** | Inviscid, incompressible, steady flow; a thin flat lifting surface (camber appears only through the zero-lift angle, thickness not at all); a flat-planar wake trailing straight downstream (no roll-up, no wake deflection); a flat ground stationary relative to the air far ahead; small angles. It misses the venturi suction under thick or flapped wings very near the ground, viscous effects of a ground boundary layer, the dynamics of a descending wing, waves, and the fuselage and propeller slipstream. |
| **Analytical work** | The image construction and its proof for 2-D vortices and 3-D filaments; the single-horseshoe Trefftz downwash ratio (4h/b_v)²/(1 + (4h/b_v)²); why lift at fixed α rises; the Trefftz-plane drag with an image wake; the far-field downwash gradient 2CL_α/(πAR) and why a rectangular wing falls below it; the tail-moment change ΔCm = −V_H a_t Δε. |
| **Computational deliverable** | `starter/main.py`: `mirror`, `influence_matrix`, `solve_lattice`, `strip_circulation`, `lift_coefficient`, `trefftz_downwash`, `induced_drag_coefficient`, `wing_in_ground_effect`, `ground_plane_residual`, `single_horseshoe_trefftz_ratio`, `wieselsberger_ratio`, `tail_downwash`, `downwash_gradient_far_field` (the Biot–Savart helpers `segment_velocity`, `semi_infinite_velocity`, `horseshoe_velocity`, the `Lattice` class, `rectangular_lattice`, `observed_order` and `main` are given). Tests (14) prove: mirroring is an involution; the image system leaves a ground-plane residual below 10⁻¹²; h → ∞ recovers free-air CL, CDi and span load; the span load is symmetric; lift is linear in sin α; a very long wing tends to 2π; lift rises and the induced-drag factor falls monotonically as h falls; a single horseshoe matches its closed form; the Trefftz downwash equals the 3-D Biot–Savart sum far downstream; the lattice agrees with the empirical ratio within 0.05; the tail downwash tends to the wake value and falls with tail height and in ground effect; the lattice converges monotonically at a positive order. |
| **Practical activity** | **computational experiment**: induced-drag ratio and lift change against h/b compared with the single-horseshoe closed form and the empirical estimate; lattice convergence; the trainer in the flare; a wing-in-ground-effect craft; downwash at the tail against distance, height and ground proximity. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `ground_effect.png`, the image proof and the single-horseshoe derivation, the five experiment tables, and a half-page note on the trainer's flare: lift, drag, tail moment and the assumption most likely to fail. |
| **Assessment** | `assessments/phase-10/grading-rubric.md`; competencies 1–4. Pass: the image condition is proved and checked numerically; all public tests pass; the lattice's ground-effect ratio is compared with the closed form and the empirical estimate with the differences explained; the note states the flare numbers and their limits. |
| **References** | `@mit-ocw-16-100-aerodynamics` Lecture 33 "Ground Effect Using Single Vortex Model"; Lecture 34 "Three-Dimensional Wall Effects" (https://ocw.mit.edu/courses/16-100-aerodynamics-fall-2005/resources/16100lectre33_cj/) (required-free): the image horseshoe with w = 0 at z = 0, and the same idea for tunnel walls. `@ntrs-naca-tm-77-wieselsberger-ground` Wieselsberger, NACA TM-77 *Wing Resistance Near the Ground* (supplementary-free; read the page images): the image method with Prandtl's wing theory, the origin of the classical estimate. `@ntrs-naca-tn-705-recant-ground-effect-flaps` Recant, NACA TN-705: Introduction; Results and Discussion; Conclusions (supplementary-free): measured ground effect on flapped wings. `@leishman-intro-aerospace-flight-vehicles` "Ground Effect Vehicles" (https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/hovercraft/): Anatomy of GEVs; WIGE Vehicle Performance (required-free). `@ntrs-naca-tr-116-prandtl-modern-hydrodynamics` Prandtl, NACA TR-116 Part II, E. "Improved Theory of Airplanes Having More Than One Wing" (supplementary-free): interference between lifting surfaces. Optional textbooks for deeper study (paid, never required): Anderson, *Fundamentals of Aerodynamics*, and Katz & Plotkin, *Low-Speed Aerodynamics*, for lifting-surface methods and images; library or used copies are fine. |
| **Extensions** | ★ Replace the flat lattice by a cambered one with a deflected flap and compare the lift change with Recant's flapped-wing data (TN-705). ★ Add a horizontal tail as a second lattice in the same influence matrix (with images) and compute the trimmed elevator angle in and out of ground effect. ★ Let the wing descend at constant rate with the wake laid down at the heights where it was shed (a time-marching lattice) and compare the cushion with the steady result. |
| **Time** | unmeasured |

## 1. The engineering problem

Every pilot learns that an aircraft "floats" in the last metres before touchdown: it keeps flying when the speed should have let it settle, it wants less nose-up attitude for the same lift, and, on many types, it pitches nose-down just as the pilot is trying to raise the nose. All three are ground effect. Designers exploit the same physics deliberately: wing-in-ground-effect (WIG) craft skim a few metres above water to fly with half the induced drag of an aeroplane, and the Caspian Sea ekranoplans were built around it. And the downwash a wing throws back onto its own tail, which sets the aircraft's stability (`P17.03`), is the same kind of interference between vortex systems.

The tool is the vortex lattice of `P10.11` plus one idea. The ground is a plane the air cannot cross, so the vertical velocity there must vanish. Placing a mirror-image wing at the same distance below the ground, with the opposite circulation, enforces that condition exactly and everywhere. The wing then sees the image's induced flow as an upwash: it needs more circulation to meet the tangency condition at a given angle (more lift), and its wake has less downwash for a given circulation (less induced drag).

Success criteria, written first: the image system leaves a normal velocity on the ground plane at round-off level; the free-air result returns as h/b → ∞; the lattice reproduces the closed form for a single horseshoe; its induced-drag ratio agrees with the classical empirical estimate within that estimate's scatter (a few hundredths); the lattice converges; and the trainer's flare is answered with lift, drag and pitching-moment changes.

## 2. Understand and derive

### 2.1 The image condition

A 2-D point vortex of strength Γ along +x at (y₀, z₀) induces, at (y, z),

w = Γ(y − y₀)/(2πr²),  v = −Γ(z − z₀)/(2πr²),  r² = (y − y₀)² + (z − z₀)².

Put one vortex at height h and another of strength −Γ at −h. On the plane z = 0 both have the same r and the same y − y₀, so their w cancel exactly; their v add, which is allowed, because the ground only forbids flow *through* it. In three dimensions the same holds element by element for Biot–Savart: reflect every point of a filament in the plane and reverse its circulation, and on the plane the image's velocity is the mirror image of the original's, with the opposite vertical component (the proof with the reflection matrix is in the worked solutions, Problem 2). So the ground is modelled by a mirror-image lattice at z → −z with Γ → −Γ. The lesson's 16.100 reading (Lecture 33, "Ground Effect Using Single Vortex Model") builds exactly this for one horseshoe; Wieselsberger did it with Prandtl's wing theory in 1922 (`@ntrs-naca-tm-77-wieselsberger-ground`).

![Front view, looking upstream, of a wing at height h above the ground plane z = 0, its trailing vortices of circulation Γ seen end-on at the tips, and its mirror image at −h with circulation −Γ and the opposite sense of rotation. At a point on the ground the velocities induced by a tip vortex and its image have equal and opposite vertical parts, so their sum lies along the ground.](figures/ground-image-vortices.svg)
*Figure 1. The image system. Reflecting the lattice in z = 0 and reversing Γ cancels the normal velocity everywhere on the plane, so the ground needs no boundary condition of its own.*

### 2.2 One horseshoe, in closed form

A horseshoe of span b_v sheds two trailing vortices. In the Trefftz plane far downstream (`P10.09`) the downwash at the centre from the real pair is w₀ = 2Γ/(π b_v). The image pair sits 2h below, at distance r = √((b_v/2)² + (2h)²), and induces an upwash whose vertical part is the fraction (b_v/2)/r of Γ/(2πr) each. The ratio is

**w/w₀ = (4h/b_v)² / (1 + (4h/b_v)²).**   (1)

At h/b_v = 0.25 the downwash is halved; at h/b_v = 1 it is 94 % of free air. With b_v = πb/4 (the span of the horseshoe equivalent to an elliptic wing), the trainer's h/b = 0.1 gives 0.206. That number is far too pessimistic, as part 4 shows: one horseshoe concentrates all the vorticity at two points directly above their images. It is the code's closed-form test, not the answer.

### 2.3 What the image does to lift and drag

**Lift at a fixed angle.** The tangency condition at each control point i is Σⱼ w_ij Γⱼ = −V sin α, with w_ij the normal velocity at i from a unit horseshoe j. The images subtract from each w_ij (they induce upwash), so the circulations needed grow: **CL at fixed α rises**. At fixed circulation the lift would not change; what the ground changes is the circulation a given angle produces.

**Induced drag at a fixed lift.** In the Trefftz plane, Di = (ρ/2)∫Γ(y) w(y) dy with w the far-wake downwash. With the ground, w includes the image wake's upwash and is smaller for the same Γ(y). At the same lift, **the induced drag falls**. The lesson measures it with the factor K = CDi/CL² relative to its free-air value, which is what matters at a fixed weight and speed.

**The classical estimate.** A widely quoted empirical fit, attributed to Wieselsberger's work, gives the induced-drag ratio as

**K/K_free ≈ 1 − (1 − 1.32 h/b)/(1.05 + 7.4 h/b),** usually stated for 0.033 < h/b < 0.25.   (2)

It is a fit, not a derivation, and it is used here only as a comparison.

### 2.4 Downwash at the tail

A wing's trailing vortices induce a downwash ε behind it. For an elliptic wing the downwash angle at the wing is CL/(πAR), and far behind it is twice that, uniform across the span (`P10.09`), so the classical estimate of the downwash gradient at a tail is

**dε/dα ≈ 2CL_α/(πAR).**   (3)

The tail's angle of attack is α − ε, so dε/dα reduces the tail's contribution to stability (`P17.03`). The lattice computes ε directly by Biot–Savart at points across the tail span. Two things make the trainer differ from (3): the tail is close behind the wing (4.5 m, under four chords), where the downwash has not reached its far value; and a rectangular wing sheds its vorticity strongly at the tips, so the downwash at the centre of its wake, where the tail sits, is below the span average that (3) assumes. In ground effect the images reduce ε at the tail, the tail's angle of attack rises, and the tail lifts more: a **nose-down** moment

**ΔCm ≈ −V_H a_t Δα_t,  Δα_t = −Δε,  V_H = S_t l_t/(S c).**   (4)

### 2.5 Where the model stops

- **Thickness and flaps.** A thick or flapped wing very close to the ground forms a channel with it; the flow accelerates under the wing and can *reduce* the lift (Recant's flapped-wing tests, TN-705, are the kind of data to compare with). A flat lattice cannot see this.
- **The ground boundary layer.** In a wind tunnel a fixed ground board grows its own boundary layer; tunnels use moving belts. In flight the runway is at rest relative to the undisturbed air, as the model assumes.
- **Unsteadiness.** A descending wing's wake was shed higher up; the steady result is an upper bound on the cushion during the descent.
- **Wake shape.** The model's wake is flat and does not roll up; near the ground real tip vortices spread outward along it.

## 3. Build the computational model

Open `starter/main.py`. The Biot–Savart functions, the `Lattice` dataclass with its rectangular builder, `observed_order` and `main` are given (copied in compact form from `P10.11`). The lattice has nc chordwise and ns spanwise panels, bound vortices at panel quarter chords, control points at three-quarter chords, trailing legs to +∞ along x, and cosine spanwise spacing.

### Step 1: images in the influence matrix

`mirror(points)` reflects z → −z. `influence_matrix(lat, ground)` is the z-velocity at every control point from every unit horseshoe, minus the same from the mirrored horseshoes when `ground` is True:

```python
V = horseshoe_velocity(lat.cp, lat.A, lat.B)
if ground:
    V = V - horseshoe_velocity(lat.cp, mirror(lat.A), mirror(lat.B))
return V[..., 2]
```

`solve_lattice` solves for the circulations with right-hand side −V sin α; `strip_circulation` sums each spanwise strip; `lift_coefficient` is CL = 2ΣΓΔy/(VS).

### Step 2: the Trefftz plane with an image wake

`trefftz_downwash(y_edges, G, y_eval, h)`: the trailing vortex at edge k has strength s_k = G_{k−1} − G_k (zero beyond the tips) along +x; add the 2-D velocities of §2.1 from the real vortices on the wake line and, when h is given, from their images 2h below with opposite sign. `induced_drag_coefficient` is CDi = Σ G w Δy/(V²S). `wing_in_ground_effect` puts it together and returns CL, CDi, e and the span load.

### Step 3: checks and the tail

`ground_plane_residual` evaluates w on z = 0 from wing plus image, relative to the wing alone. `single_horseshoe_trefftz_ratio` is (1), `wieselsberger_ratio` is (2), `tail_downwash` averages −w/V across the tail span (with images if asked), `downwash_gradient_far_field` is (3).

### Verify

`main()` prints, among others,

```
free air: CL = 0.4088 (CL_alpha = 4.6850 /rad), CDi = 0.00642, e = 0.9945
ground-plane residual max|w|/max|w_wing| at 200 random points: 1.33e-15
trainer in the flare, wing at h = 1.0 m (h/b = 0.1): CL_alpha up 16.6 %, induced-drag factor K ratio 0.531 (46.9 % less induced drag at the same CL)
```

The tests (14) check everything listed in the metadata. Two are the independent checks of the image construction: the residual on the ground plane (10⁻¹⁵ relative), and the Trefftz-plane formula against the full 3-D Biot–Savart sum, wing plus image, evaluated 10⁶ m downstream (agreement to 10⁻⁶ relative). A third compares a single strip with the closed form (1) to 10⁻¹².

**Lattice convergence (free air, nc = 4, α = 5°).**

| ns | CL | CDi |
|---|---|---|
| 10 | 0.422075 | 0.0063607 |
| 20 | 0.413589 | 0.0064098 |
| 40 | 0.408845 | 0.0064200 |
| 80 | 0.406353 | 0.0064211 |

CL converges at observed order 0.84–0.93, about first order, which is typical of a lattice with a load that falls to zero at the tips like a square root; CDi converges faster (order 2.3–3.2). The *ratio* K/K_free that ground effect is about moves less: 0.523, 0.531, 0.535 for ns = 20, 40, 80. The absolute values are not converged at ns = 40: Richardson extrapolation of CL (p = 0.93) gives CL_α = 4.625 per radian instead of 4.685 and e = 0.969 instead of 0.9945, so absolute CL, CDi and drags from the 4 × 40 lattice are 1–3 % off. The lesson's results are the ratios; where it quotes an absolute drag it gives the extrapolated value beside it.

### Compare

The lattice against the empirical estimate (2), h/b from 0.05 to 0.3 at α = 5°:

| h/b | CL/CL_free (fixed α) | CDi/CDi_free (fixed α) | K/K_free (lattice) | Wieselsberger (2) | single horseshoe (1), b_v = πb/4 |
|---|---|---|---|---|---|
| 0.05 | 1.391 | 0.686 | 0.354 | 0.342 | 0.061 |
| 0.075 | 1.237 | 0.695 | 0.454 | 0.439 | 0.127 |
| 0.1 | 1.166 | 0.722 | 0.531 | 0.515 | 0.206 |
| 0.15 | 1.101 | 0.779 | 0.643 | 0.629 | 0.369 |
| 0.2 | 1.070 | 0.825 | 0.720 | 0.709 | 0.509 |
| 0.3 | 1.041 | 0.885 | 0.817 | (outside range) | 0.700 |
| 0.5 | 1.019 | 0.943 | 0.909 | (outside range) | 0.866 |
| 1.0 | 1.005 | 0.982 | 0.972 | (outside range) | 0.963 |

The lattice sits 0.01–0.016 above the empirical fit across its stated range, the size of the lattice's own convergence spread (it rises with ns). The single horseshoe is far below both, for the reason in §2.2. At fixed α the induced drag still falls, despite the extra lift, because K falls faster than CL² rises.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction before running it.

1. **Ground effect against height.** The table above. Predict the h/b at which the induced drag at fixed lift is halved. *Finding:* about h/b = 0.09 for the lattice (0.531 at 0.1, 0.454 at 0.075).
2. **The trainer in the flare** (h = 1.0 m, h/b = 0.1, wing assumed 1 m above the runway, a low-wing layout). Predict the drag saving at 30 m/s. *Finding:* the weight needs CL = 1.112; the induced drag falls from 314 N to 167 N on the 4 × 40 lattice (323 N to 173 N extrapolated to a converged lattice); the angle of attack (from zero lift, flat lattice) for that CL falls from 13.60° to 11.66°; the effective span efficiency becomes 1.87.
3. **A WIG craft.** A 12 m × 4 m wing (AR 3), invented for this study, at h = 0.5, 1 and 2 m. Predict the lift-to-induced-drag ratio change at α = 5°. *Finding:* CL × 2.53, 1.68, 1.27; K × 0.267, 0.442, 0.648; L/D_i from 34.7 in free air to 51.5, 46.7 and 42.2. Low aspect ratio and small h/c make the lift change large, the reason WIG craft use short, broad wings.
4. **Downwash at the tail.** The trainer's tail (2.0 m², AR 4 assumed, span 2.83 m) 4.5 m behind the quarter chord. Predict dε/dα against the far-field 2CL_α/(πAR) = 0.358. *Finding:* 0.318 in the wing plane, 0.298 and 0.276 with the tail 0.5 and 1.0 m above it; along the centreline 0.484 at x = 2 m, 0.315 at 5 m, 0.284 at 10 m, 0.276 at 20 m, 0.273 at 100 m: the far value is below the elliptic-wing formula, because the rectangular wing's centre downwash is below its average.
5. **The tail in ground effect.** Wing at h = 1 m, tail at the same height. *Finding:* dε/dα falls from 0.318 to 0.124; at α = 5° the tail sees 0.97° less downwash; at α = 12°, 2.34° less, and (4) with V_H = 0.625 and a_t = 3.73 per radian (lattice) gives ΔCm = −0.095, nose-down.

Raw data format: CSV tables `h_over_b,CL_ratio,CDi_ratio,K_ratio,wieselsberger,single_horseshoe` (1); `quantity,free,ground` (2); `h_m,CL_ratio,K_ratio,L_over_Di_free,L_over_Di_ge` (3); `x_m,z_m,deps_dalpha` (4, 5); plus `ground_effect.png`.

## 5. Validate and improve

| check | expected | lattice | difference | cause |
|---|---|---|---|---|
| w on the ground plane | 0 | 1.33 × 10⁻¹⁵ relative | round-off | exact image construction |
| h → ∞ (h = 1000 m) | free air | CL, CDi, Γ(y) to 10⁻⁵ | → 0 | images far away |
| single horseshoe | (1) | to 10⁻¹² | round-off | same formula, discrete |
| Trefftz formula vs Biot–Savart at 10⁶ m | equal | within 10⁻⁶ relative (test tolerance) | the bound vortices still act at a finite distance | — |
| K/K_free at h/b = 0.1 | Wieselsberger 0.515 | 0.531 | 0.016 | empirical fit vs flat lattice; lattice not converged (0.523–0.535 for ns = 20–80) |
| dε/dα at the tail | 2CL_α/(πAR) = 0.358 | 0.318 | −11 % | tail close behind; rectangular load |

The engineering answer: **at 1 m above the runway the trainer's wing has 16.6 % more lift slope and 46.9 % less induced drag at the same lift (about 173 N instead of 323 N at 30 m/s, extrapolated to a converged lattice; 167 N and 314 N on the 4 × 40 lattice), and the tail's downwash gradient falls from 0.32 to 0.12, giving a nose-down ΔCm of about −0.10 at 12°, so the pilot must hold more up-elevator as the aircraft settles.** The induced-drag number is supported within a few hundredths by an independent empirical estimate; the moment change carries the larger uncertainty (tail slope, tail efficiency, the real flap and fuselage).

What to improve: flaps and the fuselage (the trainer lands with flaps, which load the inboard wing and change both lift and downwash), and the wing's real height (dihedral raises the tips). What not to improve: lattice density beyond ns = 40 for the ratios, which move by less than 0.01.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the image proof and the derivation of (1); `ground_effect.png`; the five experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) on the trainer's flare.

### Adversarial self-checks

- **What if h goes to infinity?** All ground-effect ratios → 1; the tests check CL, CDi and Γ(y) at h = 1000 m.
- **What if h goes to zero?** K → 0 in the model (the images cancel the wake), and CL grows without bound at fixed α: the flat-plate model has no thickness, no ground boundary layer and no channel flow, so it is wrong in this limit; real wings show a finite (and sometimes reversed) effect.
- **Which assumption fails first?** For the trainer, the flap and the fuselage; for a WIG craft, the flat stationary ground (waves); for both, a descent rate that makes the problem unsteady.
- **What measurement could falsify the model?** Flight-test drag polars flown at several heights over a runway (or a tunnel test with a moving belt): K/K_free against h/b should follow the lattice within a few hundredths; a lift *loss* near the ground would show the channel effect the model lacks.
- **Does the result survive a change of lattice?** K/K_free changes by 0.012 between ns = 20 and 80; nc from 2 to 8 changes it by less than 0.001.
- **Which decision changes if the uncertainty doubles?** A ±0.03 uncertainty in K does not change the flare technique; a ±0.2 uncertainty in ΔCm could, for an aircraft whose elevator authority at forward centre of gravity is marginal in the flare.

## Check yourself

1. Why does a mirror-image wing with reversed circulation model a flat ground exactly?
2. Explain, from the tangency condition and from the Trefftz plane, why ground effect raises lift at fixed α and lowers induced drag at fixed lift.
3. Roughly what fraction of its free-air induced drag does the trainer's wing keep at h/b = 0.1?
4. Name three independent checks of an image-system code.
5. What is the classical estimate of dε/dα, and why does it overestimate the downwash at the trainer's tail?
6. Which way does ground effect pitch a conventional aircraft, and why?
7. Is the steady ground-effect result an upper or a lower bound on the cushioning during a descent?
