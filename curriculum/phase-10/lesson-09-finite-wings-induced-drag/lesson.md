# P10.09  Finite Wings & Induced Drag

> A wing that makes lift must shed vortices from its tips; the vortices push air down, the lift tilts back, and that tilt is a drag no streamlining can remove. Only span reduces it.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the Biot–Savart velocity of a straight vortex segment and of a semi-infinite line, and check them against the infinite line Γ/(2πh). 2. Build the horseshoe vortex and verify Stokes' theorem and Kelvin's and Helmholtz's theorems numerically (circulation round each leg, none round the whole wake, a divergence-free field). 3. Assemble a lifting line from a stack of horseshoes, show that an elliptic load has uniform downwash Γ0/(2b) and CDi = CL²/(πAR), and compute induced drag two independent ways (near field and Trefftz plane) with an observed order of convergence. 4. Build the parabolic drag polar and find the best-L/D and minimum-power speeds of the reference trainer, stating the assumed CD0 and e. |
| **Prerequisites** | `P01.11` (vector fields, line integrals, circulation), `P01.12` (Stokes' theorem), `P10.03` (circulation, Kutta–Joukowski L′ = ρV∞Γ), `P10.04`, `P03.11` (observed order), `P08.08` by id only (potential flow; recap in section 2.1). Skill check: you can compute a cross product and evaluate ∫ dx/(x² + h²)^{3/2}. |
| **Engineering question** | At its 55 m/s cruise, how much of the reference trainer's drag is induced drag, what is the least induced drag its 10 m span allows, and at what speeds does it fly at best lift-to-drag ratio and at minimum power? Answer: an elliptic load would give CDi = CL²/(πAR) = 0.00418 (93 N); with an assumed whole-aircraft CD0 = 0.030 and Oswald efficiency e = 0.80, the cruise drag is 783 N of which 116 N (14.8 %) is induced; (L/D)max = 13.2 at 35.5 m/s and minimum power (17.4 kW) at 27.0 m/s, only just above the clean stall speed of `P10.08`. |
| **Assumptions and limits** | Incompressible, inviscid, steady; thin vortex filaments with a planar wake trailing straight downstream (no roll-up, no wake deflection); small induced angles (lifting-line view: the downwash only tilts the local lift). The drag polar is parabolic with constant CD0 and e, both assumed here (CD0 is built up in `P10.13`, e computed for the wing alone in `P10.10`–`P10.11`); no Reynolds-number or compressibility effects, no propeller. |
| **Analytical work** | The Biot–Savart law for a segment and its two limits (infinite and semi-infinite lines); Helmholtz's theorems and the horseshoe; the downwash of a single horseshoe and of the elliptic load; induced angle and induced drag; the Trefftz-plane theorem (far-field downwash twice the lifting-line value; D_i = ½ρ∫Γ w_T dy); CDi = CL²/(πeAR); the polar's optima: CL* = √(πeAR CD0), (L/D)max = ½√(πeAR/CD0), CL_mp = √3 CL*, V_mp = V_md/3^{1/4}. |
| **Computational deliverable** | `starter/main.py`: `segment_velocity`, `semi_infinite_velocity`, `horseshoe_velocity`, `loop_circulation`, `lifting_line_downwash`, `trefftz_downwash`, `lift_and_drag`, `span_efficiency`, `elliptic_gamma0`, `drag_polar`, `best_glide`, `speed_for_cl`, `drag_and_power` (`cosine_layout`, `elliptic_gamma` and `main` are given). Tests prove: the segment law equals the exact symmetric-segment result and tends to Γ/(2πh); the right-hand rule, linearity and reversal; the semi-infinite line is the limit of a long segment and gives half the infinite line at its end; circulation ±Γ round each leg, Γ round the bound leg, zero round the whole wake; zero divergence; a single horseshoe's downwash doubles downstream; the elliptic load's downwash is uniform, e = 1, and near-field drag equals Trefftz drag, converging at order 2; a two-term load gives Glauert's e = 1/(1 + δ); a uniform load's drag diverges; the polar's optima. |
| **Practical activity** | **computational experiment**: Biot–Savart verification, a convergence study of the horseshoe-stack elliptic wing against its closed form with near-field and far-field drag, and the trainer's drag and power curves with a sensitivity study of the assumed CD0 and e. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `induced_drag.png`, the derivations of the segment law, the elliptic downwash and the polar optima, the convergence table with its observed order, and a half-page note on the trainer's drag at cruise and its best speeds. |
| **Assessment** | `assessments/phase-10/grading-rubric.md`; competencies 1–4. Pass: the segment law and the elliptic results are derived; all public tests pass; induced drag is computed both ways with the observed order reported; the note states the assumed CD0 and e and how the best speeds change with them. |
| **References** | `@mit-ocw-16-unified-fluids` Fluids (spring) Lecture 5 "Introduction to 3D Wings, 2D and 3D Coefficients" and Lecture 6 "3D Vortex Filaments, Lifting Line Theory" (notes: General 3-D vortex; Wing vortex model; Downwash distribution) (https://ocw.mit.edu/courses/16-01-unified-engineering-i-ii-iii-iv-fall-2005-spring-2006/resources/f06_sp/) (required-free). `@mit-ocw-16-100-aerodynamics` Lecture 19 "Trefftz Plane Analysis of Induced Drag" and Lecture 32 "Single Horseshoe Vortex Wing Model" (https://ocw.mit.edu/courses/16-100-aerodynamics-fall-2005/resources/16100lectre19_cg/) (required-free). `@leishman-intro-aerospace-flight-vehicles` "Aerodynamics of Finite Wings" (https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/finite-wing-characteristics/): Origin of Trailing Vortices; Effects of Trailing Vortices; Drag Polar for a Finite Wing; Calculating the Lift & Drag on a Finite Wing (required-free). `@nasa-grc-beginners-guide-aeronautics` "Downwash Effects on Lift" and "Induced Drag Coefficient" (https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/induced-drag-coefficient/) (required-free, qualitative). `@ntrs-naca-tr-116-prandtl-modern-hydrodynamics` Part II "Applications": C. "The Finite Wing" (https://ntrs.nasa.gov/citations/19930091180) (supplementary-free): Prandtl's own account. Optional textbooks for deeper study (optional-paid, not required; library or used copy): Anderson, *Fundamentals of Aerodynamics*, for the horseshoe and the elliptic wing worked at length; Katz & Plotkin, *Low-Speed Aerodynamics*, for vortex filaments and the Trefftz plane in panel-method form. |
| **Extensions** | ★ Add a finite vortex core (a Rankine or Lamb–Oseen profile) to `segment_velocity` and show what it does to the self-induced velocity near a filament. ★ Put a second wing behind the first (a tail 4.5 m back) and compute the downwash it sees; compare with the far-field value 2α_i. ★ Let the wake roll up: march the trailing vortices downstream with their own induced velocity and watch the tip vortices form. |
| **Time** | unmeasured |

## 1. The engineering problem

A two-dimensional airfoil in inviscid flow makes lift with no drag at all (d'Alembert). A real wing has ends. At the tips the high pressure under the wing spills round to the low pressure above it, and the air leaves the trailing edge with a sideways twist that rolls up into the two **tip vortices** that anyone who has watched a crop-sprayer or an airliner on a humid day has seen. The vortices induce a downward velocity, the **downwash**, over the whole wing. Each section therefore meets the air at a slightly smaller angle than the wing's geometric angle, and its lift, perpendicular to the *local* flow, is tilted backwards. The backward component is **induced drag**. It is not a viscous loss and no surface finish removes it; it is the price of the kinetic energy the wing leaves behind in the wake every second.

For a slow aircraft induced drag dominates: on the trainer's approach it is two thirds of the drag. For a fast one it is small but never negligible: it sets the best-glide speed, the endurance speed, the climb rate and, through the wing's span, much of the aircraft's shape. This lesson builds the vortex model from the Biot–Savart law up, verifies it against the conservation theorems of vorticity, computes induced drag in two independent ways, and turns the result into the trainer's drag and power curves.

Success criteria, written first: the segment law reproduces the exact symmetric-segment result to round-off; circulation round each trailing leg equals Γ and round the whole wake is zero; the horseshoe stack reproduces the elliptic wing's uniform downwash Γ0/(2b) and e = 1, with near-field and far-field drag equal and converging to CL²/(πAR) at a measured order; the drag polar's best speeds are reported with the assumptions that set them.

## 2. Understand and derive

### 2.1 Vortex filaments (a recap of potential flow)

From `P08.08` (stated here so the lesson stands alone): in an incompressible flow with vorticity confined to thin filaments, the velocity is fixed by the vorticity alone, through the **Biot–Savart law**. A filament element dl of strength Γ induces at a point r away from it

dV = (Γ/4π) dl × r/|r|³.

Helmholtz's theorems constrain filaments: the strength Γ is constant along a filament, and a filament cannot end in the fluid; it must close on itself, run to infinity or end on a boundary. Kelvin's theorem adds that the circulation round a material loop does not change in time: a wing that starts from rest with zero circulation round a loop enclosing it and its wake must still have zero.

### 2.2 The straight segment

Integrate along a straight segment from A to B. With r0 = B − A, r1 = P − A, r2 = P − B:

**V = (Γ/4π) (r1 × r2)/|r1 × r2|² · r0·(r1/|r1| − r2/|r2|).**   (1)

Direction: r1 × r2 follows the right-hand rule about the segment. Size: for a point at distance h from the middle of a segment of length 2L, (1) gives V = (Γ/2πh) L/√(L² + h²), which tends to the **infinite line** Γ/(2πh) as L → ∞ with a relative error h²/(2L²) (5.0 × 10⁻⁵ at L/h = 100). Letting B run to infinity along a unit vector d gives the **semi-infinite line**

**V = (Γ/4π) (d × r1)/|d × r1|² (1 + d·r1/|r1|),**   (2)

which in the plane through its start is exactly half the infinite line: the factor (1 + cos) equals 1 there. Units: Γ in m²/s, distances in m, V in m/s.

### 2.3 The horseshoe and the conservation theorems

Prandtl's simplest wing is a **horseshoe**: a bound vortex along the span (vector along +y, positive Γ, so the Kutta–Joukowski force ρV × Γ points up) and two trailing legs that carry the same Γ to infinity downstream, as Helmholtz requires. The axes: x downstream, y to starboard, z up. The right leg runs from the bound vortex's end to +∞ with vector +x; the left leg comes in from +∞ (vector −x). Three checks follow from the theorems, and the code makes all three:

- **Stokes:** the circulation round a small loop encircling one leg is that leg's strength: +Γ round the right leg (counter-clockwise about +x), −Γ round the left, Γ round the bound vortex.
- **Kelvin:** a loop in a plane downstream that encloses both legs has zero circulation, as the starting wing had.
- **Continuity:** ∇·V = 0 everywhere off the filaments.

At the centre of the bound vortex each leg contributes Γ/(4π(b/2)) downwards, so a single horseshoe of span b has **w = Γ/(πb)** there; far downstream each leg is a full infinite line and the downwash doubles to 2Γ/(πb).

### 2.4 The lifting line, induced angle and induced drag

Replace the single horseshoe by many, one per strip of span, each with its own Γ_i; the trailing legs at a strip edge add up to the jump in Γ there, a **trailing vortex sheet** of strength −dΓ/dy. On the bound line the bound segments are collinear with the point and induce nothing; the downwash comes from the semi-infinite legs:

w(y) = (1/4π) ∫ (−dΓ/dη) dη/(y − η).   (3)

A section at y meets the air at a reduced angle α − α_i, α_i = w/V∞, and its lift ρV∞Γ, normal to the local flow, has a rearward component ρΓw per unit span. So

**D_i = ρ ∫ Γ w dy  (near field).**   (4)

Far downstream (the **Trefftz plane**) the legs are complete infinite lines, the downwash is twice (3), w_T = 2w for a flat wake, and the momentum and energy left in the wake give the same drag:

**D_i = ½ρ ∫ Γ w_T dy  (far field).**   (5)

The two expressions are computed by different formulas in the code (3-D Biot–Savart at the wing for (4); 2-D point vortices in a cross-flow plane for (5)), so their agreement is a check.

### 2.5 The elliptic wing (worked example)

Take Γ = Γ0√(1 − (2y/b)²). With y = −(b/2) cos θ, Γ = Γ0 sin θ, and (3) becomes Glauert's integral, whose value is constant: **w = Γ0/(2b)** across the whole span. Lift L = ρV∞Γ0 πb/4, so Γ0 = 2CL V∞ S/(πb), and

α_i = w/V∞ = CL/(πAR),  **CDi = CL α_i = CL²/(πAR).**   (6)

For any other load the drag is larger: CDi = CL²/(πeAR) with the **span efficiency** e ≤ 1 (equality only for the elliptic load; `P10.10` derives e = 1/(1 + δ)). For the trainer in cruise (W = 7357.5 N, V = 55 m/s, S = 12 m², b = 10 m, AR = 8.333): CL = 0.3309, Γ0 = 13.90 m²/s, **w = 0.695 m/s, α_i = 0.724°, CDi = 0.00418**, an induced drag of 93 N even for a perfect elliptic load. Notice D_i = L²/(½ρV²πeb²): at fixed weight and speed only the **span**, not the area, controls it.

### 2.6 The drag polar and the best speeds

With a parasite drag coefficient CD0 (skin friction and pressure drag of the whole aircraft at zero lift, `P10.13`) and an Oswald efficiency e for the whole aircraft:

**CD = CD0 + CL²/(πeAR).**   (7)

L/D = CL/CD is largest where d(CD/CL)/dCL = 0, i.e. where induced drag equals parasite drag: **CL* = √(πeAR CD0), (L/D)max = ½√(πeAR/CD0)**. Power P = DV ∝ CD/CL^{3/2} is least where CDi = 3CD0: **CL_mp = √3 CL***, so V_mp = V_md/3^{1/4}. Assumed for the trainer: CD0 = 0.030 (a fixed-gear light aircraft) and e = 0.80. They give (L/D)max = 13.21 at CL = 0.793, V_md = 35.5 m/s, and V_mp = 27.0 m/s. The best-L/D speed is the best-glide and best-range speed (for a propeller aircraft at constant efficiency); the minimum-power speed is the best-endurance speed.

### 2.7 Where the model stops

- The wake is assumed flat and straight; real wakes roll up within a few spans and descend. The drag computed in the Trefftz plane does not depend on that (it depends on the circulation shed), but the downwash at a tail does.
- A single horseshoe has infinite Trefftz drag: two point vortices of finite strength carry infinite kinetic energy per unit length. Real tip vortices have viscous cores.
- Lifting-line downwash assumes a high aspect ratio: each section sees a slowly varying downwash. Low-aspect-ratio and swept wings need the lattice of `P10.11`.
- CD0 and e in (7) are assumptions here; the polar's numbers are only as good as they are (part 4 tests that).

## 3. Build the computational model

Open `starter/main.py`. Given: `cosine_layout` (strip edges at y = −(b/2) cos φ_k with control points at the φ-midpoints), `elliptic_gamma`, `main`.

### Step 1: Biot–Savart

`segment_velocity(P, A, B, gamma)` is (1), vectorised to return the contribution of every segment at every point, shape (points, segments, 3); a point on a filament's line gets zero. `semi_infinite_velocity` is (2).

```python
c = np.cross(r1, r2)
k = np.sum(r0 * (r1 / n1[..., None] - r2 / n2[..., None]), axis=-1) / np.sum(c * c, axis=-1)
v = gamma / (4 * np.pi) * k[..., None] * c
```

### Step 2: horseshoes and circulation

`horseshoe_velocity` adds a bound segment and two semi-infinite legs per horseshoe. `loop_circulation` integrates V·dl round a circle by the midpoint rule, which converges spectrally for a smooth periodic integrand.

### Step 3: the lifting line from horseshoes

`lifting_line_downwash(edges, gamma, y_control)` evaluates (3) by 3-D Biot–Savart at the control points; `trefftz_downwash` sums 2-D vortices of strength Γ_left − Γ_right at the strip edges; `lift_and_drag` returns L, (4) and (5). Place the control points with care: with edges on the cosine layout and control points at the φ-midpoints (the "semicircle" rule), a sampled elliptic load gives *exactly* uniform ratios and e = 1; at the y-midpoints the control point of the tip strip sits almost on the tip vortex.

### Step 4: polar and speeds

`drag_polar`, `best_glide`, `speed_for_cl` and `drag_and_power` are section 2.6.

Running `main()` prints (abridged):

```
finite segment +-L at distance h = 0.5 m against the infinite line Gamma/(2 pi h):
  L/h =   100: v / (Gamma/(2 pi h)) = 0.99995000  (1 - ratio = 5.000e-05; exact L/sqrt(L^2 + h^2) gives 5.000e-05)
circulation round the right trailing leg 1.000000, left leg -1.000000, bound leg 1.000000, whole wake 2.78e-17  (Gamma = 1)
reference trainer cruise: CL = 0.3309, AR = 8.333, elliptic Gamma_0 = 13.904 m^2/s, w = Gamma_0/(2b) = 0.6952 m/s, alpha_i = CL/(pi AR) = 0.724 deg
  n strips   max|w/w_exact - 1|   CDi near     CDi Trefftz   CL^2/(pi AR)   e (Trefftz)
     10      4.107e-03         0.0041485    0.0041485     0.0041828      1.000000000000
     20      1.028e-03         0.0041742    0.0041742     0.0041828      1.000000000000
     40      2.570e-04         0.0041806    0.0041806     0.0041828      1.000000000000
     80      6.425e-05         0.0041823    0.0041823     0.0041828      1.000000000000
    160      1.606e-05         0.0041827    0.0041827     0.0041828      1.000000000000
  observed order of the CDi error: [1.996 1.999 2.    2.   ]
  same 40 strips with control points at the y-midpoints instead: e = 1.0314, max|w/w_exact - 1| = 8.178e+00
single horseshoe, same lift: Gamma = 10.920, centre downwash 0.3476 m/s (Gamma/(pi b) = 0.3476) vs elliptic 0.6952
```

and writes `induced_drag.png`: the stepped and the exact elliptic load, the downwash at the wing and in the Trefftz plane, and the trainer's drag curves.

### Verify

`python -m pytest curriculum/phase-10/lesson-09-finite-wings-induced-drag` runs 13 tests, listed in the metadata. The two-term test is worth noticing: for Γ ∝ sin φ + a₃ sin 3φ the horseshoe stack, which knows nothing of Fourier series, returns e = 1/(1 + 3a₃²) to 2 × 10⁻⁴, the result `P10.10` derives.

### Compare

Two independent routes agree. The near-field drag (3-D Biot–Savart at the bound line) and the far-field drag (2-D vortices in the Trefftz plane) are equal to 12 digits for every layout, as the planar-wake theorem w_T = 2w says they must. Both converge to the closed form (6) at order 2.00, and the ratio CL²/(πAR CDi) is 1 to 12 digits at every resolution, because the semicircle layout integrates the elliptic load's lift and drag with errors that cancel in the ratio. The cautionary row is the y-midpoint layout: e = 1.031, *above* the theoretical maximum. A span efficiency above 1 for a planar wing is always a numerical error, and a useful alarm.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction before running it.

1. **Segment to line.** v/(Γ/2πh) for L/h = 10, 100, 1000. Predict the error. *Finding:* 4.963 × 10⁻³, 5.000 × 10⁻⁵, 5.000 × 10⁻⁷, matching 1 − L/√(L² + h²) ≈ h²/(2L²).
2. **Convergence of the elliptic wing.** n = 10–160 strips. *Finding:* the table above; observed order 1.996, 1.999, 2.000, 2.000; the downwash error falls from 4.1 × 10⁻³ to 1.6 × 10⁻⁵.
3. **One horseshoe vs elliptic.** Same lift, same span. *Finding:* the single horseshoe's centre downwash is 0.348 m/s, half the elliptic wing's 0.695 m/s, but its Trefftz drag is infinite; a uniform load on 10, 40, 160 strips shows the drag growing without bound (the test checks e falling). Concentrated tip vortices are expensive.
4. **The trainer's drag and power** (assumed CD0 = 0.030, e = 0.80):

| V (m/s) | CL | D (N) | parasite (N) | induced (N) | induced share | P (kW) | L/D |
|---|---|---|---|---|---|---|---|
| 55.0 (cruise) | 0.3309 | 783.3 | 667.0 | 116.2 | 14.8 % | 43.08 | 9.39 |
| 35.5 (best L/D) | 0.7927 | 556.9 | 278.5 | 278.5 | 50.0 % | 19.79 | 13.21 |
| 30.0 (approach) | 1.1122 | 589.2 | 198.4 | 390.7 | 66.3 % | 17.68 | 12.49 |
| 27.0 (min power) | 1.3729 | 643.1 | 160.8 | 482.3 | 75.0 % | 17.36 | 11.44 |

Predict the induced share at cruise first.

5. **Sensitivity to the assumptions.** *Finding:* CD0 0.025 → (L/D)max 14.47, V_md 37.2 m/s; CD0 0.035 → 12.23, 34.2 m/s; e 0.70 → 12.36, 36.7 m/s; e 0.90 → 14.01, 34.5 m/s. A ±17 % change in CD0 moves (L/D)max by 7–10 % (it scales as CD0^(−1/2)) and V_md by about 4 % (CD0^(−1/4)).

Raw data format: CSV tables `L_over_h, ratio` (1); `n, w_error, CDi_near, CDi_trefftz, e` (2); `V, CL, D, D0, Di, P, LD` (4); `CD0, e, LDmax, Vmd` (5); plus `induced_drag.png`.

## 5. Validate and improve

| quantity | theory (§2) | code | difference | cause |
|---|---|---|---|---|
| segment, L/h = 100 | 0.99995000 | 0.99995000 | < 10⁻¹² | — |
| circulation round a leg | Γ = 1 | 1.000000 | < 10⁻⁹ | — |
| wake circulation | 0 | 2.8 × 10⁻¹⁷ | round-off | — |
| elliptic CDi, 160 strips | 0.0041828 | 0.0041827 | −1.1 × 10⁻⁷ | discretisation, order 2 |
| near vs far-field drag | equal | equal to 12 digits | — | planar wake |
| e, elliptic | 1 | 1.000000000000 | < 10⁻¹² | semicircle layout |

The vortex model is verified. Its validity for the trainer's *wing* is the subject of `P10.10`–`P10.11`, which compute e for the rectangular planform; its validity for the *aircraft* polar rests on the assumed CD0 and e, which `P10.13` replaces with a drag build-up. The engineering answer: **in cruise the trainer's drag is about 783 N, of which 116 N (15 %) is induced; even a perfect elliptic load would leave 93 N, set by the 10 m span. Best L/D (13.2) is at 35.5 m/s and minimum power (17.4 kW) at 27.0 m/s, which is only 2 % above the clean stall speed of `P10.08` (26.4 m/s): the trainer cannot actually fly at its minimum-power speed clean, which is common for light aircraft and one more reason for flaps.**

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (1), (2), (6) and the polar optima; `induced_drag.png`; the five experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) giving the cruise drag split, the best speeds, and the assumptions they depend on.

### Adversarial self-checks

- **What if the span doubles at constant area?** At fixed weight and speed induced drag falls to a quarter (D_i ∝ 1/b²); at cruise that saves about 87 N of 783 N, and the wing's root bending moment grows (`P10.12`, `P05.12`).
- **What if the speed goes to zero?** CL → ∞ and induced drag ∝ 1/V² grows without bound; the polar is meaningless below the stall speed, where CL cannot exceed CL_max.
- **Which assumption fails first?** The assumed e = 0.80 for the whole aircraft (the wing alone computes to about 0.93–0.97 in `P10.10`–`P10.11`; the fuselage and the trim load lower it).
- **What measurement could falsify the model?** A glide test: measured sink rate against speed gives CD(CL), and the slope of CD against CL² gives πeAR directly.
- **Does the result survive a change of method or resolution?** The elliptic drag agrees by two methods and converges at order 2; the trainer's polar depends on assumptions, not on resolution.
- **Which decision changes if the uncertainty doubles?** If CD0 is uncertain by ±0.010 instead of ±0.005, the best-glide speed spans about 33–39 m/s; a glide-speed placard would need a flight test.

## Check yourself

1. State Helmholtz's two vortex theorems used in building a horseshoe.
2. What is the velocity at distance h from an infinite line vortex, and from the end-plane of a semi-infinite one?
3. Why does an elliptic load have the same downwash at every station?
4. Write CDi for a wing of span efficiency e and explain why at fixed weight and speed only the span matters.
5. Why is the Trefftz-plane downwash twice the downwash at the lifting line?
6. At best L/D, how do induced and parasite drag compare? At minimum power?
7. A code reports e = 1.03 for a planar wing. What do you conclude?
8. Why does a single horseshoe have infinite induced drag?
