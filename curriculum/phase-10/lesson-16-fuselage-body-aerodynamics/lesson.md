# P10.16  Fuselage & Body Aerodynamics

> A closed body in potential flow feels no lift, only a nose-up couple: the fuselage pushes the neutral point forward, and the tail has to pay for it.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the slender-body loading dN/dx = 2qα dS/dx from the apparent mass of a cross-section, integrate it to C_Nα = 2 on the base area and to zero for a closed body, and derive Munk's moment dM/dα = ρV²·Vol. 2. Implement and test the slender-body integrals, the exact potential moment of a prolate spheroid from Lamb's added-mass coefficients, the Allen–Perkins viscous cross-flow correction, and body drag from skin friction, wetted area and an empirical form factor. 3. Show numerically that slender-body theory is the limit of the exact spheroid solution as the fineness ratio grows, and quantify its error at a real fuselage's fineness. 4. Compute the reference trainer's fuselage Cm_α and the forward shift of its neutral point, bracketing it between a potential-flow and a separated-tail-cone estimate. 5. Judge which fineness ratio minimises drag for a given volume and how much the answer depends on the empirical correlation. |
| **Prerequisites** | `P10.01` (moment reference and coefficients), `P10.09` and `P10.10` (wing lift slope, downwash), `P10.13` (skin friction and form factors), `P03.04` (numerical integration), `P01.11` (circulation, momentum), `P00.14`. From Phase 08 (recapped in §2.1): `P08.08` (potential flow, added mass). Skill check: you can integrate by parts and compute the volume of a solid of revolution. |
| **Engineering question** | How much does the reference trainer's 7 m fuselage (maximum diameter 1.1 m) destabilise it in pitch, and how far does it move the neutral point forward of the wing-and-tail value? The answers: the slender-body (potential) fuselage contributes Cm_α = +0.56 /rad about the wing quarter chord (0.49 /rad with the exact-spheroid correction k₂ − k₁ = 0.88), moving the neutral point from 0.505 c to 0.412 c, 11 cm forward; a model that drops the separated tail cone's potential loading gives +0.22 /rad and a shift of 0.048 c. The fuselage's drag is ΔC_D = 0.0049. |
| **Assumptions and limits** | Inviscid, incompressible flow; slender body (radius ≪ length, small slopes), small α; the body is a solid of revolution (the fuselage's cross-section area is represented by a circle of the same area); no wing–body interference (the wing's upwash ahead and downwash behind it change the local α on the fuselage); the viscous cross-flow coefficients are empirical. Slender-body theory fails near a blunt nose and wherever the boundary layer separates (a steep tail cone, high α). |
| **Analytical work** | Apparent mass and the slender-body loading; integration by parts for N and M; Munk's moment and its independence of the reference point; the paraboloid nose (uniform loading, centre of pressure at L/2); the exact spheroid moment q·Vol(k₂ − k₁) sin 2α and its slender limit; the neutral-point formula with a fuselage term; body drag scaling with fineness at fixed volume. |
| **Computational deliverable** | `starter/main.py`: `section_area`, `body_volume`, `wetted_area`, `slender_body_loading`, `normal_force_and_moment`, `slender_body_closed_form`, `lamb_coefficients`, `spheroid_munk_moment`, `allen_perkins_cn`, `spheroid_drag_area`, `best_fineness`, `fuselage_pitch_derivatives`, `trainer_neutral_point` (`ellipsoid_radius`, `paraboloid_radius`, `trainer_fuselage_radius`, `turbulent_cf`, `helmbold_lift_slope`, `body_form_factor` and `main` are given). Tests prove: volume and wetted area against the sphere and spheroid closed forms; zero normal force and a reference-independent couple 2qα·Vol on a closed body; linearity in α and q; the paraboloid's uniform loading, C_Nα = 2 and centre of pressure at L/2; numerical loading against the closed form for an open body; Lamb's coefficients for a sphere and their limits; the exact moment approaching slender-body theory from below as fineness grows; Allen–Perkins reducing to potential theory at small α and growing as sin²α; an interior drag optimum at fixed volume; the fuselage destabilising and moving the neutral point forward. |
| **Practical activity** | **computational experiment**: slender-body integrals against closed forms; slender-body theory against the exact spheroid as fineness varies; the trainer's fuselage Cm_α under two loading models and the neutral-point shift; the fineness ratio of least drag; the Allen–Perkins correction at larger α. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `fuselage_slender_body.png`, the derivations of §2, the experiment tables of §4 and a half-page note giving the trainer's neutral-point shift with its bracket and the assumption that sets the bracket. |
| **Assessment** | `assessments/phase-10/grading-rubric.md`; competencies 1–4. Pass: the loading and Munk's moment are derived; all public tests pass; the spheroid convergence is reported with the error at the fuselage's fineness; the neutral-point shift is given as a bracket with its cause. |
| **References** | `@ntrs-naca-tr-184-munk-airship-hulls` (Munk, *The Aerodynamic Forces on Airship Hulls*, https://ntrs.nasa.gov/citations/19930091249): I. General Properties of Aerodynamic Flows; II. The Aerodynamic Forces on Airship Hulls; III. Some Practical Conclusions (required-free: apparent mass, the unstable moment of a hull). `@ntrs-naca-tr-1048-slender-inclined-bodies` (Allen & Perkins, https://ntrs.nasa.gov/citations/19930090962): Introduction; Potential Flow Theory; Comparisons with Experiment and Discussion of the Effects of Viscosity; Estimation of Cross-Flow Drag Coefficient (required-free: the viscous cross-flow addition). `@ntrs-naca-tr-116-prandtl-modern-hydrodynamics` Part II, A. "Distribution of Pressure on Airship Bodies" (supplementary-free). `@leishman-intro-aerospace-flight-vehicles` "Airships, Blimps, & Aerostats" (https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/lth/): Aerodynamic Effects (supplementary-free, link only). `@nasa-glenn-beginners-guide-aeronautics` "Fuselage" (https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/fuselage/) (supplementary-free). `@drela-avl` AVL User Primer: General Description; the `BODY` keyword in "Geometry Input File -- xxx.avl" (supplementary-free: how a vortex-lattice code represents a fuselage by slender-body theory). Optional textbooks for deeper study: M. Drela, *Flight Vehicle Aerodynamics* (slender-body theory in modern notation, including wing–body interaction); S. F. Hoerner, *Fluid-Dynamic Drag* and D. P. Raymer, *Aircraft Design: A Conceptual Approach* (the origins of the body form-factor correlations used here). |
| **Extensions** | ★ Add the wing's upwash ahead of it and downwash behind it to the fuselage's local α (Multhopp's idea) using the lifting-line downwash of `P10.10`, and recompute Cm_α. ★ Replace the spheroid correction by a numerical axisymmetric source-panel solution of the actual fuselage profile and compare its moment with slender-body theory. ★ Compute the body's apparent mass in heave and its effect on the short-period mode of `P09`-style dynamics. |
| **Time** | unmeasured |

## 1. The engineering problem

A wing alone is stable in pitch if its centre of gravity is ahead of its aerodynamic centre; a tail moves the **neutral point**, the centre-of-gravity position at which the aircraft is neutrally stable, far aft, so that the aircraft can carry a useful range of centre-of-gravity positions. The fuselage works against the tail. It produces little lift, but in a flow at incidence it develops a nose-up moment that grows with α: a destabilising Cm_α. The designer has to know how big it is, because the tail must be sized to overcome it and the centre-of-gravity limits are set from the resulting neutral point. A 5 % error in the neutral point is a 5 % error in the static margin, the quantity on which handling, trim drag and certification limits depend.

The tool is **slender-body theory**: for a long thin body, each cross-section sees a two-dimensional cross-flow V sin α, and the change of that section's apparent mass along the body gives the local force. It is due to Munk, who developed it for airship hulls in 1924; it is still how vortex-lattice codes such as AVL represent fuselages. This lesson derives it, checks it against the exact potential flow round a prolate spheroid, adds the empirical viscous correction of Allen and Perkins, and applies it to the reference trainer: fuselage length 7.0 m, maximum diameter 1.1 m, wing area 12.0 m², chord 1.2 m, a 2.0 m² tail 4.5 m behind the wing quarter chord, the wing's leading edge assumed at 2.0 m from the nose.

Success criteria, written first: the slender-body integrals reproduce the closed forms to 10⁻⁵ (zero lift and 2qα·Vol on a closed body; C_Nα = 2 and centre of pressure at L/2 on a paraboloid); the exact spheroid moment approaches the slender-body moment as the fineness grows, with the error falling below 0.5 % by L/D = 64; the neutral-point shift is reported as a bracket with its cause.

## 2. Understand and derive

### 2.1 Recap: apparent mass (`P08.08`)

A body accelerating through a fluid must also accelerate some of the fluid; the extra inertia is the body's **apparent (added) mass**. For a circular cylinder moving sideways in 2-D potential flow, the apparent mass per unit length is exactly ρπr², the mass of fluid the cylinder displaces. For a body of revolution moving along its axis or sideways, the added mass is k₁ρ·Vol or k₂ρ·Vol, with k₁, k₂ depending on the shape (§2.4).

### 2.2 Slender-body loading

Put the body at a small angle α in a stream V; x runs from the nose aft, S(x) = πr(x)² is the cross-section area. Split the flow into an axial part V cos α ≈ V, which produces no normal force in potential flow, and a cross-flow w = V sin α ≈ Vα. In the cross-flow plane each section is a circle moving sideways at w, carrying vertical momentum ρS(x)w per unit length (its apparent mass times w).

Now ride with a slab of air as it passes the body at speed V. At station x it carries ρS(x)w of downward momentum; a moment later, at x + dx, it carries ρS(x + dx)w. The body must have supplied the change: the force per unit length on the body is

**dN/dx = d/dx(ρS w)·V = ρV²α dS/dx = 2qα dS/dx**,   (1)

positive (up) where the body grows and negative where it shrinks. Units: Pa × m²/m = N/m.

**Normal force.** Integrating over the length L,

**N = 2qα [S(L) − S(0)] = 2qα S_base,  C_Nα = N/(qαS_base) = 2.**   (2)

A body ending in a full section (a nose, a missile without a boat tail) carries lift 2qαS_base; a **closed body** (S(L) = 0) carries **none**: the aft body's negative loading cancels the nose's positive loading exactly.

**Moment.** About a point x_ref, nose-up positive, an upward force at x gives N(x_ref − x). Integrating (1) by parts,

M = 2qα ∫₀ᴸ S′(x)(x_ref − x) dx = 2qα {[S(x)(x_ref − x)]₀ᴸ + ∫₀ᴸ S dx} = **2qα (Vol − S_base(L − x_ref)).**   (3)

For a closed body the base term vanishes and

**M = 2qα·Vol,  dM/dα = ρV²·Vol (Munk's moment),**   (4)

a pure couple, the same about every point, nose-up for positive α: **destabilising**. The physical picture: the nose region pushes up, the aft region pushes down by the same amount, and the two make a couple.

**A worked example: the paraboloid nose.** r = R√(x/L) gives S = πR²x/L, so dS/dx = πR²/L is constant: the loading is uniform, N = 2qαπR² (C_Nα = 2 on the base) and the centre of pressure is at L/2. For L = 2.0 m, R = 0.55 m, q = 1852.8 Pa (55 m/s) and α = 4°, N = 245.85 N and the code's centre of pressure is 1.0000 m.

### 2.3 Body drag and fineness

The body's profile drag is estimated as in `P10.13`: D/q = C_f·FF·S_wet with the turbulent flat-plate C_f at Re = VL/ν, the wetted area S_wet = ∫2πr√(1 + r′²) dx, and an **empirical form factor** FF = 1 + 60/f³ + f/400 with f = L/D the **fineness ratio** (the 60/f³ term is the pressure drag of bluntness, the f/400 term the extra friction of a longer body). At fixed volume, a short fat body has little wetted area but a large form factor; a long thin one the reverse; there is a fineness of least drag (Experiment 4).

### 2.4 The exact check: a prolate spheroid

For a prolate spheroid (an ellipse of revolution) potential flow is known exactly. With e = √(1 − (D/L)²),

a₀ = (2(1 − e²)/e³)(artanh e − e),  b₀ = 1/e² − ((1 − e²)/e³) artanh e,  k₁ = a₀/(2 − a₀),  k₂ = b₀/(2 − b₀),   (5)

are Lamb's classical added-mass coefficients along and across the axis (a sphere has k₁ = k₂ = ½). The moment on the spheroid is

**M = q·Vol (k₂ − k₁) sin 2α,**   (6)

which for small α is 2qα·Vol(k₂ − k₁). As L/D → ∞, k₁ → 0 and k₂ → 1, and (6) becomes Munk's (4). Slender-body theory therefore over-predicts the potential moment by the factor 1/(k₂ − k₁): 13 % at the trainer's fineness 6.36 (k₂ − k₁ = 0.882).

### 2.5 Viscosity: the Allen–Perkins correction

Real bodies at incidence shed vortices from their lee side, and the cross-flow behaves like the flow past a cylinder, with a drag coefficient c_dc ≈ 1.2. Allen and Perkins (NACA TR-1048) added that cross-flow drag to the potential normal force:

C_N = (S_base/S_ref) sin 2α cos(α/2) + η c_dc (A_plan/S_ref) sin²α,   (7)

with A_plan the planform (side-view) area and η ≈ 0.6–0.7 a finite-length factor. Both coefficients are empirical. The added term grows as α², so it matters above a few degrees (Experiment 5).

A related effect matters even at small α: the thick, possibly separated boundary layer on a tail cone cannot carry the negative potential-flow loading of (1). A common engineering approximation keeps the loading only ahead of the maximum section and drops it behind. That is the second, "cut" model of the trainer's fuselage.

### 2.6 The trainer's neutral point

Lift slopes on the wing area: the wing a_w by Helmbold's formula at AR = 8.33 (4.954 /rad); the tail a_t at an assumed AR of 4 (3.883 /rad), efficiency η_t = 0.9, the downwash gradient dε/dα = 2a_w/(πAR) = 0.378, and the tail volume V_H = S_t l_t/(S c) = 0.625. With positions h in chords aft of the wing leading edge, the neutral point is where dCm/dα = 0:

**h_n = (Σ a_i h_i − Cm_α,fus) / Σ a_i,**   (8)

where a_i are the wing's, the tail's (η_t(S_t/S)a_t(1 − dε/dα)) and the fuselage's normal-force slopes and Cm_α,fus is the fuselage moment derivative about the wing quarter chord. A positive Cm_α,fus moves h_n forward by Cm_α,fus/Σa_i.

### 2.7 Where the model stops

- **Wing–body interference.** The wing's upwash raises the fuselage's local α ahead of it and its downwash lowers it behind: the forward fuselage is more destabilising than the free-stream α suggests (Extension 1).
- **Non-circular sections**, canopies and a propeller in front are represented only through the equivalent area.
- **Separation** on the tail cone and at high α: the two loading models bracket the first; (7) estimates the second.

## 3. Build the computational model

Open `starter/main.py`. The geometry helpers (`ellipsoid_radius`, `paraboloid_radius`, `trainer_fuselage_radius`), `turbulent_cf`, `helmbold_lift_slope`, `body_form_factor` and `main` are given. The trainer's profile is an assumption: an elliptic nose from 0 to 1.8 m, a cabin of constant radius 0.55 m to 3.2 m, and a tail cone r = 0.55 − 0.47((x − 3.2)/3.8)^{1.3} ending in a 0.08 m base.

### Step 1: geometry integrals

`section_area`, `body_volume` (trapezoid rule on S) and `wetted_area` (the chain of frustums through the profile points, exact for the polyline).

### Step 2: slender-body loading and its integrals

```python
def slender_body_loading(x, r, alpha, q):
    return 2 * q * alpha * np.gradient(np.pi * r**2, x)      # (1)
```

`normal_force_and_moment` integrates dN/dx and dN/dx·(x_ref − x); `slender_body_closed_form` is (2)–(3). Use a cosine-clustered grid for a closed body, so that the steep ends are resolved.

### Step 3: the exact spheroid and Allen–Perkins

`lamb_coefficients(f)` is (5) with the sphere as a special case; `spheroid_munk_moment` is (6); `allen_perkins_cn` is (7).

### Step 4: drag, fineness and the trainer

`spheroid_drag_area(f, volume, V)` sizes a spheroid of the given volume and fineness (Vol = (π/6)LD²) and returns C_f, FF, S_wet and D/q; `best_fineness` minimises it with `minimize_scalar`. `fuselage_pitch_derivatives(x, r, x_ref, viscous_cut)` returns C_Nα and Cm_α per radian on S and c; `trainer_neutral_point` is (8).

Running `main()` prints (abridged):

```
spheroid L = 7.0, D = 1.1: x_ref = 0.0: N = 7.958e-13 N, M = 1147.31 N m; closed form 2 q alpha Vol = 1147.31 N m
spheroid L = 7.0, D = 1.1: x_ref = 3.5: N = 7.958e-13 N, M = 1147.31 N m; closed form 2 q alpha Vol = 1147.31 N m
paraboloid nose 2.0 m: N = 245.85 N (closed 245.85), centre of pressure 1.0000 m (closed 1.0000 m = L/2); C_N_alpha on base area = 2.0000
prolate spheroid: exact potential moment / slender-body moment = k2 - k1:
  L/D =   2.00: k1 = 0.2100, k2 = 0.7042, k2 - k1 = 0.4942
  L/D =   6.36: k1 = 0.0414, k2 = 0.9235, k2 - k1 = 0.8821
  L/D =  20.00: k1 = 0.0068, k2 = 0.9866, k2 - k1 = 0.9798
trainer fuselage: volume 4.131 m^3, wetted area 18.21 m^2, fineness 6.36, base area 0.0201 m^2
slender body, about the wing quarter chord (x = 2.3 m): C_N_alpha = 0.0034, C_m_alpha = 0.5607 /rad; 2 Vol/(S c) = 0.5738
  with the spheroid correction k2 - k1 = 0.8822: C_m_alpha = 0.4946 /rad
  loading cut behind the maximum section: C_N_alpha = 0.1584, C_m_alpha = 0.2244 /rad
neutral point (chords aft of wing LE): wing + tail 0.5054; + fuselage (spheroid-corrected slender body) 0.4123 (shift -0.0931 c = -11.2 cm); + fuselage (cut loading) 0.4570 (shift -0.0484 c)
fuselage drag: Re = 2.64e+07, Cf = 0.00258, FF = 1.2487, D/q = 0.0588 m^2, Delta C_D = 0.00490
spheroid of the same volume: least drag at L/D = 8.75
```

and writes `fuselage_slender_body.png`: the fuselage profile, its slender-body loading at 4° and 55 m/s, and k₂ − k₁ against fineness.

### Verify

Run `python -m pytest curriculum/phase-10/lesson-16-fuselage-body-aerodynamics`. The tests check the sphere's 4πR³/3 and 4πR², the spheroid's surface area; zero N and the same couple about five reference points; linearity; the paraboloid's uniform loading, C_Nα = 2 and L/2; the trainer's open body against (2)–(3); Lamb's coefficients at the sphere, their monotone limits; the exact moment below the slender value with the gap shrinking to under 0.5 % at L/D = 64; Allen–Perkins at small α and its sin²α growth; the drag optimum; the fuselage's positive Cm_α near 2Vol/(Sc) and the forward neutral-point shift; the bookkeeping of (8) (a pure moment shifts h_n by −ΔCm_α/ΣCLα; no tail gives h_n = 0.25).

### Compare

The trainer's fuselage is almost closed (base 0.020 m²), so its Cm_α of 0.5607 /rad is within 2.3 % of 2Vol/(Sc) = 0.5738: the base term of (3) is small. On the 7 m spheroid the numerical couple equals 2qα·Vol = 1147.31 N·m about the nose, the middle and the tail. The exact spheroid at the same fineness gives 0.882 of that: slender-body theory is a 13 % over-estimate here, and that correction is applied before the neutral-point estimate.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction before running.

1. **Slender-body theory against the exact spheroid.** k₂ − k₁ at L/D = 1–50. Predict the fineness at which slender-body theory is within 5 %. *Finding:* k₂ − k₁ = 0, 0.494, 0.778, 0.882, 0.940, 0.980, 0.996 at L/D = 1, 2, 4, 6.36, 10, 20, 50: within 5 % only beyond L/D = 11.3. Typical fuselages (L/D ≈ 6–8) are in the range where the correction is 10–15 %.

2. **The trainer's fuselage Cm_α.** Three models about the wing quarter chord. Predict the spread. *Finding:* slender body 0.561 /rad, spheroid-corrected 0.495 /rad, cut behind the maximum section 0.224 /rad (with C_Nα = 0.158, the nose's lift that is no longer cancelled). Neutral point: 0.505 c (wing and tail), 0.412 c (corrected slender body, 11.2 cm forward), 0.457 c (cut, 5.8 cm forward). The truth lies between: the tail cone's boundary layer carries part of the potential loading, not none of it. A designer would take the larger shift for sizing the tail.

3. **Drag of the fuselage.** At 55 m/s: Re = 2.64 × 10⁷, C_f = 0.00258, FF = 1.249, S_wet = 18.21 m², D/q = 0.0588 m², ΔC_D = 0.0049, about 15 % of the trainer's C_D0 in `P10.13`.

4. **Fineness of least drag at fixed volume.** A spheroid of the fuselage's volume (4.131 m³), L/D = 2–10. Predict the optimum. *Finding:* D/q = 0.333, 0.133, 0.0566, 0.0593, 0.0569 m² at L/D = 2, 3, 8.75, 6.36, 10: the minimum is at **L/D = 8.75** and is flat: the trainer's 6.36 costs only 5 % more drag. Problem 5 shows that another published correlation puts the optimum near 6, so the location of a flat minimum is the least reliable output.

5. **Viscous cross-flow.** C_N on the maximum section for a body carrying its potential loading on that section, with A_plan = 5.72 m², at α = 2–30°. *Finding:* potential 0.070, 0.174, 0.341, 0.633, 0.837; with Allen–Perkins 0.076, 0.209, 0.482, 1.182, 2.011 at 2, 5, 10, 20, 30°. The viscous term adds 8 % at 2° and exceeds the potential term above 22.7°.

Raw data format: CSV tables `fineness, k1, k2, k2_minus_k1` (1); `model, CN_alpha, Cm_alpha, h_n` (2); `component, Re, Cf, FF, S_wet, D_over_q` (3); `fineness, L, D, S_wet, FF, Cf, D_over_q` (4); `alpha_deg, CN_potential, CN_allen_perkins` (5).

## 5. Validate and improve

| quantity | closed form / exact | numerical | difference | cause |
|---|---|---|---|---|
| spheroid couple at 4°, 55 m/s | 2qα·Vol = 1147.31 N·m | 1147.31 N·m at x_ref = 0, 3.5, 7 m | < 10⁻⁵ | cosine grid |
| paraboloid C_Nα on base; centre of pressure | 2; L/2 = 1.0000 m | 2.0000; 1.0000 m | < 10⁻⁶ | none |
| spheroid moment, L/D = 6.36 | exact: 0.882 × slender | slender body | +13 % | finite fineness |
| trainer Cm_α | 2Vol/(Sc) = 0.574 /rad | 0.561 /rad | −2.3 % | small open base |

The engineering answer: **the fuselage contributes a destabilising Cm_α between +0.22 and +0.49 /rad about the wing quarter chord, moving the neutral point forward by 0.048–0.093 of the chord (6–11 cm), from 0.505 c to between 0.457 c and 0.412 c.** The bracket comes from the tail cone's boundary layer, not from the numerics. To narrow it: add the wing's upwash and downwash (Extension 1), which raises the forward fuselage's contribution; and, best of all, compare with a wind-tunnel test of the wing–body combination, the validation route of `P10.22`.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (1)–(4) and of the paraboloid's centre of pressure; `fuselage_slender_body.png`; the five experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) giving the neutral-point shift with its bracket and the assumption that sets it.

### Adversarial self-checks

- **What if the fuselage is longer?** Its volume grows and so does Cm_α; Problem 4 lengthens it to 7.75 m: the corrected Cm_α becomes 0.557 /rad and the neutral point moves to 0.401 c.
- **What if α → 0?** Everything is linear in α: the moment and its derivative stay finite and the derivative is what matters. At large α, (6) peaks at 45° and (7) adds the cross-flow drag.
- **Which assumption fails first in the real world?** Attached flow on the tail cone, and the free-stream α on a fuselage that sits in the wing's flow field.
- **What measurement could falsify the model?** A wind-tunnel pitching-moment curve of the fuselage alone and of the wing–body combination: the fuselage-alone slope tests (3)–(6), the difference with the wing tests interference.
- **Does the result survive a change of grid or method?** Grid refinement changes the integrals by less than 10⁻⁵; the change of loading model changes the answer by a factor of 2.
- **Which decision changes if the uncertainty doubles?** Tail sizing: a shift uncertain by ±0.05 c instead of ±0.025 c would have to be covered by a larger tail or a narrower centre-of-gravity range.

## Check yourself

1. What is the apparent mass per unit length of a circular section moving sideways in potential flow?
2. Write the slender-body loading dN/dx and state where on a body it is positive and negative.
3. Why does a closed body in potential flow carry no lift, and what force does it feel instead?
4. Why is Munk's moment independent of the moment reference point?
5. What is C_Nα of a slender body that ends in a full section, and on which area is it based?
6. Why does slender-body theory over-predict the moment of a spheroid of fineness 6, and by roughly how much?
7. Which physical effect does the Allen–Perkins cross-flow term add, and how does it scale with α?
8. Does a fuselage move the neutral point forward or aft, and what has to be enlarged to compensate?
