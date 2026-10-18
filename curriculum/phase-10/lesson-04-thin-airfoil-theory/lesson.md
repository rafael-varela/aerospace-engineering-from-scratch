# P10.04  Thin Airfoil Theory

> Replace the airfoil by a vortex sheet on its chord, make the camber line a streamline, and three Fourier coefficients of the camber slope give the lift, the zero-lift angle and the pitching moment.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the thin-airfoil integral equation, solve it with the Glauert substitution as a Fourier series, and obtain c_l = 2π(α − α_L0), α_L0, c_m,c/4 = (π/4)(A₂ − A₁) and the centre of pressure. 2. Treat a plain flap as a camber line and derive its effectiveness τ = 1 − (θ_h − sin θ_h)/π and its moment. 3. Implement and test a thin-airfoil solver for any camber-line slope, with Gauss–Legendre quadrature split at kinks, and check it against closed forms, a discrete-vortex model with second-order convergence, and the vortex panel method as the thickness goes to zero. 4. Use it to size the trainer's cruise attitude and a landing flap, and judge where the linear theory stops. |
| **Prerequisites** | `P02.09` (Fourier series, orthogonality of cosines), `P10.03` (circulation, Kutta condition, L′ = ρVΓ), `P10.02` (the NACA camber lines), `P03.04` (Gauss–Legendre quadrature), `P03.11` (observed order), `P01.13` (small linear systems, for the discrete-vortex check). Skill check: you can compute a Fourier cosine coefficient of a given function by integration. |
| **Engineering question** | At what angle of attack does the trainer's NACA 2412 section give its cruise c_l of 0.331, what is its pitching moment, and how much deflection of a 25 % chord plain flap gives the approach c_l of 1.112 at α = 4°? Thin-airfoil answers: α_L0 = −2.077°, cruise at α = 0.94°, c_m,c/4 = −0.0531 at every α; the flap needs δ = 6.7° (τ = 0.609), and the trimmed section then has c_m,c/4 = −0.129. The panel method (thickness 12 %) gives 11 % more lift at the same angle. |
| **Assumptions and limits** | Inviscid, incompressible, 2-D; thin section (thickness ignored) with small camber and small angles, so the boundary condition is applied on the chord and slopes replace angles; attached flow with the Kutta condition at the trailing edge. Fails for thick sections (lift slope above 2π, `P10.03`, `P10.05`), for large flap deflections and near stall (separated flow, `P10.08`, `P10.15`), and gives no drag. |
| **Analytical work** | The vortex-sheet integral equation; the Glauert substitution and integral; the Fourier coefficients A₀, Aₙ; c_l, α_L0, c_m,c/4, c_m,le, x_cp; the parabolic camber line in closed form; the plain flap's τ and dc_m/dδ; the small-tab limit τ ≈ 4√E/π. |
| **Computational deliverable** | `starter/main.py`: `theta_of_x`, `glauert_integrals`, `fourier_coefficients`, `thin_airfoil`, `flap_camber_slope`, `flap_effectiveness`, `flap_cm_derivative`, `lumped_vortex` (the NACA camber-line slopes, `combine_slopes`, `observed_order`, the compact panel method `naca4_coords` + `panel_cl_cm` copied from `P10.05`, and `main` are given). Tests prove: flat plate exact; parabolic camber closed forms; slope 2π and constant c_m,c/4; c_l(α_L0) = 0; x_cp moves forward with c_l and c_m,le = c_m,c/4 − c_l/4; split quadrature equal to adaptive integration to 10⁻¹²; NACA 2412 α_L0 against a hand-integrated closed form (−2.077°); τ and dc_m/dδ at four flap sizes; flap limits; superposition; the 230 line's ideal c_l of 0.3; the lumped-vortex model exact for plate and parabola and second order for the 2412; panel-method lift tending to thin-airfoil lift in proportion to thickness. |
| **Practical activity** | **computational experiment**: quadrature with and without the split; the trainer's lift line and flap sizing; a camber-position sweep; the discrete-vortex convergence study; thin-airfoil theory against the panel method as thickness → 0. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `thin_airfoil.png`, the derivations of §2.2–2.5, the experiment tables, and a half-page note on the trainer's cruise attitude and flap with the theory's limits. |
| **Assessment** | `assessments/phase-10/grading-rubric.md`; competencies 1–4. Pass: the Fourier coefficients and c_l, α_L0, c_m,c/4 derived; all public tests pass; the discrete-vortex order and the thickness trend reported; the note says which number thin-airfoil theory gets wrong first and by how much. |
| **References** | `@mit-ocw-16-01-unified-engineering-fluids` 16.01 Fluids (spring) Lecture 2 "Airfoil Vortex Sheet Models, Thin Airfoil Analysis Problem", Lecture 3 "Thin Airfoil Analysis Problem (cont.)", Lecture 4 "Thin Airfoil Theory Application: Analysis Example" (https://ocw.mit.edu/courses/16-01-unified-engineering-i-ii-iii-iv-fall-2005-spring-2006/resources/f02_sp/, …/f03_sp/, …/f04_sp/) (required-free). `@mit-ocw-16-100-aerodynamics` Lecture 15 "Thin Airfoil Theory Summary", Lecture 16 "Important Concepts in Thin Airfoil Theory" (https://ocw.mit.edu/courses/16-100-aerodynamics-fall-2005/resources/16100lectre15_cg/, …/16100lectre16_cg/) (required-free). `@leishman-intro-aerospace-flight-vehicles` "Classic Airfoil Theory" (https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/classic-airfoil-theory/): Thin Airfoil Assumption; Theory of Vortex Sheets; Thin Airfoil Theory (required-free, link only). `@ntrs-naca-tr-824-summary-of-airfoil-data` "Theoretical Considerations": Angle of Zero Lift; Moment Coefficients (https://ntrs.nasa.gov/citations/19930090976) (supplementary-free: NACA's own numerical methods and examples). Optional textbooks for deeper study: J. D. Anderson, *Fundamentals of Aerodynamics*, and I. H. Abbott and A. E. von Doenhoff, *Theory of Wing Sections* (both paid; libraries), the classic derivations with many worked camber lines. |
| **Extensions** | ★ Compute the chordwise load Δcp(x) = 2γ/V from the Fourier series and compare with the panel method's upper-minus-lower Cp. ★ Add a leading-edge flap (droop) and show it changes c_l,ideal but hardly α_L0. ★ Invert the theory: given a target load, find the camber line (the start of `P10.07`). |
| **Time** | unmeasured |

## 1. The engineering problem

An aircraft designer needs three numbers from an airfoil before anything else: the angle at which it gives no lift, how fast the lift grows with angle, and the pitching moment it brings to the trim of the aeroplane. Thin-airfoil theory, worked out in the 1920s, gives all three from the camber line alone, with a pencil. It is still the first tool for camber-line design, flap sizing and the quick check of any section code, and it underlies the lifting-line and vortex-lattice methods of `P10.10` and `P10.11`.

For the reference trainer (NACA 2412, chord 1.2 m) the questions are concrete. At cruise (55 m/s) it needs c_l = 0.331; on approach (30 m/s), 1.112. At what angle does the section give the cruise value, what nose-down moment must the tail balance, and how much plain flap would let it reach the approach value without a nose-high attitude? (These are section numbers; the finite wing needs more, `P10.09`.)

Success criteria, written first: the solver reproduces the flat plate and the parabolic camber line exactly; c_l is linear in α with slope 2π and c_m,c/4 is independent of α; the flap reproduces its closed form; an independent discrete-vortex model converges to the same answers at second order; and the panel method approaches thin-airfoil theory as the thickness goes to zero.

## 2. Understand and derive

### 2.1 The model

Replace the airfoil by a **vortex sheet** of strength γ(x) per unit length lying on the chord, from x = 0 to c. By `P10.03`, each element γ dx is a point vortex, and the total Γ = ∫γ dx gives the lift L′ = ρVΓ. The sheet must make the camber line z(x) a streamline. For a thin section at a small angle, the velocity normal to the camber line is V(α − dz/dx) from the free stream plus the sheet's own induced velocity, and the boundary condition can be applied on the chord instead of on the camber line (an error of second order in the camber). With γ clockwise positive, an element at ξ induces a downwash γ dξ/(2π(x − ξ)) at x, so

**(1/2π) ∫₀^c γ(ξ) dξ/(x − ξ) = V(α − dz/dx),**   (1)

plus the **Kutta condition** γ(c) = 0: the two surfaces have the same speed at the trailing edge, so there is no load there.

### 2.2 The Glauert substitution

Put x = (c/2)(1 − cos θ), so θ runs from 0 (leading edge) to π (trailing edge). The solution of (1) that satisfies the Kutta condition is

**γ(θ) = 2V [A₀ (1 + cos θ)/sin θ + Σₙ Aₙ sin nθ].**   (2)

The first term is the flat plate's load, infinite at the leading edge (the suction peak of a sharp nose) and zero at the trailing edge; the series adds the camber's load. Substituting (2) into (1) with the Glauert integral ∫₀^π cos nφ dφ/(cos φ − cos θ) = π sin nθ/sin θ turns (1) into dz/dx = (α − A₀) + Σ Aₙ cos nθ: a Fourier cosine series (`P02.09`) for the camber slope. Orthogonality gives

**A₀ = α − (1/π)∫₀^π (dz/dx) dθ,  Aₙ = (2/π)∫₀^π (dz/dx) cos nθ dθ.**   (3)

Only A₀ depends on α. The details are Problem 2.

### 2.3 Lift, zero-lift angle, moment

Integrating (2): Γ = cV(πA₀ + πA₁/2), so

**c_l = π(2A₀ + A₁) = 2π(α − α_L0),  α_L0 = −(1/π)∫₀^π (dz/dx)(cos θ − 1) dθ.**   (4)

The **lift slope is 2π per radian** for every thin section; camber only shifts the line. The moment about the leading edge is c_m,le = −(c_l/4 + (π/4)(A₁ − A₂)); moved to the quarter chord,

**c_m,c/4 = (π/4)(A₂ − A₁),  x_cp/c = ¼ [1 + π(A₁ − A₂)/c_l].**   (5)

c_m,c/4 depends only on the camber line: **the quarter chord is the aerodynamic centre**, as `P10.01`'s synthetic data assumed. The centre of pressure moves forward towards c/4 as c_l grows and runs off to infinity as c_l → 0: it is a poor reference point, which is why designers quote c_m,ac.

**Worked example: a parabolic camber line.** z = 4h x(1 − x) has dz/dx = 4h(1 − 2x) = 4h cos θ, a single Fourier mode: A₀ = α, A₁ = 4h, Aₙ = 0 for n ≥ 2. Then α_L0 = −2h, c_m,c/4 = −πh and the ideal c_l (A₀ = 0) is 4πh. For h = 0.025: α_L0 = −2.865°, c_m,c/4 = −0.0785, and at α = 3°, c_l = 0.643 with the centre of pressure at 0.372c.

### 2.4 The NACA 2412 and the 230 lines

For the 4-digit line (`P10.02`) the slope is piecewise linear in x with a kink at p, so the integrals (3) have closed forms on each piece, and numerically the quadrature must be split at θ_p = arccos(1 − 2p). With Gauss–Legendre points on each smooth piece the integrals converge to round-off with 8 points per piece; across the kink with 64 points the error is still 2.5 × 10⁻⁶ (§4). For the 2412: **A₁ = 0.08150, A₂ = 0.01386, α_L0 = −2.077°, c_m,c/4 = −0.05312.** The 230 line gives α_L0 = −1.094°, c_m,c/4 = −0.01284 and ideal c_l 0.300 (its design value): forward camber, a quarter of the pitching moment.

### 2.5 A plain flap is a camber line

A flap of chord E·c deflected δ (trailing edge down) on a flat plate is a camber line with slope 0 ahead of the hinge x_h = 1 − E and −δ behind it. Integrating (3) from θ_h = arccos(1 − 2x_h) to π:

**τ ≡ −dα_L0/dδ = 1 − (θ_h − sin θ_h)/π,  dc_m,c/4/dδ = −½ sin θ_h (1 − cos θ_h).**   (6)

So Δc_l = 2πτδ. A 25 % flap has τ = 0.609: every degree of flap is worth 0.61° of incidence. A small tab still has τ ≈ 4√E/π, 0.13 for a 1 % tab, because a deflection at the trailing edge, where the Kutta condition acts, changes the circulation of the whole section. Being linear, flap and camber superpose.

### 2.6 Independent checks

**A discrete-vortex model.** Cut the chord into n panels; put a point vortex at each panel's quarter point and a collocation point at its three-quarter point, and enforce tangency there:

Σⱼ −Γⱼ / (2π(x_ci − x_vj)) = V(dz/dx(x_ci) − α).   (7)

The "1/4–3/4 rule" makes a single panel exact for a flat plate (c_l = 2πα for any n), and the model is also exact in c_l for a parabola; for the 2412 it converges at second order (§4). It shares the linearisation of (1) but not the Fourier solution.

**The panel method.** `P10.05`'s linear-strength vortex panel method solves the full potential flow round the real thick section. A compact copy is given here. As the thickness goes to zero its c_l must tend to (4).

### 2.7 Where the model stops

- **Thickness.** Adds lift slope: 11 % at 12 % thickness (§4), and a more negative c_m. Thin-airfoil theory's 2π is a lower bound for the inviscid section; viscosity pulls the real slope back down (`P10.03`, Problem 6).
- **Large angles and flap deflections.** Slopes are no longer small and, more important, the boundary layer separates at the hinge or the nose. Thin-airfoil theory predicts Δc_l = 2.0 for a 30° plain flap; no plain flap achieves that (`P10.08`).
- **The leading edge.** (2) is infinite at x = 0; real noses have a finite suction peak, which sets leading-edge stall (`P10.15`).
- **No drag.** Ever.

## 3. Build the computational model

Open `starter/main.py`. The camber-line slopes, `combine_slopes`, `observed_order` and the compact panel method are given.

### Step 1: the integrals

`glauert_integrals(dzdx, breaks, n_max, n_gauss)` computes I_n = ∫₀^π (dz/dx) cos nθ dθ for n = 0..n_max with Gauss–Legendre points on each interval between the breaks (converted to θ):

```python
edges = [0.0] + sorted(theta_of_x(b) for b in breaks if 0 < b < 1) + [np.pi]
for lo, hi in zip(edges[:-1], edges[1:]):
    th = 0.5 * (hi - lo) * (xg + 1) + lo
    s = dzdx(0.5 * (1 - np.cos(th)))
    for n in range(n_max + 1):
        I[n] += np.sum(0.5 * (hi - lo) * wg * s * np.cos(n * th))
```

`fourier_coefficients` is (3) from those integrals.

### Step 2: the coefficients

`thin_airfoil(dzdx, alpha, breaks)` returns A, c_l, α_L0, c_m,c/4, c_m,le and x_cp from (4) and (5).

### Step 3: flaps

`flap_camber_slope(E, delta)` returns the slope function and the hinge as a break; `flap_effectiveness` and `flap_cm_derivative` are (6).

### Step 4: the discrete-vortex check

`lumped_vortex(dzdx, alpha, n, spacing)` builds the n × n influence matrix of (7) and solves it; c_l = 2ΣΓ/(Vc), c_m,c/4 = −2ΣΓⱼ(x_vj − ¼)/(Vc²).

Running `main()` prints

```
NACA 2412 camber line: A1 = 0.08150, A2 = 0.01386, alpha_L0 = -2.077 deg, cm_c/4 = -0.05312
trainer cruise cl = 0.3309 at alpha = 0.940 deg (section); x_cp/c = 0.4105; approach cl = 1.1122 needs alpha = 8.07 deg without flap
   alpha =  4.0 deg: cl = 0.6664, cm_c/4 = -0.05312, x_cp/c = 0.3297
NACA 230 camber line: alpha_L0 = -1.094 deg, cm_c/4 = -0.01284, ideal cl (A0 = 0) = 0.3000
   E = 0.25: tau = 0.6090, dcm_c/4/d delta = -0.6495 /rad
approach cl 1.1122 at alpha = 4.0 deg with a 25 % plain flap: delta = 6.68 deg (check cl = 1.1122); cm_c/4 = -0.1288
approach cl 1.1122 at alpha = 6.0 deg with a 25 % plain flap: delta = 3.39 deg (check cl = 1.1122); cm_c/4 = -0.0916
   flat plate, n = 7: cl = 0.4386490845 (2 pi alpha = 0.4386490845), cm_c/4 = 3.5e-17
```

(abridged; the quadrature, flap, lumped-vortex and panel tables of §4 are printed too) and writes `thin_airfoil.png`: the lift lines with 0°, 10° and 20° of flap, τ and the flap moment against E, and the discrete-vortex convergence.

### Verify

`python -m pytest curriculum/phase-10/lesson-04-thin-airfoil-theory` runs 16 tests (listed in the metadata).

### Compare

Three independent routes: the closed forms (flat plate, parabola, hand-integrated 2412 α_L0, flap), the discrete-vortex model, and the panel method at small thickness. TR-824's "Theoretical Considerations" tabulate NACA's own α_L0 and c_m values for the standard camber lines; compare yours with them.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Predict, then run.

1. **Split the quadrature.** Predict the α_L0 error of an unsplit 64-point rule. *Finding:* with the split, 3.0 × 10⁻⁸ with 4 points per piece and round-off with 8; without it, 3.3 × 10⁻⁵ (8 points), 2.7 × 10⁻⁵ (16), 3.7 × 10⁻⁶ (32), 2.5 × 10⁻⁶ (64): a kink costs four orders of magnitude.
2. **The trainer.** *Finding:* cruise c_l 0.331 at α = 0.940° (x_cp = 0.411c); c_l = 0.228, 0.447, 0.666, 1.105 at 0°, 2°, 4°, 8°, with c_m,c/4 = −0.0531 throughout and x_cp moving from 0.483c to 0.298c. Approach without flap: α = 8.07°. With a 25 % plain flap: δ = 6.68° at α = 4° (c_m,c/4 = −0.129), or 3.39° at α = 6° (−0.092). The flap buys attitude and costs pitching moment, which the tail must pay.
3. **Camber position** (Problem 4). *Finding:* for m = 0.02, α_L0 = −1.80°, −1.92°, −2.08°, −2.29°, −2.59° and c_m,c/4 = −0.037, −0.045, −0.053, −0.063, −0.075 for p = 0.2–0.6; both exactly proportional to m.
4. **Discrete vortices.** Predict the order on uniform panels. *Finding:* c_l errors 1.46 × 10⁻⁴ (n = 10) to 1.56 × 10⁻⁷ (n = 320), orders 1.93, 1.97, 1.98, 1.99, 2.00; c_m errors 6.07 × 10⁻⁴ to 5.91 × 10⁻⁷, order 2.00. On cosine spacing the c_l errors jump about (from 2.7 × 10⁻⁴ to 2.2 × 10⁻⁷, orders −5.1 to 5.3) because the kink at x = 0.4 falls inside a panel (Problem 5). The flat plate is exact with n = 7 (c_l = 0.4386490845 = 2π × 4°), the parabola too.
5. **Thickness → 0.** Predict the panel method's lift excess at 12 % thickness. *Finding:* (2 % camber at 40 %, α = 4°, 240 panels) c_l = 0.74160, 0.70303, 0.68411, 0.67474 for t/c = 0.12, 0.06, 0.03, 0.015, against 0.66644: +11.28 %, +5.49 %, +2.65 %, +1.24 %, halving with the thickness. c_m,c/4 = −0.0612, −0.0569, −0.0550, −0.0541 against −0.0531.

Raw data format: CSV `n_gauss, err_split, err_unsplit` (1); `alpha_deg, cl, cm_c4, x_cp` (2); `m, p, alpha_L0_deg, cm_c4` (3); `n, err_cl_uniform, err_cm_uniform, err_cl_cosine` (4); `t_over_c, cl_panel, cl_thin, cm_panel, cm_thin` (5).

## 5. Validate and improve

| quantity | closed form / independent route | solver | difference |
|---|---|---|---|
| flat plate c_l | 2πα | `thin_airfoil`; lumped vortex n = 7 | < 10⁻¹² |
| parabola α_L0, c_m | −2h, −πh | `thin_airfoil` | < 10⁻¹² |
| 2412 α_L0 | hand-integrated: −2.077° | −2.077° | < 10⁻¹⁰ relative |
| flap τ, E = 0.25 | (6): 0.6090 | quadrature | < 10⁻¹² |
| 2412 c_l at 4° | lumped vortex n = 320: 0.66644 | 0.66644 | 1.6 × 10⁻⁷, order 2 |
| 2412-camber c_l at 4° | panel, t → 0: 0.6747 at 1.5 % | 0.66644 | 1.24 %, ∝ t |

Every route agrees within its expected error, and the panel method's excess falls in proportion to thickness, the expected thin-airfoil limit. The engineering answer: **the trainer's NACA 2412 section reaches the cruise c_l of 0.331 at α = 0.94° (α_L0 = −2.077°) with c_m,c/4 = −0.0531; a 25 % plain flap deflected 6.7° reaches the approach c_l of 1.112 at α = 4°, at the cost of c_m,c/4 = −0.129.** Thickness adds about 11 % to the lift at a given angle (inviscid) and the boundary layer takes most of it back (`P10.03`, `P10.06`); the finite wing needs more angle (`P10.09`).

What to improve: the flap's effectiveness at real deflections (`P10.08`); the thick-section and viscous corrections (`P10.05`, `P10.06`). What not to improve: the quadrature, which is at round-off.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (3)–(6); `thin_airfoil.png`; the five experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) on the trainer's cruise attitude, flap and pitching moment, naming the first assumption to fail.

### Adversarial self-checks

- **What if the camber doubles?** α_L0 and c_m,c/4 double exactly (linearity); the lift slope does not change.
- **What if the flap chord → 0?** τ → 4√E/π, not zero as fast as the chord: trailing-edge devices (Gurney flaps, tabs) are powerful for their size.
- **Which assumption fails first?** Attached flow at a deflected flap's hinge and at the nose at high α; then the thin-section assumption for sections thicker than about 12 % (the panel method's +11 %).
- **What measurement could falsify it?** The zero-lift angle and the pitching moment at moderate Re are close to thin-airfoil values for thin sections; a measured c_m that changes strongly with α at small α would mean the a.c. is not at c/4 (thickness, viscosity), and a lift slope well above 2π would mean thickness matters.
- **Does the result survive a change of method or grid?** Yes: three independent routes, and a quadrature converged to round-off.
- **Which decision changes if the flap is 30 % less effective than τ?** The approach deflection at α = 4° rises from 6.7° to about 9.5°, still small; the decision that changes is the flap type, if the real section cannot keep the flow attached at the larger deflection.

## Check yourself

1. What is the Glauert substitution, and why does it help?
2. Write c_l in terms of the Fourier coefficients. What is the lift slope?
3. Write c_m,c/4 in terms of the coefficients. Why does it not depend on α, and what does that make the quarter chord?
4. Write α_L0 as an integral. What is it for the NACA 2412 camber line?
5. How is a plain flap modelled in thin-airfoil theory, and what is its effectiveness for a 25 % flap?
6. Name two independent checks of a thin-airfoil solver, and what each should show.
7. Why is the quadrature split at x = p for a NACA 4-digit camber line?
