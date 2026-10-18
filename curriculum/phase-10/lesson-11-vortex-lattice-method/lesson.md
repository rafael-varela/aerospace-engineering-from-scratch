# P10.11  Vortex Lattice Method

> Cover the wing with horseshoe vortices, make the flow tangent at one point on each panel, and solve one linear system: the vortex lattice gets lift, induced drag, span load and pitching moment for swept, tapered and dihedral wings that lifting-line theory cannot handle.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the vortex lattice method: horseshoe vortices with bound legs on the panels' quarter-chord lines, control points at three-quarter chord (and why that 1/4–3/4 rule reproduces 2π), the influence-coefficient matrix and the flow-tangency system. 2. Implement it from scratch for straight-tapered wings with sweep and dihedral, with lift and moment from Kutta–Joukowski on the bound legs, induced drag in the Trefftz plane, the span load and the aerodynamic centre. 3. Run a lattice-convergence study, measure the observed order, and explain why Hough's quarter-panel tip inset raises it from 1 to 2. 4. Verify the lattice against lifting-line theory for unswept wings (a difference that falls like 1/AR), the 2-D limit, symmetry, and the Helmbold–DATCOM formula for swept wings. 5. Report the trainer wing's CL_α, e, CDi and aerodynamic centre. |
| **Prerequisites** | `P10.09` (Biot–Savart, horseshoes, Trefftz plane), `P10.10` (lifting line, e, τ), `P03.08` (direct linear solvers: the dense LU solve and its cost), `P03.11` (observed order without an exact solution), `P10.01` (moment transfer and aerodynamic centre), `P04.04` (rotations, for dihedral). Skill check: you can write the normal velocity of a vector field on a tilted panel. |
| **Engineering question** | What lift-curve slope, span efficiency, induced drag and aerodynamic centre does the reference trainer's rectangular wing have when its chord is resolved, and how much do the lifting-line answers of `P10.10` change? Answer: CL_α = 4.631 per radian (0.0808 per degree), 5.1 % below the lifting line's 4.878; e = 0.970 (lifting line 0.934); at cruise (CL 0.331, α = 2.02°) CDi = 0.00431, an induced drag of 95.7 N; the aerodynamic centre sits at 24.3 % of the chord. The lattice converges at observed order 2 with Hough's tip inset; its difference from the lifting line falls as −0.43/AR. |
| **Assumptions and limits** | Inviscid, incompressible, steady, attached flow; thin lifting surface (thickness ignored; camber only as a zero-lift-angle shift); small angles (boundary conditions on the mean surface, trailing legs along the x axis rather than the local flow); flat wake; linear in α, so no stall and no leading-edge or tip vortex lift (which matters below AR ≈ 2). Dense LU solve: O(N³) cost, fine for a few thousand panels. |
| **Analytical work** | The 1/4–3/4 rule from a single vortex on a flat plate; the flow-tangency condition on a panel; Kutta–Joukowski on a bound segment; the Trefftz-plane drag of trailing vortices in the y–z plane; moment transfer and x_ac = x_ref − (dCm/dCL)c_ref; the observed order from three solutions without an exact answer; the Helmbold–DATCOM estimate. |
| **Computational deliverable** | `starter/main.py`: `make_lattice`, `horseshoe_unit_velocity`, `solve_vlm`, `vlm_forces`, `analyse`, `lift_slope_and_ac`, `observed_order` (`segment_velocity` and `semi_infinite_velocity` from `P10.09`, a compact `P10.10` lifting line `lifting_line_CL_e`, `helmbold_datcom`, `mac` and `main` are given). Tests prove: the lattice tiles the planform with upward normals and its points where the 1/4–3/4 rule puts them; Hough's inset; a unit horseshoe's far wake is 2-D; a strip in the middle of a very long wing gives 2π and its aerodynamic centre the quarter chord; the span load is symmetric; lift ∝ sin α and drag ∝ lift²; zero lift at α_L0; observed order 2 (Hough) and 1 (plain) in span, 2 in chord; VLM below lifting line by about 0.43/AR; swept slopes follow the Helmbold–DATCOM trend; moment transfer is exact and dCm/dCL = 0 about x_ac; dihedral and taper act the right way. |
| **Practical activity** | **computational experiment**: lattice convergence with observed orders for three spanwise spacings and in the chordwise direction; VLM against lifting line for AR 4–50; swept wings against the Helmbold–DATCOM formula, with a check against measured lift slopes (NACA TN-1669) from `curriculum/phase-10/data/`; the trainer wing's coefficients. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `vlm.png`, the derivation of the 1/4–3/4 rule and of the Trefftz drag, the convergence table with observed orders, the VLM–lifting-line table, the swept-wing table, and a half-page note on the trainer's wing. |
| **Assessment** | `assessments/phase-10/grading-rubric.md`; competencies 1–4. Pass: the 1/4–3/4 rule is derived; all public tests pass; the observed orders are reported and Hough's inset explained; the note states the trainer's CL_α, e and x_ac and compares them with the lifting line. |
| **References** | `@ntrs-nasa-tn-d-6142-vortex-lattice-program` (Margason & Lamar, https://ntrs.nasa.gov/citations/19710007949): Basic Concepts and Limitations; Program Description, Part I "Geometry Computation" (Section 3. Horseshoe Vortex Lattice), Part II "Vortex Strength Computation", Part III "Aerodynamic Computation"; Effect of Vortex-Lattice Arrangement on Computed Aerodynamic Characteristics (supplementary-free, read first). `@ntrs-nasa-sp-405-vortex-lattice-utilization` (https://ntrs.nasa.gov/citations/19760021075): paper 1 "Historical Evolution of Vortex-Lattice Methods" (DeYoung); paper 17 "Arrangement of Vortex Lattices on Subsonic Wings" (DeJarnette); paper 18 "Lattice Arrangements for Rapid Convergence" (supplementary-free): where the tip inset of this lesson comes from. `@leishman-intro-aerospace-flight-vehicles` "Lifting Line & Finite Wing Theory" (https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/lifting-line-theory/): Introduction to Lifting Surface Theory (required-free). `@drela-avl` AVL User Primer "Vortex-Lattice Modeling Principles"; "Vortex Lattice Spacing Distributions" (https://web.mit.edu/drela/Public/web/avl/avl_doc.txt) (supplementary-free): the established VLM tool and its spacing advice. `@ntrs-naca-tn-1669-untapered-wings-aspect-ratio-sweep` (Goodman & Brewer, https://ntrs.nasa.gov/citations/19930082295) (required-free): the measured lift slopes used in part 5. Optional textbooks for deeper study (optional-paid, not required; library or used copy): Katz & Plotkin, *Low-Speed Aerodynamics*, for the VLM in full generality (camber, ring vortices, unsteady wakes); Bertin & Cummings, *Aerodynamics for Engineers*, for a worked four-horseshoe example. |
| **Extensions** | ★ Add camber: tilt each panel's normal by the camber slope at its three-quarter chord and recover the 2412's α_L0. ★ Add a horizontal tail as a second lifting surface 4.5 m behind and compute the trainer's neutral point (with downwash), the bridge to static stability. ★ Use `scipy.linalg.lu_factor` once and reuse it for many α and for the transposed (adjoint) system, and time it against `np.linalg.solve` (`P03.08`). |
| **Time** | unmeasured |

## 1. The engineering problem

Lifting-line theory puts all of a wing's lift on one line. That works for straight wings of high aspect ratio, but aircraft designers need more: swept wings (all jets), tapered and twisted wings with dihedral, tails, canards and the interaction between them, and the *chordwise* position of the lift, which sets the pitching moment and the aerodynamic centre, the reference point of every stability calculation. The **vortex lattice method** (VLM) spreads the lift over the planform as a grid of horseshoe vortices. It is the workhorse of conceptual design and stability-derivative estimation (AVL is the best-known free implementation), and it takes one linear system.

This lesson builds one from the Biot–Savart law of `P10.09`. Because a lattice is an approximation with two resolutions (chordwise and spanwise), it must be verified like any numerical method: by a convergence study with an observed order, by limiting cases (the 2-D strip, the infinite wing), by agreement with an independent method where both are valid (lifting line for straight high-aspect-ratio wings), and against a semi-empirical formula and measured data for swept wings.

Success criteria, written first: a single chordwise panel on a very long wing gives the 2-D slope 2π to 1 %; the lattice converges with a measured order; the VLM and lifting-line lift slopes agree for straight wings with a difference that vanishes as AR → ∞; the swept-wing slopes follow the Helmbold–DATCOM trend within 3 %; the moment about the computed aerodynamic centre does not change with α.

## 2. Understand and derive

### 2.1 The 1/4–3/4 rule

Model a 2-D flat plate of chord c at angle α by a single point vortex Γ at a distance a behind the leading edge, and require the flow to be tangent to the plate at one control point at distance d behind the vortex. The vortex induces a downwash Γ/(2πd) there, which must cancel the free stream's normal component V∞α: Γ = 2πdV∞α. The lift is ρV∞Γ, so cl = 2Γ/(V∞c) = 4πdα/c. Thin-airfoil theory says cl = 2πα and that the flat plate's centre of pressure is at the quarter chord, so **d = c/2 and a = c/4**: the vortex at the quarter chord, the control point at the three-quarter chord. Applied to every panel, this **1/4–3/4 rule** makes each chordwise strip reproduce thin-airfoil theory's lift and moment exactly in the 2-D limit; it is what makes a coarse lattice useful.

### 2.2 The lattice and the linear system

Divide the planform into nx chordwise × ny spanwise panels. Each panel j carries a horseshoe of unknown strength Γ_j: a bound leg along its own quarter-chord line from A_j to B_j (left to right, so positive Γ lifts) and two trailing legs to x = +∞ (`P10.09`). At each panel's control point P_i (three-quarter chord, mid-span of the strip) the normal velocity must vanish:

Σ_j AIC_ij Γ_j = −V∞·n_i,  AIC_ij = v_j(P_i)·n_i,   (1)

where v_j is the velocity induced by a unit horseshoe j and n_i the panel's unit normal (tilted by the dihedral). With V∞ = (cos α, 0, sin α), (1) is an N × N dense linear system (N = nx·ny), solved by LU decomposition (`P03.08`) at a cost of N³/3 operations; for N = 240 that is a few milliseconds. Camber would enter as a tilt of each n_i; here the surface is flat, and a section's camber is represented by shifting α to α − α_L0.

### 2.3 Forces and moments

The force on bound leg j is Kutta–Joukowski, F_j = ρΓ_j V∞ × l_j, l_j = B_j − A_j. Lift is its component normal to V∞; the pitching moment about (x_ref, 0, 0), positive nose-up, is M = Σ[(x_ref − x_j)F_z,j + z_j F_x,j]. With q = ½ρV∞², CL = L/(qS) and Cm = M/(qS c_ref), c_ref the mean aerodynamic chord (MAC) (2/3)c_r(1 + λ + λ²)/(1 + λ). Moving the reference point by Δx changes Cm by CL cos α Δx/c_ref (exactly, a test). The **aerodynamic centre** is the reference point about which Cm does not change with CL:

**x_ac = x_ref − (dCm/dCL) c_ref.**   (2)

For a wing alone it is also the wing's neutral point.

### 2.4 Induced drag in the Trefftz plane

Using the bound-leg force with the *induced* velocity added would give induced drag too, but it converges poorly. The far-field route of `P10.09` is better: far downstream the trailing legs at the strip edges (y_k, z_k) are infinite 2-D vortices of strength s_k = G_left − G_right, where G is the strip circulation (the sum of its panels' Γ). They induce a velocity v in the y–z plane, and

**D_i = −½ρ Σ_strips G_j (v·n_j) Δs_j,**   (3)

with n_j the strip's normal in the y–z plane and Δs_j its width there (dihedral included). Then e = CL²/(πAR CDi).

### 2.5 Tips and convergence: Hough's inset

A continuous span load falls to zero at the tip like √(distance). A lattice that puts its last trailing vortex exactly at the geometric tip effectively concentrates the tip's loss of lift there; the error per strip near the tip is first order in the strip width, and the whole solution converges at **order 1**. Hough (summarised in SP-405's "Lattice Arrangements for Rapid Convergence", `@ntrs-nasa-sp-405-vortex-lattice-utilization`) showed that insetting the lattice by a quarter of a panel width at each tip places the tip vortex where the continuous load's vorticity centroid is, and restores **order 2**. With uniform strips of width d = b/(ny + ½) the lattice spans b − d/2. The code offers uniform, cosine and Hough spacings, and part 4 measures all three. Because no exact solution exists for a rectangular wing, the order is measured from three successive lattices: p = ln[(f₁ − f₂)/(f₂ − f₃)]/ln 2 (`P03.11`).

### 2.6 Swept wings and the Helmbold–DATCOM estimate

On a swept wing the downwash of the (swept) bound vortices and of the inboard trailing legs loads the outboard wing: the lift moves towards the tips and aft. The lift slope falls roughly with the cosine of the sweep. A widely used semi-empirical estimate, a fit to lifting-surface results and experiments (the Helmbold–DATCOM form):

**CL_α ≈ 2πA / (2 + √(A²(1 + tan²Λ) + 4))**,   (4)

with Λ the mid-chord sweep (for untapered wings every chord line has the same sweep). It is a formula, not a law, and its agreement with a lattice is a consistency check on the *trend*, with the absolute level differing by a few per cent.

### 2.7 Where the method stops

No thickness (pressure distributions are not surface pressures), no viscosity (no drag beyond induced, no stall), a wake along x instead of along the local flow (fine at small α), and no separated vortex lift (the leading-edge and side-edge vortices of low-aspect-ratio wings, which add lift nonlinearly below AR ≈ 2). The VLM is a linear lifting-surface theory: excellent for slopes, loads, moments and induced drag of attached flow, and silent about everything else.

## 3. Build the computational model

Open `starter/main.py`. Given: `segment_velocity`, `semi_infinite_velocity` (from `P10.09`, unit strength), `lifting_line_CL_e` (compact `P10.10`), `helmbold_datcom`, `mac`, `main`.

### Step 1: the lattice

`make_lattice(b, c_root, taper, sweep_c4, dihedral, nx, ny, spacing)` builds for every panel the bound leg A → B on its quarter-chord line, the control point at three-quarter chord of the strip's mid-span section, the normal (from the cross product of the panel's diagonals) and the area. The leading edge is placed from the quarter-chord sweep: x_LE(y) = c_r/4 + |y| tan Λ − c(y)/4.

### Step 2: influence and solve

`horseshoe_unit_velocity(P, A, B)` = bound segment + leg from B along +x − leg from A along +x, all panels at all points at once (an array of shape (N, N, 3)). `solve_vlm` forms AIC_ij = v·n_i with `np.einsum("ijk,ik->ij", ...)` and solves (1).

### Step 3: forces

`vlm_forces` is section 2.3 and (3); `analyse` chains them with x_ref at the root quarter chord and c_ref = MAC; `lift_slope_and_ac` differences ±2° for CL_α and (2); `observed_order` is section 2.5.

```python
F = gamma[:, None] * np.cross(np.broadcast_to(V, l.shape), l)        # rho = V = 1
My = np.sum((x_ref - xm[:, 0]) * F[:, 2] + xm[:, 2] * F[:, 0])
```

Running `main()` prints (abridged):

```
reference trainer wing (rectangular, AR 8.333), VLM nx = 6, ny = 40 (Hough spacing):
  CL_alpha = 4.6305 /rad (0.08082 /deg); lifting line 4.8779 (-5.07 %)
  aerodynamic centre x_ac = 0.2913 m = 0.2427 MAC behind the LE
  cruise CL 0.3309 at alpha 2.017 deg: CDi = 0.004306, e = 0.9702 (lifting line e = 0.9339); induced drag 95.7 N
lattice convergence, trainer wing, alpha = 5 deg: CL for ny = 10, 20, 40, 80 (nx = 4)
  uniform : CL [0.425492 0.41558  0.409914 0.406903], observed orders [0.81, 0.91]; e [1.0697 1.0186 0.994  0.982 ]
  cosine  : CL [0.422075 0.413589 0.408845 0.406353], observed orders [0.84, 0.93]; e [1.0698 1.0193 0.9945 0.9823]
  hough   : CL [0.400471 0.402951 0.403586 0.403737], observed orders [1.97, 2.06]; e [0.9723 0.9707 0.9703 0.9702]
  hough ny 40: CL(nx = 1, 2, 4, 8, 16) = [0.401524 0.403165 0.403586 0.40368  0.403701], observed orders in nx [1.96, 2.16, 2.21]
```

and writes `vlm.png`: the lattice of a tapered swept wing, cl/CL along the span for three sweeps, and the swept-wing slopes against (4).

### Verify

`python -m pytest curriculum/phase-10/lesson-11-vortex-lattice-method` runs 13 tests, listed in the metadata. Two of them guard the method's foundations: a strip in the middle of a 2000 m wing of 1 m chord gives cl = 2π sin α to 1 % (the 1/4–3/4 rule), and on a 400 m wing the aerodynamic centre is at 0.250 MAC.

### Compare

The established tool is AVL (`@drela-avl`); this lesson compares instead with an independent method and with measurements. The classic check of a coarse lattice (an untapered wing of AR 5 swept 45°, one chordwise row of 8 uniform horseshoes, the four-horseshoe-per-side example of many textbooks) gives CL_α = 3.4435 per radian with the uniform layout, and 3.1906 with Hough's inset; refined to 64 strips the two converge to 3.2020 (still falling) and 3.1707. The coarse textbook number is 8.6 % above the converged value: it is a demonstration, not an answer.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment** (with one comparison against authentic measured data). Write each prediction before running it.

1. **Lattice convergence.** The table above. Predict the orders. *Finding:* uniform and cosine spacing converge at order 0.8–0.9 (tending to 1), from above, and their e passes through values above 1 on coarse lattices (1.070 at ny = 10); Hough's inset converges at 1.97–2.06 from below, with e steady at 0.970. Chordwise refinement converges at order 2 (1.96–2.21). A lattice of nx = 4 and ny = 40 with the inset is within 4 × 10⁻⁴ of the extrapolated CL.
2. **VLM against lifting line** (unswept, untapered, nx = 2, ny = 100):

| AR | VLM CL_α (/rad) | lifting line | difference | AR × difference |
|---|---|---|---|---|
| 4 | 3.5999 | 4.0284 | −10.64 % | −42.5 % |
| 6 | 4.2048 | 4.5304 | −7.19 % | −43.1 % |
| 8.33 | 4.6262 | 4.8775 | −5.15 % | −42.9 % |
| 12 | 5.0154 | 5.1991 | −3.53 % | −42.4 % |
| 20 | 5.4286 | 5.5443 | −2.09 % | −41.7 % |
| 50 | 5.8760 | 5.9264 | −0.85 % | −42.6 % |

Predict how the difference scales. *Finding:* the lattice is always lower, by **0.43/AR** to within ±0.01: the lifting line ignores that near each tip the load falls over a distance of the order of the chord, a fraction c/b = 1/AR of the span. The two methods agree in the limit where lifting-line theory is exact, and the VLM's CL_α/2π rises towards 1 (0.935 at AR 50).
3. **Sweep** (untapered AR 6, nx = 6, ny = 40):

| sweep | VLM | Helmbold–DATCOM | difference | x_ac (MAC behind MAC LE) |
|---|---|---|---|---|
| 0° | 4.2112 | 4.5287 | −7.0 % | 0.239 |
| 15° | 4.1261 | 4.4218 | −6.7 % | 0.200 |
| 30° | 3.8474 | 4.0928 | −6.0 % | 0.176 |
| 45° | 3.3301 | 3.5174 | −5.3 % | 0.174 |
| 60° | 2.5293 | 2.6613 | −5.0 % | 0.207 |

*Finding:* the formula sits 5–7 % above the lattice at every sweep, but the trend agrees: VLM(45°)/VLM(0°) = 0.791 against 0.777 from the formula. The untapered swept wing's aerodynamic centre lies ahead of the MAC quarter chord because its MAC (at b/4) is outboard of the lift centroid. An airliner-like planform (taper 0.4, sweep 30°, AR 8.33) gives CL_α = 4.392, x_ac = 0.280 MAC, e = 0.984; 5° dihedral changes them to 4.385, 0.279 and 0.985.
4. **Measured slopes** (authentic data: NACA TN-1669, untapered wings with NACA 0012 sections, R ≈ 0.8–1.6 × 10⁶, values digitised from the report's figure with a reading uncertainty of ±0.0005 per degree; `curriculum/phase-10/data/finite_wing_cla.csv`). VLM (nx 6, ny 40, a0 = 2π implicit) in per degree, with the measured value in brackets: AR 5.16: 0.0698 (0.0689) at 0°, 0.0561 (0.0560) at 45°, 0.0430 (0.0434) at 60°; AR 2.61: 0.0508 (0.0486), 0.0448 (0.0463), 0.0367 (0.0360); AR 1.34: 0.0323 (0.0280), 0.0310 (0.0315), 0.0282 (0.0286). The Helmbold–DATCOM values at AR 5.16 are 0.0751, 0.0592, 0.0452.
5. **The trainer.** *Finding:* section above: CL_α = 4.631, e = 0.970, x_ac = 0.243 c, cruise α = 2.02°, CDi = 0.00431 (95.7 N).

Raw data format: CSV tables `spacing, ny, CL, e` and `nx, CL` (1); `AR, CLa_vlm, CLa_ll` (2); `sweep, CLa_vlm, CLa_formula, x_ac_mac` (3); `AR, sweep, CLa_vlm_per_deg, CLa_measured_per_deg` (4); plus `vlm.png`.

## 5. Validate and improve

| quantity | reference | VLM | difference | cause |
|---|---|---|---|---|
| 2-D strip slope | 2π sin α | within 1 % | < 1 % | finite span (AR 2000) |
| 2-D aerodynamic centre | 0.25 c | 0.250 MAC (AR 400) | < 0.003 | — |
| trainer CL, ny 40 → 80 | Richardson | 0.403586 → 0.403737 | 1.5 × 10⁻⁴ | order 2 |
| unswept vs lifting line | lifting line | −0.43/AR | 5.2 % at AR 8.33 | chordwise tip effect (model) |
| swept trend 45°/0° | formula 0.777 | 0.791 | +1.8 % | semi-empirical fit |
| TN-1669, AR 5.16 | measured 0.0689, 0.0560, 0.0434 /° | 0.0698, 0.0561, 0.0430 | +1.3 %, +0.2 %, −0.9 % | a0, viscosity, reading |

The lattice is verified (2-D limits, convergence at order 2) and, for this family of wings, validated: at AR 5.16 it matches the measured slopes within 1.3 %, inside the combined uncertainty of the digitised data (±0.0005 per degree, about ±0.8 %) and the unknown section slope. At AR 2.61 the differences grow to 2–5 %, and at AR 1.34 the unswept wing is overpredicted by 15 % while the swept ones agree: at very low aspect ratio the measured lift curve is nonlinear (tip vortices separate and add lift at higher α, lowering the small-α slope reading), which a linear VLM cannot represent. At AR 5.16 the Helmbold–DATCOM formula is 4–9 % high.

The engineering answer: **the trainer's wing has CL_α = 4.63 per radian (0.081 per degree), e = 0.970 and its aerodynamic centre at 24.3 % chord; at cruise CDi = 0.0043 (96 N) at 2.0°. The lifting line overestimates the lift slope by 5 % and underestimates e by 0.036: for the trainer's induced drag the difference is 4 N, for its gust response and stability derivatives 5 %.** To improve: camber (extension 1), the section slope from `P10.06`, and the fuselage.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of the 1/4–3/4 rule, (2) and (3); `vlm.png`; the five experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) on the trainer's wing with the convergence evidence.

### Adversarial self-checks

- **What if α doubles?** CL doubles (∝ sin α), CDi quadruples (∝ CL²), x_ac stays. Near the stall the real wing departs; the VLM never stalls.
- **What if AR → 0?** The linear lattice tends to the slender-wing slope (π/2)AR but misses the separated vortex lift that dominates real low-AR wings (part 4: +15 % at AR 1.34 unswept).
- **Which assumption fails first?** For the trainer, thickness and viscosity in the section slope (a few per cent); for highly swept or low-AR wings, attached flow.
- **What measurement could falsify the model?** A span-load survey (pressure taps along the span) of a swept wing showing less tip loading than the lattice predicts, or a measured neutral point far from x_ac.
- **Does the result survive a change of lattice?** With Hough's inset, yes: ny 40 → 80 changes CL by 0.04 %, and nx 4 → 16 by 0.03 %. Without it, the same refinement still changes CL by 0.7 %, and e by more.
- **Which decision changes if the uncertainty doubles?** If the lift slope were uncertain by ±10 % rather than ±5 %, a gust-load or stability margin set within that band would need a wind-tunnel or flight measurement; the induced-drag estimate (±4 N) would still be good enough for performance.

## Check yourself

1. Where does the 1/4–3/4 rule put the bound vortex and the control point, and what does it reproduce?
2. What is in row i, column j of the influence-coefficient matrix?
3. Why compute induced drag in the Trefftz plane rather than from the bound-leg forces?
4. Define the aerodynamic centre and give the formula the code uses.
5. Why does a lattice with its tip vortices at the geometric tips converge only at first order?
6. How can you measure an observed order when no exact solution exists?
7. Why is the VLM lift slope below the lifting line's, and how does the difference scale with AR?
8. Name two phenomena the VLM cannot represent.
