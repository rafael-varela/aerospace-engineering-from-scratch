# P07.14  Heat Transfer: Radiation

> In space there is nothing to convect to. Every watt a spacecraft makes leaves as light, and how much leaves depends on two numbers for the coating and one for the geometry: how strongly it emits, how much sunlight it absorbs, and how much of the sky it can see.

| | |
|---|---|
| **Learning outcomes** | 1. State Planck's law, integrate it to σT⁴, find the Wien peak, and compute the fraction of black-body energy in a wavelength band; use it to explain why a white paint can absorb little sunlight and still emit well in the infrared. 2. Define the view factor as a double area integral, prove reciprocity A_i F_ij = A_j F_ji and the summation rule Σ_j F_ij = 1, and use closed forms (parallel and perpendicular rectangles, coaxial disks, two-dimensional strips, crossed strings). 3. Implement a numerical double-area integration of the view-factor kernel and show that it converges to the closed forms, at second order for parallel plates and first order when the surfaces share an edge. 4. Set up and solve the radiosity equations for an enclosure of N gray, diffuse surfaces with mixed temperature and heat-flow conditions, and check them against the black-body limit, the parallel-plate formula and the second law. 5. Size a spacecraft radiator in eclipse and in sunlight, and decide how much a neighbouring solar array costs. |
| **Prerequisites** | `P07.13` (convection, the linearised radiation coefficient), `P07.12` (thermal resistance networks), `P01.13` (linear systems), `P03.04` (numerical integration), `P01.05` (integration), `P03.02` (root finding). Skill check: you can solve Ax = b with `numpy.linalg.solve` and write a double sum over two grids with NumPy broadcasting. |
| **Engineering question** | A spacecraft radiator painted with zinc-oxide white paint (infrared emittance 0.93, solar absorptance 0.12–0.18, Lienhard Table 10.5) must reject 200 W from electronics at 20 °C. How big must it be in eclipse, how big if sunlight can fall on it, and how much hotter does it run beside a deployed solar array? The answers: 0.514 m² with no sun; 1.079 m² in normal sunlight at α = 0.15 (0.884–1.384 m² across the paint's range); a 1.0 m wide array face at 330 K along one edge (view factor 0.269) heats a 0.6 m² panel from 8.8 °C to 35.1 °C, and holding 20 °C beside it takes 0.759 m². |
| **Assumptions and limits** | Opaque, gray, diffuse surfaces (emittance and absorptance independent of wavelength and direction within the infrared; sunlight handled with a separate solar absorptance); each surface isothermal with uniform radiosity; non-participating medium (vacuum or dry air over short paths); no obstruction between surfaces unless built into the enclosure; space as a black sink at about 3 K. The radiator has an insulated back and no conduction losses; the array temperature is given, not solved. Real coatings degrade (solar absorptance rises with ultraviolet exposure), real surfaces are partly specular, and real orbits add Earth infrared and albedo loads that this lesson's radiator model takes only as an input. |
| **Analytical work** | Planck's law, its integral and its peak; the band fraction as a function of λT; the view-factor integral, reciprocity, summation, and the closed forms; Hottel's crossed-string rule; radiosity J = εE_b + (1 − ε)G and the two equations per surface; the three-resistance two-surface formula and its limits; the reradiating surface; the linearised coefficient h_rad and its error; the radiator energy balance in eclipse and in sunlight. |
| **Computational deliverable** | `starter/main.py`: `blackbody_emissive_power`, `planck_spectral`, `blackbody_fraction`, `vf_parallel_rectangles`, `vf_perpendicular_rectangles`, `vf_coaxial_disks`, `vf_parallel_strips`, `vf_perpendicular_strips`, `vf_crossed_strings`, `vf_numerical`, `solve_enclosure`, `two_surface_exchange`, `h_rad_linear`, `h_rad_mean`, `radiator_temperature`, `radiator_area` (`box_view_factors`, `radiator_with_array`, the constants and `main` are given). Tests (19) prove: Planck integrates to σT⁴ to 10⁻⁶ at three temperatures and peaks at λT = 2898 μm K; band fractions in [0, 1] and functions of λT only; every closed form reciprocal and correct in its limits; the six faces of a box summing to one to 10⁻¹²; numerical integration converging to the closed forms at order 2 (parallel) and 1 (shared edge), with discrete reciprocity to round-off; crossed strings equal to the strip formulas and the triangle rule; enclosure heat flows summing to zero, the black limit, the parallel-plate and concentric-sphere formulas, zero exchange in an isothermal enclosure, the reradiating-wall network formula; the linearisation error; and the radiator balances. |
| **Practical activity** | **computational experiment**: a convergence study of the numerical view factor; a radiator sizing study across coatings and sun angles; the radiator beside a solar array of growing width; a three-surface enclosure with reradiating walls; a multilayer-insulation (radiation-shield) stack; and the accuracy of linearised radiation over the temperature range of the home lab `P07.15`. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `radiation.png`, the derivations of reciprocity, the radiosity equations and the radiator area, the experiment tables, and a half-page note answering the engineering question with the coating and view-factor assumptions that matter most. |
| **Assessment** | `assessments/phase-07/grading-rubric.md`; competencies 1–4. Pass: reciprocity, summation and the radiosity equations derived; all public tests pass; the observed convergence orders of the numerical view factor reported with the reason for the difference; the radiator area in eclipse and in sun reproduced by hand; the note states which assumption most changes the area. |
| **References** | `@lienhard-heat-transfer-textbook` (required-free; https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf): §1.3 Modes of heat transfer, subsection "Radiation"; §2.3, subsection "Resistance for thermal radiation"; Ch. 10 "Radiative heat transfer": 10.1 The problem of radiative exchange (Black body; Table 10.1 Total emittances for a variety of surfaces); 10.2 Kirchhoff's law (The gray body approximation); 10.3 Radiant heat exchange between two finite black bodies (View factor reciprocity; Calculation of the black-body view factor, F1–2; Table 10.2; Table 10.3); 10.4 Heat transfer among gray bodies (The electrical circuit analogy when more than two gray bodies are involved in heat exchange; Algebraic solution of multisurface enclosure problems); 10.6 Solar energy (Selective emitters, absorbers, and transmitters; Table 10.5 Solar absorptance and infrared emittance for several surfaces near 300 K). `@mit-16-unified-thermodynamics-propulsion-notes` Ch. 19 "Radiation Heat Transfer": 19.1 Ideal Radiators; 19.2 Kirchhoff's Law and "Real Bodies"; 19.3 Radiation Heat Transfer Between Planar Surfaces; 19.4 Radiation Heat Transfer Between Arbitrary Surfaces (required-free; https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node136.html; link only). `@mit-ocw-16-050-thermal-energy` Part 3: 9.0 Radiation Heat Transfer; 9.3 Radiation Heat Transfer Between Planar Surfaces; 9.4 Radiation Heat Transfer Between Black Surfaces of Arbitrary Geometry (supplementary-free). `@openstax-university-physics-vol2` 1.6 Mechanisms of Heat Transfer (Radiation) (supplementary-free). |
| **Extensions** | ★ Replace the midpoint rule in `vf_numerical` with Gauss–Legendre points on each surface and a subtraction of the edge singularity; measure the new orders. ★ Add Earth infrared and albedo to the radiator for a 500 km orbit with your own view factor to the Earth (a sphere seen from a plate), and find the hot-case area. ★ Build the enclosure of a radiator between two body-mounted panels (a U-shaped bay) from `vf_perpendicular_rectangles` and `vf_parallel_rectangles`, and find the best width-to-depth ratio. |
| **Time** | unmeasured |

## 1. The engineering problem

A satellite in orbit is a box of electronics in a vacuum. Conduction carries heat from the boards to the structure, but from the structure outward there is only radiation. The usual solution is a **radiator**: a panel on a face that never sees the Sun, coated to emit strongly in the infrared, connected to the electronics by heat pipes or a conductive plate. Its size is a mass and packaging decision made early in a design, and it is set by three numbers: the coating's infrared emittance ε, its solar absorptance α_s if sunlight can reach it, and the fraction of the panel's view that is cold space rather than warm spacecraft.

The last of these is a **view factor**: the fraction of the radiation leaving one surface that arrives at another. It depends only on geometry, it is a double integral over two surfaces, and it is the part of radiation calculations most often done wrong. A solar array deployed beside a radiator takes part of its view; the array is warm (it converts only a fraction of the sunlight it absorbs into electricity), so the radiator both loses part of its sink and gains a heat source.

This lesson builds the tools in order: black-body emission and its spectrum (why a coating can be "white" to the Sun and "black" in the infrared); view factors by closed form and by numerical integration, checked against each other and against the two rules every set of view factors must obey; the radiosity method for an enclosure of gray surfaces; and the radiator. Success criteria, written first: Planck's law integrates to σT⁴ to 10⁻⁶; the six view factors of a box sum to one to 10⁻¹²; numerical integration converges to the closed forms at the order its error analysis predicts; every enclosure's net heat flows sum to zero; the black and parallel-plate limits are reproduced exactly.

## 2. Understand and derive

### 2.1 Black bodies and the spectrum

A **black body** absorbs all radiation that falls on it; it is also the best possible emitter at every wavelength (Lienhard 10.1). Its spectral emissive power, the power per unit area per unit wavelength, is Planck's law:

**e_b,λ(λ, T) = 2πhc² / (λ⁵ [exp(hc/(λkT)) − 1])**,

with Planck's constant h, the speed of light c and Boltzmann's constant k (exact SI values in the code). Two consequences follow, and the code checks both numerically. Integrating over all wavelengths gives the **Stefan–Boltzmann law** e_b = σT⁴, σ = 5.670374419 × 10⁻⁸ W/(m² K⁴): at 300 K, 459.30 W/m², and `quad` of Planck's law returns the same to six figures. The peak sits at **λ_max T = 2898 μm K** (Wien's law): 9.66 μm at 300 K, 0.50 μm for the Sun's 5772 K effective temperature (Lienhard 10.6).

With x = hc/(λkT), the fraction of e_b emitted below λ is a function of λT only. The code computes it as (15/π⁴)∫_x^∞ t³/(eᵗ − 1) dt. For the Sun, 97.9 % of the energy lies below 3 μm; for a 300 K surface, 8.7 × 10⁻⁵. Sunlight and a radiator's own emission occupy almost disjoint parts of the spectrum.

### 2.2 Real surfaces, Kirchhoff and the gray body

A real surface emits εσT⁴ with an **emittance** ε < 1 that depends on the material and its finish: polished aluminium 0.04–0.06, heavily oxidised aluminium 0.20–0.33, black gloss paint 0.90 (Lienhard Table 10.1). **Kirchhoff's law** (10.2) says that at each wavelength and direction a surface absorbs as well as it emits: α_λ = ε_λ. If ε_λ is the same at all wavelengths that matter (a **gray** surface), then α = ε for any incoming radiation.

Because sunlight and infrared emission are at such different wavelengths, a surface can be gray within the infrared and still have a solar absorptance α_s quite different from its infrared emittance. Lienhard's Table 10.5 lists pairs: zinc-oxide white paint α_s = 0.12–0.18 with ε_IR = 0.93; acrylic white paint 0.26 and 0.90; carbon black in acrylic 0.94 and 0.83; pure aluminium 0.09 and 0.1. The radiator model uses α_s for sunlight and ε for everything in the infrared.

### 2.3 The view factor

Two surfaces A₁ and A₂ exchange radiation along every line joining a point on one to a point on the other. A diffuse surface emits into each direction in proportion to the cosine of the angle from its normal, and a patch dA₂ at distance s, tilted by θ₂ from that line, subtends a solid angle cos θ₂ dA₂/s². Adding up (Lienhard 10.3, "Calculation of the black-body view factor"),

**F₁₂ = (1/A₁) ∫_A₁ ∫_A₂ cos θ₁ cos θ₂ / (π s²) dA₂ dA₁**,

the fraction of the radiation leaving A₁ that strikes A₂. The integrand is symmetric in 1 and 2, so

**A₁F₁₂ = A₂F₂₁** (reciprocity),

and since everything leaving surface i goes somewhere, for an enclosure of N surfaces

**Σ_j F_ij = 1** (summation),

where F_ii = 0 for a flat or convex surface. These two rules are the first check on any set of view factors.

**Closed forms.** The integral can be done for simple shapes (Lienhard Tables 10.2 and 10.3): directly opposed equal rectangles a × b a distance c apart; two rectangles sharing an edge of length l at right angles, widths w and h; coaxial parallel disks; and, for infinitely long two-dimensional surfaces, opposed and perpendicular strips. For any two-dimensional pair, **Hottel's crossed-string rule** gives F₁₂ = [(sum of the crossed strings) − (sum of the uncrossed strings)]/(2L₁), where the strings are stretched between the ends of the two surfaces. For opposed strips of width 1 and spacing 0.5 it gives 0.618034, the same as √(1 + (h/w)²) − h/w; the tests check every 2-D closed form this way.

**The summation rule as a test.** The inside of a rectangular box has six faces: each sees the opposite face (parallel rectangles) and four neighbours (perpendicular rectangles). Built from the two closed forms, the six rows must each sum to one. For a 1.0 × 0.8 × 0.5 m box they do, to 10⁻¹²: an independent check of both formulas, because neither was derived with the other in mind. For a unit cube, opposite faces have F = 0.19982 and adjacent faces 0.20004.

### 2.4 The radiosity method

In an enclosure of gray, diffuse surfaces, each surface i has an emissive power E_b,i = σT_i⁴, receives an **irradiation** G_i (W/m²) from the others, and sends out a **radiosity** J_i: what it emits plus what it reflects,

J_i = ε_i E_b,i + (1 − ε_i) G_i,  A_i G_i = Σ_j A_j F_ji J_j = A_i Σ_j F_ij J_j (by reciprocity).

The net heat leaving surface i is what goes out minus what comes in:

**Q_i = A_i (J_i − G_i) = A_i (J_i − Σ_j F_ij J_j)**.

Eliminating G gives two forms, depending on what is known (Lienhard 10.4, "Algebraic solution of multisurface enclosure problems"):

- known temperature: **J_i − (1 − ε_i) Σ_j F_ij J_j = ε_i σ T_i⁴**;
- known net heat flow: **J_i − Σ_j F_ij J_j = Q_i / A_i**.

That is N linear equations for N radiosities (`P01.13`). Afterwards each Q follows from the definition, and for a known-Q surface the temperature from σT⁴ = J + Q(1 − ε)/(εA). A surface with Q = 0 is **reradiating** (an insulated wall): it absorbs and re-emits everything, and its emittance does not matter.

Three checks follow from the equations. Summing Q_i over all surfaces and using reciprocity and summation gives **Σ Q_i = 0**: energy is conserved. If every ε = 1, J = σT⁴ and Q_i = Σ_j A_i F_ij σ(T_i⁴ − T_j⁴), the black-body exchange. And if all temperatures are equal, every Q is zero for any emittances: the second law forbids net exchange in an isothermal enclosure.

For **two surfaces** the equations collapse into three resistances in series (Lienhard's circuit analogy):

**Q₁₂ = σ(T₁⁴ − T₂⁴) / [(1 − ε₁)/(ε₁A₁) + 1/(A₁F₁₂) + (1 − ε₂)/(ε₂A₂)]**.

For large parallel plates (F₁₂ = 1, A₁ = A₂) the flux is σ(T₁⁴ − T₂⁴)/(1/ε₁ + 1/ε₂ − 1): 3594.5 W/m² between plates at 600 K (ε 0.8) and 300 K (ε 0.6).

**A worked enclosure.** Two 1 m² plates 1 m apart, one at 600 K (ε = 0.8) and one at 300 K (ε = 0.6), with the four sides of the cube between them insulated (reradiating, ε = 0.5). The radiosity solution gives Q = +2666.65 W from the hot plate, −2666.65 W into the cold one, zero from the walls (which settle at 529.6 K), and ΣQ = 0. The walls more than double the exchange: with only the direct view (F₁₂ = 0.1998) the plates would exchange 1163.6 W. In circuit terms the reradiating walls add a parallel path 1/(A₁F₁R) + 1/(A₂F₂R) beside the direct resistance 1/(A₁F₁₂), and the tests check the solver against that network formula to 10⁻¹⁰.

### 2.5 Linearised radiation

Writing q = εσ(T⁴ − T_sur⁴) = h_rad(T − T_sur) gives **h_rad = εσ(T² + T_sur²)(T + T_sur)** exactly (Lienhard §2.3, "Resistance for thermal radiation"). The common approximation 4εσT_m³ with T_m the mean temperature is low by ((T − T_sur)/(T + T_sur))² relative to (1 + that): 0.22 % for a body at 50 °C in a 21 °C room, 0.82 % at 80 °C. Using 4εσT_sur³ instead (the surroundings' temperature) is 13.6 % low at 50 °C: always use the mean or the exact form. Lesson `P07.15` relies on this.

### 2.6 The spacecraft radiator

A flat radiator of area A with an insulated back, carrying Q from the electronics, in equilibrium:

**ε σ A (T⁴ − T_sink⁴) = Q + α_s G_s cos θ_sun A + ε q_IR A**,

where G_s = 1361 W/m² is the solar irradiation at 1 AU (Lienhard 10.6), θ_sun the angle of the Sun from the panel's normal, and q_IR any infrared arriving from a planet (absorbed with α = ε by Kirchhoff). Space is a sink at about 3 K, whose σT⁴ (5 × 10⁻⁶ W/m²) is negligible. Solving for the area at a required temperature:

**A = Q / [ε σ (T⁴ − T_sink⁴) − α_s G_s cos θ_sun − ε q_IR]**.

At 20 °C, ε = 0.93, a ZnO-painted panel emits 389.5 W/m². In eclipse, **A = 200/389.5 = 0.514 m²**. In normal sunlight with α_s = 0.15 it absorbs 204.2 W/m², leaving 185.3 W/m² net: **A = 1.079 m²**, more than twice. With the Sun 60° from the normal, 0.696 m². Across the paint's α_s range (0.12–0.18) the sunlit area spans 0.884–1.384 m²; with acrylic white (α_s = 0.26) it is 8.685 m²; carbon black (α_s = 0.94) cannot reach 20 °C in normal sunlight at all, since it absorbs more than it emits. That is why radiators face away from the Sun and why the coating's solar absorptance, which rises as the paint ages in ultraviolet, is a sizing parameter.

**A neighbour that is warm.** Put a solar array face, 1.0 m wide, at 330 K (typical, illustrative for an operating array) at right angles to a 1.0 × 0.6 m radiator along its 1.0 m edge. The radiator sees the array with F = 0.269 (perpendicular rectangles) and space with the rest. The three-surface enclosure (radiator with Q = 200 W given, array at 330 K, space black at 3 K) gives a radiator temperature of **35.1 °C instead of 8.8 °C**. To hold 20 °C beside the array the panel must grow to 0.759 m², 48 % more than in open space.

### 2.7 Where the model stops

- **Non-gray and specular surfaces.** Coatings vary with wavelength inside the infrared; metal and glass surfaces reflect specularly, so the diffuse view factors misplace their reflections.
- **Non-uniform surfaces.** A radiator fed by heat pipes is warmer near them; one radiosity per surface averages that out. Split surfaces into smaller ones when temperatures vary.
- **Participating media.** Water vapour and CO₂ absorb and emit in bands (Lienhard 10.5); over metres of humid air or in combustion gases, the enclosure is not transparent.
- **Orbit environment.** Earth infrared, albedo and the changing Sun angle make the radiator's load time-dependent; the hot and cold cases bracket it.

## 3. Build the computational model

Open `starter/main.py`. The constants (σ, h, c, k, G_s, the surface data from Tables 10.1 and 10.5), `box_view_factors`, `radiator_with_array` and `main` are given; they call your functions.

### Step 1: black body

`blackbody_emissive_power` is σT⁴. `planck_spectral` is 2.1; use `np.expm1` and suppress the harmless overflow far below the peak with `np.errstate(over="ignore")`. `blackbody_fraction` integrates t³/(eᵗ − 1) from x = hc/(λkT) upward with `quad` and multiplies by 15/π⁴.

### Step 2: closed-form view factors

Each formula of 2.3 is a few lines. For `vf_perpendicular_rectangles(l, w, h)` take logs of each factor separately (the W² and H² powers overflow otherwise for large ratios). `vf_crossed_strings(a1, b1, a2, b2)` takes four points; label them so the uncrossed strings join a1–a2 and b1–b2.

### Step 3: the numerical view factor

`vf_numerical(o1, u1, v1, o2, u2, v2, n)` describes each planar surface as o + s u + t v, 0 ≤ s, t ≤ 1. Divide each into n × n cells and evaluate the kernel between cell centres:

```python
R = P2[None, :, :] - P1[:, None, :]            # every pair of cell centres
s2 = np.einsum("ijk,ijk->ij", R, R)
cos1 = np.clip(R @ n1, 0, None) / np.sqrt(s2)
cos2 = np.clip(-(R @ n2), 0, None) / np.sqrt(s2)
F12 = np.sum(cos1 * cos2 / (np.pi * s2)) * dA1 * dA2 / A1
```

Orient each normal toward the other surface's centroid. The double sum is symmetric in the two surfaces, so A₁F₁₂ = A₂F₂₁ holds for the discrete sum to round-off, before it has converged.

### Step 4: the enclosure

`solve_enclosure(A, F, eps, T, Q)` builds the N × N matrix row by row from the two forms of 2.4, solves with `np.linalg.solve`, and returns J, Q and T. `two_surface_exchange` is the three-resistance formula; `h_rad_linear` and `h_rad_mean` are 2.5.

### Step 5: the radiator

`radiator_temperature` and `radiator_area` are 2.6; `radiator_area` raises `ValueError` when the absorbed flux exceeds the emitted flux.

Running `main()` prints (abridged)

```
black body 300 K: sigma T^4 = 4.5930e+02 W/m^2, integral of Planck (10 nm - 1 mm) = 4.5930e+02
fraction below 3 um: sun (5772 K) 0.979, 300 K body 8.70e-05
parallel rectangles 1.0 x 0.8 m, gap 0.5 m: closed form F = 0.376013
  numerical n = 40: 0.376128, error +1.16e-04
perpendicular rectangles, shared edge 1.0 m, widths 0.8 and 1.2 m: closed form F12 = 0.240853, F21 = 0.160569
  numerical n = 40: 0.244166, error +3.31e-03
box 1.0 x 0.8 x 0.5 m: row sums [1. 1. 1. 1. 1. 1.]; max reciprocity residual 0.0e+00 m^2
radiator 200 W at 20 C, ZnO white paint eps 0.93, alpha_s 0.15:
  area, no sun 0.514 m^2; sun at normal incidence 1.079 m^2; sun 60 deg from normal 0.696 m^2
```

and writes `radiation.png`: the Sun's and a 300 K body's spectra, the convergence of the numerical view factor, and the radiator temperature against area.

### Verify

`python -m pytest curriculum/phase-07/lesson-14-radiation` runs 19 tests, listed in the metadata row. Two deserve a comment. The numerical view factor converges at **order 2** for the parallel rectangles (errors 2.93 × 10⁻³, 7.26 × 10⁻⁴, 1.81 × 10⁻⁴ at n = 8, 16, 32: the midpoint rule on a smooth integrand) but only at **order 1** for the perpendicular pair (1.62 × 10⁻², 8.21 × 10⁻³, 4.14 × 10⁻³): along the shared edge s → 0 and the kernel is singular, so the cells next to the edge carry an error that halves, not quarters, with each refinement. That is the same lesson as the elliptic lift in `P05.12`: test a numerical method on a case that is not smooth. The second is the isothermal enclosure: whatever the emittances, every net heat flow is zero to 10⁻⁹ of σT⁴, which is the second law as a unit test.

### Compare

Four independent routes meet. Closed forms and numerical integration agree at the predicted orders. The summation rule ties the parallel and perpendicular closed forms together. The crossed-string rule reproduces the 2-D formulas exactly. The radiosity solver reproduces the black-body sum, the parallel-plate and concentric-sphere formulas and the reradiating-wall network formula to 10⁻¹⁰–10⁻¹². What none of these checks can test is the input: the emittance and absorptance of the real coating, and the temperature of the real neighbour.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction before running it.

1. **Numerical view factors.** n = 5, 10, 20, 40 for the two rectangle pairs of part 3. *Predict* the two orders. *Finding:* errors 7.67 × 10⁻³, 1.87 × 10⁻³, 4.64 × 10⁻⁴, 1.16 × 10⁻⁴ (order 2) and 2.55 × 10⁻², 1.30 × 10⁻², 6.59 × 10⁻³, 3.31 × 10⁻³ (order 1). A 1 % view factor for a shared-edge pair needs n ≈ 40; a better quadrature is Extension 1.
2. **Coatings and the Sun.** The 200 W, 20 °C radiator in normal sunlight with each coating of Table 10.5. *Predict* which cannot work. *Finding:* ZnO white 0.884 (α 0.12) to 1.384 m² (α 0.18); acrylic white 8.685 m²; carbon black impossible; no sun 0.514 m². A 0.6 m² panel runs at 8.8 °C in eclipse and 44.6 °C in normal sun: if the Sun can reach the radiator, either the panel or the electronics' limits must change.
3. **The neighbour.** A 330 K array face of width 0.5, 1.0 and 2.0 m along the 1.0 m edge of the 0.6 m² radiator. *Predict* whether the effect saturates. *Finding:* F = 0.214, 0.269, 0.294; radiator 30.6, 35.1, 37.0 °C. It saturates because a perpendicular plane can occupy at most half the hemisphere; most of the effect comes from the first half-metre.
4. **Cold case.** With the electronics off in eclipse, the radiator falls toward the 3 K sink; holding the 0.6 m² panel at −20 °C takes a 130 W heater. *Predict* first whether a smaller radiator would save heater power, and at what cost in the hot case.
5. **Radiation shields (multilayer insulation).** A wall at 293.15 K (ε 0.9) facing black space through n thin shields of polished aluminium (ε = 0.05 on both faces, Table 10.1's 0.04–0.06), treated as large parallel plates in series: each shield adds a resistance 2/ε − 1. *Finding:* 376.889, 10.440, 2.135, 1.071 and 0.536 W/m² for n = 0, 1, 5, 10, 20. One low-emittance shield cuts the loss by a factor of 36; after that each layer helps less. Real blankets do worse than this ideal (layers touch and conduct; seams and penetrations leak), so this is a lower bound.
6. **Linearisation over the home-lab range.** Exact h_rad against 4εσT_m³ and 4εσT_sur³ for a body at 50 °C and 80 °C in a 21 °C room, and at 600 K against 300 K. *Finding:* mean-temperature form low by 0.220 % and 0.824 %, and by 10.0 % at 600/300 K; the surroundings form low by 13.6 %, 25.6 % and 73.3 %.

Raw data format: CSV tables `n, F_num_parallel, err_parallel, F_num_perp, err_perp` (1); `coating, alpha_s, eps, area_m2` (2); `w_array_m, F_ra, T_rad_C` (3); `n_shields, q_Wm2` (5); `T_K, T_sur_K, h_exact, err_mean_pct, err_sur_pct` (6); plus `radiation.png`.

## 5. Validate and improve

| quantity | theory (§2) | code | difference | cause |
|---|---|---|---|---|
| ∫ Planck, 300 K | σT⁴ = 459.30 W/m² | 459.30 | < 10⁻⁶ | quadrature |
| box row sums | 1 | 1 | < 10⁻¹² | two closed forms, independent |
| parallel rectangles, n = 40 | 0.376013 | 0.376128 | 1.2 × 10⁻⁴ | midpoint rule, order 2 |
| perpendicular, n = 40 | 0.240853 | 0.244166 | 3.3 × 10⁻³ | edge singularity, order 1 |
| reradiating enclosure Σ Q | 0 | 0.0 W | round-off | conservation |
| parallel plates 0.8/0.6, 600/300 K | 3594.5 W/m² | 3594.5 | < 10⁻¹² | identity |
| radiator area, eclipse | Q/(εσT⁴) = 0.514 m² | 0.514 | — | — |

The mathematics is verified to round-off or to its predicted convergence rate. The engineering answer is limited by inputs: **the radiator needs 0.514 m² if it never sees the Sun and 1.08 m² if it can face the Sun squarely, a factor set by the coating's α_s/ε; a warm array beside it costs another 48 % of area (0.759 m²).** The largest uncertainty is the solar absorptance at end of life, which for white paints rises as ultraviolet darkens them; the next is the array's temperature, which this lesson takes as given.

What to improve: solve the array's temperature from its own energy balance (it absorbs sunlight on its front and exchanges with the radiator on its back) instead of fixing it at 330 K, and split the radiator into strips to see its temperature gradient toward the array. What not to improve: the numerical view factor beyond 1 %, which is far smaller than the uncertainty in α_s.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of reciprocity, summation, the radiosity equations and the radiator area; `radiation.png`; the experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) giving the radiator area for each case, the assumption that most changes it, and the design rule you would write for the spacecraft layout.

### Adversarial self-checks

- **What if the electronics' power doubles?** In eclipse the area doubles (A ∝ Q at fixed T). In sunlight it more than doubles only if the absorbed flux is part of Q; here it scales the same way, because the solar load is per square metre: 2.158 m² at normal incidence.
- **What if the radiator temperature limit rises to 40 °C?** Emission grows as T⁴: (313.15/293.15)⁴ = 1.30, so the eclipse area falls by 23 %. High-temperature electronics shrink radiators fast.
- **Which assumption fails first in the real world?** The coating: α_s rises over the mission, and a radiator sized on the beginning-of-life value runs hot at the end. Then the gray, diffuse surfaces: a specular, metallised neighbour reflects differently from the view-factor picture.
- **What measurement could falsify the model?** A thermal-vacuum test with a calibrated heater on the panel and a cold shroud: the panel's equilibrium temperature at known Q gives εσ(T⁴ − T_shroud⁴) directly. A panel hotter than predicted with the array deployed would show the array's view factor or temperature is larger than modelled.
- **Does the result survive a change of method?** View factors by closed form, numerical integration and crossed strings agree; the enclosure solver agrees with the circuit formulas; the radiator balance is algebra. The result changes only with the inputs.
- **Which decision changes if the uncertainty doubles?** If α_s were uncertain over 0.10–0.25 instead of 0.12–0.18, the sunlit area would span 0.79–4.06 m², and the design rule "radiators never see the Sun" would become non-negotiable rather than an optimisation.

## Check yourself

1. Write the view-factor integral and use it to prove reciprocity in one line.
2. A closed box has six faces. What must each row of its view-factor matrix sum to, and why?
3. Where does the Planck spectrum of a 300 K body peak, and why does that let a white paint emit well while absorbing little sunlight?
4. Write the two radiosity equations, one for a surface of known temperature and one for a surface of known net heat flow.
5. Why does a reradiating (insulated) wall's emittance not affect the exchange between the other surfaces?
6. Two large parallel plates have emittances 0.8 and 0.6. What factor multiplies σ(T₁⁴ − T₂⁴) to give the net flux?
7. Why did the numerical view factor converge at second order for parallel plates but only first order for perpendicular plates sharing an edge?
8. A radiator of ε = 0.93 must reject 200 W at 20 °C. Estimate its area in eclipse by hand.
