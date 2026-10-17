# P05.18  Home Lab: Cantilever Beam Deflection 🔬

> Beam theory says a 25 × 3 mm aluminium strip, loaded 290 mm from a clamp, deflects 2.1 mm for every newton. Write that number down with its uncertainty, then hang the weights and find out whether it is right.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the deflected shape of a cantilever under a point load by integrating E I w″ = M, including the straight part beyond the load and the sag under the strip's own weight, and say why differencing two readings removes the sag. 2. Predict the deflection at the load point with an uncertainty budget from E, the width, the thickness and the load position, by first-order propagation and by a seeded Monte Carlo, and identify the dominant term (thickness, which enters cubed). 3. Choose a load range that keeps the clamp stress below half the yield strength and the deflection inside the range where small-deflection theory holds to better than 1 %. 4. Measure load against deflection, fit the compliance with a weighted line, recover E with its uncertainty, and compare it with the prediction by a normalised error E_n. 5. Fit the deflected shape along the strip to the cubic plus a rigid clamp rotation, and decide from the numbers whether your clamp's compliance is detectable at your reading resolution. |
| **Prerequisites** | `P05.12` (E I w″ = M, cantilever boundary conditions, the finite-difference beam), `P05.10` (the flexure formula σ = Mc/I), `P05.05` (I = bh³/12), `P05.06` (Young's modulus and yield strength), `P05.14` (margins against yield), `P00.13` (safety), `P00.14` (measurement, calibration and uncertainty), `P00.15` (the first experiment, normalised error), `P03.10` (nonlinear boundary-value problems with `solve_bvp`). Skill check: you can propagate the uncertainty of y = c x³/z through the relative-uncertainty rule. |
| **Engineering question** | A 25 × 3 mm aluminium strip is clamped to a table with 300 mm free and loaded at a = 290 mm with up to 800 g. Does its deflection follow δ = P a³/(3EI) with a typical E of 69 ± 2 GPa, within the uncertainty of a ruler read from a phone photograph? Can you tell a compliant clamp from a low modulus? The prediction, before any load: 2.095 mm/N ± 3.7 % (16.4 ± 0.6 mm at 800 g), clamp stress 61 MPa at 800 g against a 72.5 MPa limit. |
| **Assumptions and limits** | Euler–Bernoulli beam (plane sections, shear deformation negligible for a strip 100 times longer than it is thick); linear elastic, isotropic material; small deflections (linear theory over-predicts by 1 % when δ/a reaches about 0.10); a perfectly rigid clamp (the assumption most likely to fail at home); a point load that stays vertical; a strip of uniform section. E is a typical handbook value, not a property of your piece: aluminium alloys sit within a few per cent of 69 GPa, and the alloy of a hardware-store strip is usually not marked. |
| **Analytical work** | Integration of E I w″ = P(a − x) with w(0) = w′(0) = 0; the straight tail beyond the load; self-weight deflection; the clamp stress 6Pa/(bh²) and the safe load; the uncertainty of δ and of E by the relative-uncertainty rule; the clamp rotation term θ₀x and why a load–deflection line cannot see it; the elastica and the limit of linear theory. |
| **Computational deliverable** | `starter/main.py`: `rect_second_moment`, `deflection_point_load`, `deflection_self_weight`, `clamp_stress`, `max_safe_load`, `predict_deflection`, `predict_deflection_mc`, `fit_line`, `modulus_from_slope`, `fit_profile`, `en_number`, `fd_cantilever`, `elastica_deflection` (`compliance_shape`, the SYNTHETIC data generator and writer, `read_csv` and `main` are given). Tests (15) prove: I and its uncertainty; the cubic's boundary conditions, continuity at the load and straight tail; curvature equals M/EI; second-order convergence of the finite differences to the cubic; self-weight against its closed form; the safe load puts the clamp stress at the allowable; GUM propagation against Monte Carlo; exact recovery of a line and χ² statistics on seeded noise; `modulus_from_slope` inverts the forward model; `fit_profile` recovers a planted clamp rotation; E_n; the elastica tends to linear theory at small load and falls below it at large load; the SYNTHETIC files are labelled and read back. |
| **Practical activity** | **physical experiment**: clamp a metal strip to a table, write a dated prediction with an uncertainty, hang known masses up to 800 g (never more than 2 kg in any case), read the deflection at the load and along the strip, and reduce the readings to E and to a clamp-rotation estimate. See `projects/p05-cantilever/`. A learner without a kit completes the lesson with the SYNTHETIC dataset in `projects/p05-cantilever/data/synthetic/`, labelled as synthetic. |
| **Equipment** | Minimum: an aluminium or steel strip or a flat metal ruler, a C-clamp, a ruler, a phone camera, coins or water in a bag, a kitchen scale (from `p00-calibration`), string. Recommended: a caliper. Enhanced: a dial gauge on a magnetic or clamp stand. See `projects/p05-cantilever/bom.yaml` (reuses `p00-calibration`, `p00-first-experiment` and `p04-pendulum` items). |
| **Safety** | See `projects/p05-cantilever/safety.md`. Summary box at the start of part 4. |
| **Evidence** | `data/prediction.md` (dated before any load), `data/calibration.csv`, `data/raw/*.csv`, the load–deflection plot with the prediction band, the profile plot with the cubic, E with its uncertainty budget and E_n, the clamp-rotation estimate with its uncertainty, `cantilever_lab.png`, the report. |
| **Assessment** | `assessments/phase-05/grading-rubric.md`; all four competencies. Threshold: a prediction with uncertainty recorded before loading; all public tests pass; E reported with an uncertainty budget and E_n against the prediction; the clamp rotation either detected or shown to be below the resolution, with the numbers; the safety box completed. |
| **References** | `@mit-ocw-3-11-mechanics-of-materials` Roylance, module "Beam Displacements": Multiple integration; Superposition; module "Stresses in Beams": Normal Stresses; appendix "Material Properties" (required-free): the cantilever closed form, the check that the strip stays elastic, and handbook moduli to compare with. `@libretexts-baker-haynes-engineering-statics` 10.2 Moments of Inertia of Common Shapes (required-free): I = bh³/12. `@mit-ocw-16-001-unified-materials-structures` Lab 02 "Unified Materials and Structures Lab 2", §1 "Deflections of a cantilevered beam subject to concentrated loads" (required-free): the laboratory version of this experiment, used as the model; the home procedure here is different. |
| **Extensions** | ★ Repeat with a steel rule of 0.8 mm thickness and loads up to the safe limit, and measure the departure from linear theory that the elastica predicts (Problem 3). ★ Clamp the strip on edge (bending about its strong axis) and predict how much stiffer it becomes; check the clamp first, because the stress for the same deflection is much higher. ★ Pluck the loaded strip and film its vibration: the first natural frequency of a cantilever with a tip mass M is close to (1/2π)√(3EI/((M + 0.24 m_strip)a³)); compare E from the frequency with E from the static slope (`P04.06`). |
| **Time** | unmeasured |

## 1. The engineering problem

Every structure on an aircraft is checked twice: once for strength (does it break?) and once for stiffness (does it deflect too much?). A wing that is strong enough can still be too flexible, because its twist under load changes the lift it makes, and because a control surface hinged to it can bind. Before an aircraft flies, its wing is loaded on the ground, step by step, with whiffle trees or sandbags, and the measured deflections are compared with the analysis. If they disagree, the model of the structure is wrong somewhere, and the strength analysis built on the same model is suspect too.

You will run the same test at the scale of a kitchen table. A metal strip clamped to a table edge is a cantilever, like a wing clamped at its root. Hang a mass near the free end and beam theory (`P05.12`) predicts the deflection everywhere along the strip from three numbers: the load, the distance from the clamp, and the bending stiffness EI. You know I from a caliper (or a ruler, less well), and E from a handbook, with a spread. That makes a prediction with an uncertainty you can write down before you hang anything.

Then the experiment does what ground tests do: it checks the model and finds what it left out. Two things are most often left out at home. The clamp is not rigid: a C-clamp on a wooden table lets the strip rotate a little at the root, which makes the strip look softer than it is. And the thickness is not exactly what the label says: E I depends on h³, so a 1 % error in thickness is a 3 % error in the prediction. This lesson shows how to tell those effects apart with the numbers you can actually measure, and how to decide, before starting, which of them your equipment can see at all.

Success criteria, written first. Before loading: a dated prediction table with uncertainties, and a computed safe load. After: the compliance and E with their uncertainties, E_n < 2 against the prediction or a named and tested cause, a χ² that matches the reading uncertainty, and a clamp-rotation estimate with its uncertainty. In code: the finite-difference beam agrees with the cubic to better than 1 µm, and the elastica agrees with linear theory at small load.

## 2. Understand and derive

### 2.1 The deflected shape

Put x along the strip from the clamp (x = 0) to the free end (x = L), and measure the deflection w downward, so a hanging load makes it positive. A downward load P at x = a gives, at a section x < a, a bending moment M(x) = P(a − x), counted positive when it curves the strip downward; beyond the load the moment is zero. The moment–curvature relation of `P05.12` with this sign convention is

**E I w″ = P(a − x)** for x ≤ a, and E I w″ = 0 for x > a.

The clamp holds w(0) = 0 and w′(0) = 0. Integrate twice:

E I w′ = P(ax − x²/2),  E I w = P(ax²/2 − x³/6) = P x²(3a − x)/6.

So

**w(x) = P x²(3a − x)/(6EI)** for 0 ≤ x ≤ a.

At the load, **δ = w(a) = P a³/(3EI)**, and the slope there is P a²/(2EI). Beyond the load the strip carries no moment, so it stays straight and continues at that slope:

**w(x) = P a²(3x − a)/(6EI)** for a ≤ x ≤ L.

The deflection is linear in P (double the load, double the deflection) and depends on a cubed: moving the load 1 % further out raises δ by 3 %. Write it as w(x) = P g(x, a)/(EI), where g(x, a) = x²(3a − x)/6 or a²(3x − a)/6 is the shape function, a length cubed. `compliance_shape` returns it.

### 2.2 The strip

For a rectangular strip of width b and thickness h bent about the axis along its width (flat side up, the easy direction), the second moment of area is (`P05.05`; `@libretexts-baker-haynes-engineering-statics` 10.2)

**I = b h³/12.**

For the lesson's strip, b = 25.0 mm and h = 3.00 mm: I = 56.25 mm⁴. With E = 69 GPa (a typical handbook and supplier data-sheet value for 6061 and 6063 aluminium, the alloys hardware-store strip is usually made of; Roylance's appendix "Material Properties" tabulates handbook moduli to compare with; the ± 2 GPa covers the spread between alloys and sources), EI = 3.881 N·m².

The weak-axis bending is the point: the same strip bent on edge would have I = h b³/12 = 3906 mm⁴, 69 times stiffer. Aircraft structure is shaped to put material far from the bending axis, which is why spars are deep and skins are stiffened.

### 2.3 Self-weight and why you difference

The strip sags under its own weight before you hang anything: a uniform load q = ρ g b h (1.99 N/m for this strip, with ρ = 2700 kg/m³) gives

w_q(x) = q x²(6L² − 4Lx + x²)/(24EI),

0.52 mm at the tip of the 300 mm strip. The equation is linear, so loads superpose: the deflection with a hanging mass is w_q + w_P. If you read the ruler once with no load and once with the load, and subtract, w_q cancels, and so does the offset of the ruler's zero. **Always record two readings, unloaded and loaded, for every load**, rather than one zero reading at the start: the difference also cancels any slow drift (the clamp settling, the table moving, the camera shifting), as long as the two readings are taken a few seconds apart.

### 2.4 Is the strip safe and elastic?

The largest bending stress is at the clamp edge, where the moment Pa is largest (`P05.10`; Roylance "Stresses in Beams", Normal Stresses):

**σ = M c/I = P a (h/2)/(b h³/12) = 6 P a/(b h²).**

At 800 g that is 60.7 MPa. The weakest common extruded alloy, 6063-T5, has a typical yield strength of 145 MPa (handbook and supplier data-sheet value; a 6061-T6 strip is about twice as strong). The lab's allowable is **half the yield strength**, 72.5 MPa, which gives

P_max = F_allow b h²/(6a) = **9.38 N, 956 g.**

Why half? The alloy of your strip is not marked, the clamp edge concentrates stress, and permanent set would spoil the experiment (the strip would no longer return to zero) long before anything broke. The lesson's loads stop at 800 g. The strain at the clamp is then 0.09 %, well inside the elastic range.

For any other strip, compute both: the safe load from this formula with the yield strength of its material divided by two, and the 2 kg hard limit of the safety box. Use the smaller.

### 2.5 The prediction and its uncertainty

Reading at the load point, δ = P a³/(3EI) = 4 P a³/(E b h³). Every input is a product or power, so the relative uncertainties add in quadrature with their exponents as weights (`P00.14`):

**(u_δ/δ)² = (u_E/E)² + (u_b/b)² + (3u_h/h)² + (3u_a/a)² + (u_P/P)²**.

For the lesson's strip: E 69 ± 2 GPa contributes 2.9 %; h = 3.00 ± 0.02 mm (a caliper, allowing for variation along the strip) contributes 3 × 0.67 % = 2.0 %; a = 290 ± 1 mm contributes 1.0 %; b = 25.00 ± 0.05 mm contributes 0.2 %; P from a kitchen scale to 1 g contributes less than 0.2 % above 100 g. Total **3.7 %**. The handbook E dominates, and the thickness is next, even though it is measured to 0.02 mm, because it enters cubed. With only a ruler for the thickness (± 0.25 mm on 3 mm) the thickness term would be 25 %, and the experiment would measure your ruler, not the strip: a caliper is in the recommended tier for this reason.

`predict_deflection` does this with numerical partial derivatives (it handles a reading point that is not the load point, where the formula is not a pure power law), and `predict_deflection_mc` checks it with 200 000 normal samples:

| mass [g] | δ [mm] | u (GUM) [mm] | Monte Carlo mean, sd [mm] |
|---|---|---|---|
| 100 | 2.05 | 0.08 | 2.06, 0.08 |
| 400 | 8.22 | 0.30 | 8.23, 0.30 |
| 800 | 16.44 | 0.60 | 16.46, 0.61 |

The Monte Carlo mean is 0.1 % high because δ is convex in h and E (it goes as 1/h³ and 1/E), a second-order effect that first-order propagation ignores; at 3.7 % uncertainty it is negligible. The predicted **compliance** (deflection per newton) is 2.095 mm/N ± 3.7 %.

### 2.6 From the measurement back to E

Load the strip in steps and fit a straight line δ = s P + c. The slope s is the measured compliance; inverting the formula,

**E = a³/(3 s I)**, with (u_E/E)² = (u_s/s)² + (3u_a/a)² + (u_b/b)² + (3u_h/h)².

The handbook term has gone: E is now the output. The intercept c should be zero within its uncertainty; a significant intercept means a zero error that is not removed by differencing (a hanger that was not counted, a strip that was already bent). The **χ²** of the fit, Σ(r_i/u_i)² over the residuals r_i, should be close to the degrees of freedom (number of points minus two) if your reading uncertainty u_i is honest: much larger means the scatter is bigger than you claimed, much smaller means you overestimated it.

Compare the measured compliance with the prediction, or the measured E with the handbook value, with the normalised error of `P00.15`:

**E_n = |x₁ − x₂| / √(U₁² + U₂²)**, with expanded uncertainties U = 2u.

E_n < 1 is consistent, E_n > 2 is a disagreement that needs a cause.

### 2.7 A clamp that rotates

A C-clamp on a wooden table is not a built-in end. Under the moment Pa the strip and the wood under the clamp yield slightly, and the strip rotates at the clamp by θ₀ ≈ Pa/k_θ, where k_θ is the clamp's rotational stiffness in N·m/rad. The strip then deflects by an extra rigid-body amount θ₀x:

**w(x) = P g(x, a)/(EI) + θ₀ x.**

At the load point this adds Pa²/k_θ. That is also proportional to P, so **the load–deflection line cannot see the clamp**: a compliant clamp just makes the slope steeper, and the slope turns it into a low E. To separate the two, use the *shape*. The beam term grows as x²(3a − x), the clamp term as x; with readings at several x along the strip under one load, a least-squares fit of w to the two columns P g(x, a) and x gives 1/EI and θ₀ separately (`fit_profile`). A rigid clamp gives θ₀ consistent with zero.

Whether your equipment can see it is a calculation, not a hope. With k_θ = 2000 N·m/rad (a plausible C-clamp on a wooden table), at 700 g θ₀ = 1.0 mrad and the extra deflection at the load point is 0.29 mm, 2 % of the total. The fit's uncertainty on θ₀ with seven ruler readings along the strip, each difference uncertain by √2 × 0.30 mm, is **4.7 mrad**: the clamp effect is five times smaller than the resolution, invisible. With a dial gauge (0.01 mm per reading, so 0.014 mm on each difference) the uncertainty drops to **0.16 mrad**, and a 1 mrad rotation stands out at six standard uncertainties. That calculation is part of your prediction.

### 2.8 How far does linear theory go?

E I w″ = M assumes small slopes. The exact (elastica) equation, in terms of the angle θ(s) of the strip along its arc length s, with the load staying vertical, is

E I θ″ = −P cos θ,  θ(0) = 0,  θ′(a) = 0,  w(a) = ∫₀ᵃ sin θ ds.

As the strip bends, the load's lever arm a cos θ shrinks, so the true deflection is *smaller* than the linear one. `elastica_deflection` solves it with `solve_bvp` (`P03.10`):

| mass [g] | linear [mm] | elastica [mm] | linear high by | δ/a |
|---|---|---|---|---|
| 200 | 4.11 | 4.11 | 0.02 % | 0.014 |
| 500 | 10.27 | 10.26 | 0.13 % | 0.035 |
| 800 | 16.44 | 16.38 | 0.33 % | 0.056 |
| 2000 (computed only) | 41.10 | 40.28 | 2.03 % | 0.139 |

The error depends only on the dimensionless load P a²/(EI), so it is the same function of δ/a for every strip: linear theory is **1 % high when δ/a ≈ 0.10**. The 2000 g row is computed only; that load would stress this strip to 152 MPa, past its yield. For this strip, the nonlinearity is ten times smaller than the reading uncertainty. For a thin steel rule loaded to its safe limit it is not (Problem 3).

### 2.9 Where the model fails

In order of how likely they are at home: the clamp rotates (2.7); the thickness is not what you think, or varies along the strip (2.5); the load hangs from a loop that slides, so a is not where you measured it; the reading is taken with parallax (the eye or camera not level with the strip's edge); the load swings; the strip had a permanent bend before you started (harmless if you difference, as long as it does not change). Beam theory itself, shear deformation and large deflection, comes last.

## 3. Build the computational model

Open `starter/main.py`. The constants (g, the typical E and yield strength, the density), `compliance_shape`, the SYNTHETIC generator `make_synthetic_cantilever_data`, `write_synthetic_csv`, `read_csv` and `main` are given. The generator and `compliance_shape` call `deflection_point_load`, so they run once you have written it.

### Step 1: the section and the closed forms

`rect_second_moment(b, h, u_b, u_h)` returns I and its standard uncertainty (2.2, 2.5). `deflection_point_load(x, P, a, EI)` is the two-piece cubic of 2.1, vectorised over x with `numpy.where`. `deflection_self_weight` is 2.3. `clamp_stress` and `max_safe_load` are 2.4.

### Step 2: the prediction

`predict_deflection(P, a, E, u_E, b, u_b, h, u_h, u_a, x=None, u_P=0)`: write the deflection as a function of the vector (P, a, E, b, h), take central differences with a relative step of 10⁻⁶, and sum (c_i u_i)². When `x` is None the reading is at the load point and moves with a; otherwise it is at a fixed mark. `predict_deflection_mc` draws normal samples with `numpy.random.default_rng(seed)` and returns the mean and standard deviation.

### Step 3: fitting

`fit_line(P, d, u_d)` is weighted least squares for d = sP + c with the covariance (AᵀWA)⁻¹ giving the uncertainties, χ² and the degrees of freedom. `modulus_from_slope` is 2.6; for a fixed reading mark it uses g(x, a) instead of a³/3. `fit_profile(x, w, P, a, u_w)` solves the two-column least-squares problem of 2.7 and returns EI, θ₀, their uncertainties and χ². `en_number` is 2.6.

### Step 4: two independent checks of the theory

`fd_cantilever(P, a, L, EI, n)` solves E I w″ = M on n + 1 nodes with the three-point second difference (`P05.12`); the zero slope at the clamp enters through a ghost node w₋₁ = w₁. `elastica_deflection` is 2.8. The elastica is a two-point boundary-value problem; give `solve_bvp` the linear solution as its first guess.

### Verify

`pytest curriculum/phase-05/lesson-18-home-lab-cantilever-deflection` runs 15 property tests: boundary conditions, curvature, convergence, limits and inverse problems, as listed in the metadata. Then run `python main.py`. The output starts:

```
strip 25.0 x 3.00 mm: I = 56.25 +/- 1.13 mm^4, E = 69 +/- 2 GPa (typical), EI = 3.8813 N m^2
safe load (clamp stress <= F_ty/2 = 72.5 MPa): 9.38 N = 956 g; stress at 800 g: 60.7 MPa
self weight 1.99 N/m: tip sag 0.52 mm before loading (removed by differencing)
```

### Compare

At 700 g the finite-difference solution with n = 200 agrees with the cubic to 0.197 µm, against a deflection of 14 mm. The second difference is exact for a cubic everywhere except in the interval containing the load, where M has a kink; that one interval sets the error, which falls as 1/n² (the test checks the order). This is the computational check of the closed form; the experiment is the physical one.

## 4. Experiment at home, or investigate real data

**Mode: physical experiment.** The bill of materials, the authoritative safety box, the full procedure, the report template and the data formats are in `projects/p05-cantilever/` (`bom.yaml`, `safety.md`, `procedure.md`, `report-template.md`, `data/README.md`). Read the safety box first, every session. If you have no kit, run the whole reduction on the SYNTHETIC files in `projects/p05-cantilever/data/synthetic/` (or on `make_synthetic_cantilever_data` with your own seed) and label every table and plot **SYNTHETIC**; you still complete the lesson, but you do not get to claim a measurement.

> **Safety (summary; the authoritative box is `projects/p05-cantilever/safety.md`).** Hazards: the hanging mass falling onto feet or the floor; the clamp slipping and the strip springing up; a strip snapping and throwing fragments; pinched fingers in the clamp; the strip's sharp cut ends. Energy limits: **total hanging mass ≤ 2 kg** under any circumstances, and never more than the safe load computed for your strip (956 g for the lesson's strip; the lesson stops at 800 g). At 800 g the strip stores about 0.06 J; a 2 kg mass falling 0.75 m from a table to the floor carries 15 J. Controls: **nothing hung over feet**; a box or cushion under the hanger; **clamp checked (tightened, strip not moved) before each load**; masses added gently, never dropped; cut ends deburred or taped; **safety glasses whenever loading a strip that can snap** (any plastic ruler, a hardened steel rule, a strip near its safe load). Supervision: none for an adult; children only with an adult loading. Location: a sturdy table that does not rock, clear floor below. Stop criteria: **visible permanent set** (the strip does not return to its unloaded reading within 0.5 mm), **cracking**, whitening or a creak; the clamp moves; the table rocks; the hanger loop slides. Not permitted: **loading glass, carbon-fibre rods or anything that can splinter** (including wooden rulers), hanging mass over 2 kg, hanging anything over feet, standing on or leaning on the loaded strip.

The setup, seen from the side:

```
   C-clamp (checked before each load)
      ║        x →                                       ruler held vertical
   ═══╬══════════════════════════════════════════╗      beside the strip, read
   ███████  strip, flat side up                   ║ ← end    from a phone photo
   ███████ ─────────────────────────────┬─────────┘          at the strip's height
   ███████   table                      │ a = 290 mm from the clamp edge
   ███████                              │ string loop in a filed notch
   ███████                             [bag]  coins / water, weighed
   ███████                              ⇣
   ███████                         box or cushion below; no feet underneath
```

**Prediction first (in `data/prediction.md`, dated).** Measure b at three points, h at five points along the strip with a caliper (or the best you have), the clamp edge position and a. Compute I ± u_I, the safe load for your strip, and the predicted δ ± u at 100, 200, …, up to your top load at the load point, with the uncertainty budget term by term (which term dominates?). Predict the compliance and the deflection at the mark positions you will read along the strip at your top-but-one load. Then compute, with `fit_profile` on noise-free predicted data, the uncertainty of θ₀ your reading method will give, and state whether a 1 mrad clamp rotation will be detectable. Date the file and do not edit it afterwards.

**Procedure in outline** (full version in `procedure.md`). Clamp the strip with its flat side up and a 10–30 mm overlap on the table, edge parallel to the table edge. File a small notch at a so the string loop cannot slide. Fix the ruler vertically beside the strip's end (taped to a heavy book or a second clamp, not touching the strip). Photograph the reading level with the strip's top edge; read to 0.5 mm, or better from the zoomed photo against the graduations. For each load: photograph unloaded, hang the bag gently, wait until the swinging stops, photograph loaded, remove the load, and confirm the strip returns to its unloaded reading. Go up in 100 g steps, then repeat the series in a different order. At the top-but-one load, read the deflection at five or more marks along the strip (move the ruler, keep the load). Weigh the bag, string and hook together on the kitchen scale for every load.

**Raw data format** (`data/raw/YYYY-MM-DD_<tier>_loads.csv`, one row per load, never edited; `reading_after_mm` is the unloaded reading after the load comes off, the permanent-set check):

```
run,mass_g,x_mm,reading_unloaded_mm,reading_loaded_mm,reading_after_mm,photo_unloaded,photo_loaded,note
1,101,290,151.5,153.5,151.5,IMG_0412,IMG_0413,
1,203,290,151.6,155.8,151.5,IMG_0414,IMG_0415,
```

## 5. Validate and improve

Run the reduction on your data (or the SYNTHETIC files) and put three comparisons side by side.

**Prediction against measurement.** The SYNTHETIC dataset (seed 518, `projects/p05-cantilever/data/synthetic/`) shows what a consistent result looks like. It was generated with a "truth" that differs from the nominal inputs on purpose: E = 67.5 GPa, h = 3.01 mm, b = 25.02 mm, a clamp with k_θ = 2000 N·m/rad, and reading noise of 0.30 mm on every reading.

| quantity | prediction (before loading) | SYNTHETIC result (seed 518) | E_n (k = 2) |
|---|---|---|---|
| compliance at x = a [mm/N] | 2.095 ± 0.078 | 2.090 ± 0.056 | 0.03 |
| intercept [mm] | 0 | +0.13 ± 0.26 | within 1 u |
| χ² / dof | about 1 | 5.6 / 7 | consistent |
| E from the slope [GPa] | 69 ± 2 (typical) | 69.2 ± 2.4 | 0.03 |
| E from the profile, clamp removed [GPa] | 69 ± 2 | 68.4 ± 7.1 | 0.04 |
| clamp rotation θ₀ at 700 g [mrad] | 0 (rigid clamp assumed) | 1.18 ± 4.70 | not detectable |

Everything agrees, and it would be wrong to conclude that the clamp is rigid or that E is 69.2 GPa. The generator's true compliance, with its clamp, is 2.160 mm/N, of which 0.042 mm/N (2 %) is the clamp; read through the nominal b and h, it corresponds to an apparent E of 66.9 GPa. This seed's noise happened to land the slope 1.3 standard uncertainties below the truth (another seed, 15, gives 66.1 ± 2.3 GPa). The fit is consistent with the truth and with the prediction because its uncertainty is honest, not because it found the truth exactly. And the 2 % clamp effect is buried: it is smaller than the slope's 2.7 % uncertainty, and the 1.0 mrad rotation is five times smaller than the profile's 4.7 mrad resolution for θ₀. That is the result 2.7 predicted, and with a ruler it is the correct report: "clamp rotation not detectable; upper bound about 2 × 4.7 = 9 mrad at 700 g".

| pattern in your data | likely cause | test |
|---|---|---|
| E low by 3–10 %, E_n > 2, χ² fine | clamp rotating, or h thinner than measured | profile fit with a dial gauge; re-measure h at the clamp edge |
| E high by several per cent | a measured too long (loop slid inward), h measured on a burr | re-measure a under load; deburr and re-measure h |
| significant positive intercept | uncounted hanger mass; first load seating the clamp | weigh the hanger; repeat the series in reverse order |
| χ² ≫ dof | reading uncertainty underestimated: parallax, a swinging load, the ruler touching the strip | photograph level with the edge; wait for the swing to stop |
| unloaded reading drifts after big loads | permanent set (stop) or the clamp creeping | stop criteria; re-tighten and repeat a low load |
| deflection grows over a minute at constant load | creep of a plastic strip or of the table under the clamp | use metal; clamp onto a hard surface |

**Code against closed form.** The finite differences agree with the cubic to 0.2 µm and the elastica tends to the linear result at small load (part 3, the tests). Hand-check one point: δ at 400 g from your readings, and E from that single point with a calculator.

**Improve and re-run.** The two biggest terms in the E budget from the SYNTHETIC fit are the slope (2.7 %, the ruler's 0.3 mm noise) and the thickness (2.0 %). A dial gauge cuts the slope term below 0.1 %; then thickness dominates, and measuring h at more points with a better caliper (± 0.005 mm gives 0.5 %) is the next step. Re-run one series after the improvement and report both results.

## 6. Defend your work

Submit: `data/prediction.md` (dated before loading), `data/calibration.csv`, `data/raw/*.csv`, your analysis script, the load–deflection plot with k = 2 error bars and the prediction band, the profile plot with the cubic and the fitted clamp rotation, the E and θ₀ table with E_n, `cantilever_lab.png`, and a report following `projects/p05-cantilever/report-template.md`. If you used the SYNTHETIC dataset, every table and figure says so in its caption.

### Adversarial self-checks

- **What if the thickness is 2.9 mm, not 3.0?** δ rises by (3.0/2.9)³ − 1 = 10.7 %. A strip sold as "3 mm" can easily be that far off; this is why h is measured, at several points, before the prediction.
- **What if the load goes to zero?** The deflection goes to zero and so does every systematic effect proportional to load; what remains is the reading noise and any offset. A non-zero intercept is therefore not a beam effect.
- **Which assumption fails first?** The rigid clamp. Then the load position (a loop that slides). Large deflection and shear come last for a metal strip at these loads; for a thin steel rule near its safe load, large deflection comes earlier (Problem 3).
- **What measurement could falsify the model?** A profile whose shape is not x²(3a − x) plus a straight line in x: for example, a deflection beyond the load that is not straight, which would mean the strip is not uniform or is touching something; or a load–deflection line that curves upward, which would mean yielding.
- **Does the result survive a change in method?** E from the load–deflection slope at one point, E from the profile shape, and (extension) E from the vibration frequency use different combinations of the inputs; the static methods share the clamp error, the frequency method sees it differently.
- **Which decision changes if the uncertainty doubles?** At twice the reading noise the clamp-rotation bound doubles to about 19 mrad, and E from the slope carries about 6 %, so a 3–5 % softening from the clamp, which a better kit would flag, is no longer distinguishable from a good result. The identification of the metal does not change: steel is three times stiffer than aluminium, so even a 10 % measurement of E says which one you have. No uncertainty this kit reaches can tell 6061 from 6063, whose typical moduli are the same.

## Check yourself

1. Write the deflection at the load point of a cantilever with a point load at distance a. By what factor does it change if the thickness of a rectangular strip doubles?
2. Why is the strip's thickness the most important dimension to measure? With h = 3.00 ± 0.02 mm, what does it contribute to the relative uncertainty of δ?
3. What is the bending stress at the clamp of a 25 × 3 mm strip loaded with 400 g at 290 mm?
4. Why record an unloaded and a loaded reading for every load, instead of one zero reading at the start?
5. A clamp rotates by θ₀ = Pa/k_θ. What does this add to the deflection at the load point, and why can a load–deflection line not tell it apart from a lower E?
6. At roughly what ratio δ/a does small-deflection theory over-predict the deflection by 1 %, and in which direction is the error?
7. Your measured E is 69.2 ± 2.4 GPa (k = 1) and the handbook says 69 ± 2 GPa. What is E_n, and what may you conclude?
8. Why must the prediction be written and dated before any load is hung?
