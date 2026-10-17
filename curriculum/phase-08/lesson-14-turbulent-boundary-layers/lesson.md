# P08.14  Turbulent Boundary Layers

> Near the wall a turbulent layer forgets everything but the wall shear and viscosity; the log law that follows turns one velocity profile into a skin-friction measurement, and the correlations turn a Reynolds number into drag.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the law of the wall (viscous sublayer and log law) from the inner and outer scalings, state Spalding's composite profile and Coles's wake, and derive the 1/7th-power-law boundary layer from the momentum integral. 2. Implement and test the law of the wall, a Newton inversion of Spalding's formula, a Clauser fit of the friction velocity, the flat-plate skin-friction correlations (local and integrated), two mixed laminar–turbulent plate methods, and a friction law built from the log-wake profile and the momentum integral alone. 3. Check the Clauser fit against the u_τ used to build a seeded synthetic profile, the correlations against each other and against integrated local laws, and the correlations against the log-wake friction law. 4. Estimate a wing's turbulent and transitional friction drag with a defensible range, and the near-wall mesh a CFD calculation of it needs. |
| **Prerequisites** | `P08.10` (laminar boundary layer, δ*, θ, momentum integral, Blasius C_F), `P08.11` (momentum integral with pressure gradient), `P08.12` (transition point by e^N), `P08.13` (Reynolds decomposition and Reynolds stress), `P08.05` (dimensional analysis), `P03.02` (Newton's method), `P02.18` (least squares), `P03.04` (quadrature). Skill check: you can solve x + e^x = 10 by Newton's method. |
| **Engineering question** | How much friction drag does the light aeroplane's wing (3000 m, 60 m/s, chord 1.5 m, 30 m² wetted) carry with a turbulent boundary layer, compared with the all-laminar 29.7 N of `P08.10` and with transition at the e⁹ point of `P08.12`, and how thin must a CFD mesh's first cell be? Answer: turbulent from the leading edge, C_F = 0.00331–0.00342 over four correlations, 163–168 N (166 N, Prandtl–Schlichting), 5.6 times the laminar value; with transition at the e⁹ point (0.865 m) 83–95 N; with early transition (Re_x = 5 × 10⁵) 150–153 N. The turbulent layer is 25.6 mm thick at the trailing edge (laminar 3.35 mm) and the first cell for y⁺ = 1 is 8.2 μm. |
| **Assumptions and limits** | Steady, incompressible, two-dimensional mean flow over a smooth flat plate at zero pressure gradient; the law of the wall with κ = 0.41, B = 5.0 (other fits such as κ = 0.384, B = 4.17 are compared); Coles's wake with Π = 0.55; skin-friction **correlations** (empirical fits, not measurements) within their Reynolds-number ranges (the 1/7th law below about 10⁷); transition as an abrupt switch at a given Re_tr. **Profiles used here are synthetic**; published DNS profiles are named for the learner to download, not bundled. |
| **Analytical work** | Inner and outer scaling and the overlap argument for the log law; the limits of Spalding's formula; δ*/δ, θ/δ and H of the 1/7th law; δ/x = 0.3707Re_x^(−1/5) and c_f = 0.0577Re_x^(−1/5) from the momentum integral; C_F = (5/4)c_f(L); Prandtl's mixed-plate formula; the friction velocity and the viscous length of the wing. |
| **Computational deliverable** | `starter/main.py`: `log_law`, `spalding_yplus`, `spalding_uplus`, `clauser_fit`, `power_law_bl`, `cf_local`, `CF_plate`, `integrate_cf`, `mixed_plate_prandtl`, `mixed_plate_momentum`, `log_wake_state`, `log_wake_plate` (`isa_troposphere`, `coles_wake`, `laminar_CF`, `synthetic_profile` and `main` are given). Tests prove: the log law's slope is 1/κ per e-fold; Spalding's law tends to u⁺ = y⁺ and to the log law, its Newton inversion is exact to 10⁻¹⁰ and monotone; the Clauser fit recovers the synthetic u_τ to 2 % and a clean log profile exactly, and refuses a profile with no log region; the 1/7th-law plate satisfies the momentum integral to 10⁻⁴ and C_F = (5/4)c_f; the correlations agree within 8–10 %, integrated local laws match integrated ones within 8 %, and Kármán–Schoenherr satisfies its own equation; the mixed plate lies between laminar and turbulent and reduces to each; later transition means less drag; the log-wake friction law agrees with two correlations within 8 % with no correlation in it. |
| **Practical activity** | **computational experiment** with **synthetic** profiles: the spread of skin-friction correlations, a Clauser fit's sensitivity to its window and constants (and its scatter over 40 seeded profiles), the wing's drag as a function of where transition occurs by two mixed-plate methods, and the log-wake friction law against the correlations; optionally, a comparison with published DNS profiles the learner downloads. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `turbulent_bl.png`, the derivations of the log law and of the 1/7th-law plate, the experiment tables, and a half-page note giving the wing's friction drag with a range and the reasons for each end of it. |
| **Assessment** | `assessments/phase-08/grading-rubric.md`; competencies 1–4. Pass: the log law and the 1/7th-law results are derived; all public tests pass; the Clauser fit's error is reported against the known synthetic u_τ; the note gives the drag with a range built from the transition point and the correlation spread, and labels correlations and synthetic profiles as such. |
| **References** | `@leishman-intro-aerospace-flight-vehicles` "Turbulent Flows" (https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/turbulent-flows/): Wall-Bounded Turbulent Flows; and "Boundary Layer Flows" (https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/boundary-layers/): Turbulent Boundary Layers; Skin Friction Drag on an Airfoil (required-free; link only). `@libretexts-liburdy-intermediate-fluid-mechanics` Ch. 10 "Introduction to Turbulence Effects": Turbulent Boundary Layer; Scaling for Turbulent Boundary Layer Flows; Initial Laminar Starting Length; Universal Velocity Profile (required-free). `@mit-ocw-2-27-turbulent-flow-transport` outlines "Bounded Flows I: General Near-Wall Scaling Laws; Flows in Pipes, Channels, etc." (5.1) and "Introduction to Turbulent Boundary Layers" (6.4 Coles's mean-velocity profile with wake function; 6.6 smooth and rough flat plates) (supplementary-free). Published DNS profiles for the optional comparison (supplementary-free; citation terms, no license, so downloaded by the learner and never committed): `@lee-moser-2015-channel-dns-data` (https://turbulence.oden.utexas.edu/channel2015/content/README_2015.html) and `@kth-schlatter-orlu-tbl-dns-data` (https://www.mech.kth.se/~pschlatt/DATA/); `@jhtdb-johns-hopkins-turbulence-database` also hosts channel data (ODC-By). Optional textbooks for deeper study: Pope, *Turbulent Flows*, and White, *Viscous Fluid Flow* (paid; library copy; see `references/phase-08.md`). |
| **Extensions** | ★ Download a Schlatter & Örlü profile (Re_θ ≈ 4000), fit κ, B and Π to it, and rerun the log-wake friction law with your constants. ★ Add an equivalent sand-grain roughness shift ΔB(k_s⁺) to the log law and find the roughness height at which the wing's drag rises 10 %. ★ Replace the abrupt transition with a transition length and see how much the mixed-plate drag changes. |
| **Time** | unmeasured |

## 1. The engineering problem

`P08.12` found that the light aeroplane's wing boundary layer should turn turbulent a little past mid-chord in quiet air, and much sooner in rough air or with insects on the leading edge. From the transition point back, the layer is turbulent: much thicker, with a fuller profile and a much steeper velocity gradient at the wall. Turbulent friction dominates the drag of almost every aircraft, ship and pipe, and its estimate is the first line of every drag build-up (`P10.13`).

Two things make turbulent friction tractable without solving the turbulence. Near the wall, the mean velocity depends only on the wall shear, viscosity and the distance from the wall: the **law of the wall**, with its **logarithmic** region. And integrated over a plate, the friction follows **correlations** in Reynolds number that are fits to decades of measurements. This lesson builds both, uses the log law to read the wall shear off a velocity profile (the **Clauser fit**), and checks the correlations against each other and against a friction law built from the log law and the momentum integral alone.

Success criteria, written first: Spalding's formula reproduces both limits and its inversion is exact; the Clauser fit recovers the u_τ used to build a synthetic profile within 2 %; the 1/7th-law plate satisfies the momentum integral and C_F = (5/4)c_f; the correlations agree within 10 % in their range; the mixed plate lies between laminar and turbulent; the log-wake friction law agrees with the correlations within 8 %.

## 2. Understand and derive

### 2.1 Inner variables and the law of the wall

In a turbulent layer the total shear stress near the wall is nearly constant and equal to τ_w. The only velocity scale that can be built from τ_w and ρ is the **friction velocity u_τ = √(τ_w/ρ)**, and with ν the only length is the **viscous length ν/u_τ**. In these **wall units**

u⁺ = U/u_τ,  y⁺ = y u_τ/ν,

the inner layer is universal: u⁺ = f(y⁺), whatever the outer flow (Prandtl's law of the wall). Very close to the wall the turbulence is damped and the stress is viscous, τ_w = μ dU/dy, so

**u⁺ = y⁺  (viscous sublayer, y⁺ ≲ 5).**   (1)

### 2.2 The log law

Far from the wall, the outer layer does not care about ν: the velocity *defect* scales as (U_e − U)/u_τ = g(y/δ). Between the two there is an overlap where both descriptions must hold. Writing y dU/dy from each (Problem 2a), one side depends only on y⁺ and the other only on y/δ; since δ⁺ = δu_τ/ν is large and arbitrary, both must be the same constant, written 1/κ:

**u⁺ = (1/κ) ln y⁺ + B,  κ ≈ 0.41, B ≈ 5.0  (log layer, 30 ≲ y⁺, y ≲ 0.2δ).**   (2)

κ is the von Kármán constant. Different fits to data give slightly different pairs (κ = 0.384, B = 4.17 is a modern one); the two lines cross at y⁺ = 152. Between the sublayer and the log layer lies the **buffer layer**, where neither holds. **Spalding's formula** covers all three in one explicit expression for y⁺:

**y⁺ = u⁺ + e^(−κB)[e^(κu⁺) − 1 − κu⁺ − (κu⁺)²/2 − (κu⁺)³/6],**   (3)

which tends to (1) as u⁺ → 0 and to (2) as u⁺ grows. To get u⁺(y⁺), invert it by Newton's method (`P03.02`): (3) is convex in u⁺, so Newton from below converges monotonically. Above the log layer, **Coles's wake** adds (2Π/κ)sin²(πy/(2δ)) to the profile, with Π ≈ 0.55 for a flat plate at high Reynolds number.

### 2.3 The Clauser fit

If a measured profile has a log region, (2) holds there with the unknown u_τ: U = u_τ[(1/κ)ln(yu_τ/ν) + B]. Choose the u_τ that makes the points in the log region fit best, in least squares. The selection of points depends on u_τ (through y⁺), so iterate: select, fit, reselect until the selection settles. The result is the wall shear stress τ_w = ρu_τ², obtained without a balance or a sensor at the wall. Its weaknesses are its assumptions: the log law, its constants and the window.

### 2.4 The 1/7th-law plate and the correlations

A simple outer profile, u/U = (y/δ)^(1/7), has δ*/δ = 1/8, θ/δ = 7/72 and H = 9/7 = 1.29 (against 2.59 laminar: much fuller). It cannot give the wall shear (its slope is infinite at the wall), so Prandtl borrowed Blasius's pipe-flow law τ_w = 0.0225ρU²(ν/(Uδ))^(1/4). In the momentum integral dθ/dx = c_f/2 (`P08.10`), turbulent from x = 0 (Problem 2b):

**δ/x = 0.3707 Re_x^(−1/5),  c_f = 0.0577 Re_x^(−1/5),  C_F = (5/4) c_f(L).**   (4)

The constants fitted to plate data are 0.0592 (local) and 0.074 (average). Better correlations follow the log law's logarithmic dependence and hold to higher Re (all **correlations**, labelled as such):

| | local c_f | average C_F |
|---|---|---|
| power law (Re < 10⁷) | 0.0592 Re_x^(−1/5) | 0.074 Re_L^(−1/5) |
| Schultz-Grunow / Prandtl–Schlichting | 0.370 (log₁₀ Re_x)^(−2.584) | 0.455 (log₁₀ Re_L)^(−2.58) |
| White | 0.455/ln²(0.06Re_x) | — |
| Kármán–Schoenherr | — | 0.242/√C_F = log₁₀(Re_L C_F) (implicit; Newton) |
| ITTC-1957 (ship model line) | — | 0.075/(log₁₀ Re_L − 2)² |

### 2.5 The mixed plate

Transition at x_tr. **Prandtl's subtraction**: C_F = C_F,turb(Re_L) − (Re_tr/Re_L)[C_F,turb(Re_tr) − C_F,lam(Re_tr)], which assumes the turbulent layer behind x_tr is the one that would have grown from the leading edge. **Momentum match**: carry the laminar θ across transition (momentum is conserved through it) and grow a 1/7th-law layer from a virtual origin x₀ chosen so that its θ equals the laminar θ at x_tr; the drag of one side is ρU²θ(L). The two agree for early transition and drift apart as it moves back (part 4).

### 2.6 An independent friction law

Combine (3) with Coles's wake and choose δ⁺. Then Ue⁺ = u⁺(δ⁺), c_f = 2/Ue⁺², and θ and δ* follow by integrating the profile: Re_θ = Ue⁺θ⁺. Sweep δ⁺ and integrate the momentum integral dRe_x/dRe_θ = 2/c_f: that gives c_f(Re_x) from the law of the wall, the wake and momentum conservation, with **no skin-friction correlation**. It is the check the correlations need.

### 2.7 The wing

At Re_c = 4.83 × 10⁶, White's c_f at the trailing edge is 0.00288, so u_τ = 2.275 m/s, τ_w = 4.71 Pa and ν/u_τ = 8.19 μm: a first CFD cell at y⁺ = 1 is **8.2 μm** high and the log layer starts **0.25 mm** from the skin. The 1/7th-law layer is **25.6 mm** thick (δ⁺ = 3127), eight times the laminar 3.35 mm.

### 2.8 Where the model stops being valid

- **Pressure gradient**: the log law survives a mild gradient but Π and c_f change; in an adverse gradient the turbulent layer thickens and can separate (`P08.16`); flat-plate correlations ignore it.
- **Roughness**: a rough wall shifts the log law down (ΔB) and raises c_f; insects, rivets and paint matter once k_s⁺ exceeds a few (Extension 2).
- **Low Reynolds number** just after transition: the wake and the log layer are not developed (δ⁺ of a few hundred), and the correlations and (4) are least reliable there.
- **Compressibility and heat transfer** are ignored (M = 0.18).

## 3. Build the computational model

Open `starter/main.py`. `isa_troposphere`, `coles_wake`, `laminar_CF`, `synthetic_profile` (a seeded synthetic profile generator, labelled synthetic) and `main` are given.

### Step 1: the law of the wall

`log_law` is (2); `spalding_yplus` is (3). `spalding_uplus` inverts (3) by vectorised Newton steps, u ← u − g(u)/g′(u), from the guess min(y⁺, log law), until every step is below 10⁻¹².

```python
k = kappa * u
g = u + E * (np.exp(k) - 1 - k - k**2 / 2 - k**3 / 6) - yp
dg = 1 + E * kappa * (np.exp(k) - 1 - k - k**2 / 2)
u -= g / dg
```

### Step 2: the Clauser fit

`clauser_fit(y, U, nu, delta)` starts at u_τ = 0.04 max U, selects y⁺ ≥ 30 and y/δ ≤ 0.2, minimises Σ[U − u_τ(ln(yu_τ/ν)/κ + B)]² with `minimize_scalar(method="bounded")`, and repeats until the selection stops changing; it raises `ValueError` with fewer than three points.

### Step 3: power law and correlations

`power_law_bl` is (4) with the exact coefficients. `cf_local` and `CF_plate` are the table of 2.4; Kármán–Schoenherr solves g(s) = 0.242s − log₁₀Re + 2log₁₀s = 0 for s = 1/√C_F by Newton. `integrate_cf` averages a local law over a log-spaced Re_x grid from 10³ (the part below is negligible and the correlations meaningless there).

### Step 4: mixed plate and the log-wake law

`mixed_plate_prandtl` and `mixed_plate_momentum` are 2.5, in Reynolds numbers. `log_wake_state(δ⁺)` integrates the composite profile on a log-spaced y⁺ grid; `log_wake_plate(Re_x)` sweeps δ⁺ from 200 to 2 × 10⁵, integrates dRe_x = (2/c_f)dRe_θ from a 1/7th-law start, and interpolates c_f in log Re_x.

Running `main()` prints

```
wing at 3000 m, 60.0 m/s: rho 0.9091 kg/m^3, nu 1.8630e-05 m^2/s, q 1636.4 Pa, Re_c 4.831e+06, wetted area 30 m^2

integrated skin friction C_F at Re_c (correlations, turbulent from the leading edge):
  power                C_F = 0.00341, drag 167.3 N
  prandtl-schlichting  C_F = 0.00338, drag 166.1 N
  karman-schoenherr    C_F = 0.00331, drag 162.7 N
  ittc                 C_F = 0.00342, drag 167.8 N
  integral of local 'power': C_F = 0.00340; local c_f at the trailing edge 0.00273
  integral of local 'schultz-grunow': C_F = 0.00334; local c_f at the trailing edge 0.00273
  integral of local 'white': C_F = 0.00350; local c_f at the trailing edge 0.00288
  all-laminar (P08.10): C_F = 0.00060, drag 29.7 N
  mixed plate, Re_tr = 5e5 (a common textbook default): x_tr = 0.155 m; C_F Prandtl 0.00305 (149.7 N), momentum match 0.00312 (153.1 N), virtual origin 112 mm
  mixed plate, e^9 point of P08.12 (Re_tr = 2.785e6): x_tr = 0.865 m; C_F Prandtl 0.00170 (83.4 N), momentum match 0.00193 (94.8 N), virtual origin 738 mm

turbulent trailing edge (from the leading edge): delta 25.6 mm (1/7th law), theta 2.49 mm, c_f 0.00288 (White), u_tau 2.275 m/s, tau_w 4.71 Pa, viscous length nu/u_tau 8.19 um, delta+ 3127
  first CFD cell for y+ = 1: 8.2 um; log layer starts (y+ = 30) at 0.246 mm

Clauser fit to a SYNTHETIC trailing-edge profile (Spalding + Coles wake, Pi = 0.55, 1 % noise, seed 14):
  u_tau used 2.2754 m/s, recovered 2.2633 m/s (-0.53 %) from 16 points
  with kappa = 0.384, B = 4.17 instead: 2.2634 m/s (-0.53 %)

log-wake profile + momentum integral (no correlation) against correlations, local c_f:
  Re_x 1e+06: log-wake 0.00349, White 0.00376, Schultz-Grunow 0.00361, power 0.00374
  Re_x 4.83e+06: log-wake 0.00274, White 0.00288, Schultz-Grunow 0.00273, power 0.00273
  Re_x 1e+07: log-wake 0.00247, White 0.00257, Schultz-Grunow 0.00242, power 0.00236
  Re_x 1e+08: log-wake 0.00182, White 0.00187, Schultz-Grunow 0.00172, power 0.00149
  at delta+ = 1e4: Re_theta 29488, H 1.285, c_f 0.00220
```

and writes `turbulent_bl.png`: the law of the wall with the synthetic profile on it, and the local c_f correlations, the log-wake law and the laminar line against Re_x.

### Verify

Run `python -m pytest curriculum/phase-08/lesson-14-turbulent-boundary-layers`. The tests:

1. **Law of the wall**: u⁺(1) = B, slope 1/κ per e-fold; Spalding → u⁺ = y⁺ (0.2 %) and → log law (0.05) at y⁺ = 10³–10⁵, below both in the buffer layer; inversion exact to 10⁻¹⁰ and monotone from y⁺ = 0 to 10⁶.
2. **Clauser**: u_τ of a seeded synthetic profile within 2 % with only log-region points used; a clean log-law profile to 10⁻⁸; `ValueError` for a sublayer-only profile.
3. **Power-law plate**: momentum integral to 10⁻⁴; δ/x Re_x^(1/5) = 0.3707; H = 9/7, θ/δ = 7/72; C_F = 1.25 c_f exactly and `integrate_cf` within 0.3 %.
4. **Correlations**: within 8 % (integrated) and 10 % (local) in their ranges; the power law under-reads at 10⁸; integrated local laws within 8 % of Prandtl–Schlichting and Kármán–Schoenherr; Kármán–Schoenherr satisfies its equation to 10⁻¹² and falls with Re; turbulent > 3× laminar and falling more slowly.
5. **Mixed plate**: between laminar and turbulent; virtual origin between the leading edge and transition; Re_tr = Re_L gives laminar exactly and Re_tr → 0 turbulent; the momentum match with no laminar run equals 2θ/L of the 1/7th law; later transition, less drag (both methods).
6. **Log-wake law**: c_f = 2/Ue⁺², 1.2 < H < 1.45, Ue⁺ = Spalding + 2Π/κ; within 8 % of White and Schultz-Grunow from 2 × 10⁶ to 10⁸, falling with Re.
7. `main` writes the figure.

### Compare

There is no established library for these formulas; the comparison is between independent routes. The four integrated correlations span 3.3 % at the wing's Re (162.7–167.8 N). The integrals of the three local laws span 4.8 % (0.00334–0.00350). The log-wake law, built without any correlation, sits 5 % below White and within 0.4 % of Schultz-Grunow at the wing's trailing-edge Re_x. The Clauser fit recovers the synthetic u_τ to 0.5 %; it is a check of the method, not of the log law, because the profile was built from the same law.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment** with **synthetic** profiles. Write each prediction before running it.

1. **The spread of correlations.** All of them at Re = 10⁶–10⁸. Predict which fails first at high Re. *Finding:* at the wing's Re_c, 0.00331–0.00342 integrated and 0.00273–0.00288 local at the trailing edge; the power law falls 20 % below White at 10⁸ (0.00149 against 0.00187), outside its range. Pick a log-based correlation above 10⁷.
2. **The Clauser fit.** 40 seeded synthetic profiles (u_τ = 2.2754 m/s, 1 % noise, 40 points), then window and constants varied on one. Predict the bias. *Finding:* mean error −0.61 %, scatter 0.23 %, largest 1.2 %; the noise-free profile gives −0.52 %, a bias from buffer-layer and wake points inside the window. Raising the window's upper edge from y/δ = 0.1 to 0.5 moves the noise-free result from −1.0 % to +0.6 %; κ = 0.384, B = 4.17 changes it by at most 0.42 % because the two laws cross at y⁺ = 152, inside the window.
3. **Where transition is.** The wing's drag for Re_tr = 5 × 10⁵ and the e⁹ point by both mixed-plate methods. Predict which is larger. *Finding:* 149.7 and 153.1 N at Re_tr = 5 × 10⁵; 83.4 and 94.8 N at the e⁹ point (virtual origins 112 and 738 mm). For Re_L = 4.5 × 10⁶ (Problem 4) the methods differ by 1 % at Re_tr = 3 × 10⁵ and 11 % at 2 × 10⁶. Transition moves the drag by a factor of two; the method by 10 %; the correlation by 3 %.
4. **The log-wake friction law.** c_f at Re_x = 10⁶–10⁷ for Π = 0.4, 0.55, 0.7 and the alternative κ, B. *Finding:* at Re_x = 4.83 × 10⁶, 0.002832, 0.002739, 0.002647 and 0.002645, against White's 0.002876; every 0.15 in Π moves c_f about 3 %. H at δ⁺ = 10³, 10⁴, 10⁵ is 1.408, 1.285, 1.227.
5. **Optional, published data.** Download a Schlatter & Örlü zero-pressure-gradient profile (Re_θ ≈ 4000) or a Lee & Moser channel profile from their sites (citation terms; no license, so do not commit them), plot u⁺(y⁺) over Spalding's law, fit κ and B in 30 < y⁺ < 0.15δ⁺, and run your Clauser fit with them. These are DNS (simulation) data, published, and must be labelled so.

Raw data format: CSV tables `Re, kind, cf_or_CF` (1); `seed, ut_true, ut_fit, n_points` and `window, kappa, B, ut_fit` (2); `Re_tr, method, CF, D_N, x0_m` (3); `Re_x, Pi, kappa, B, cf_logwake, cf_white` (4); each synthetic file's header says "SYNTHETIC: Spalding + Coles wake, Pi, noise, seed".

## 5. Validate and improve

| quantity | analysis | code | difference | cause |
|---|---|---|---|---|
| Spalding limits | u⁺ = y⁺; log law | 0.2 %; 0.05 at y⁺ ≥ 10³ | — | higher-order terms |
| 1/7th law momentum integral | dθ/dx = c_f/2 | residual 10⁻⁴ | — | finite differences in the test |
| C_F/c_f(L) | 5/4 | 1.25 | 0 | — |
| Clauser u_τ (synthetic) | 2.2754 m/s | 2.2633 m/s | −0.53 % | buffer layer and wake in window, noise |
| wing C_F, turbulent | correlations | 0.00331–0.00342 | 3.3 % | different fits to data |
| local c_f at Re_c | White 0.00288 | log-wake 0.00274 | −5 % | Π, κ, B; the 1/7th-law start |

**The engineering answer: turbulent from the leading edge, the wing's friction drag at cruise is 166 N (163–168 N over four correlations), 5.6 times the all-laminar 29.7 N; with transition at the e⁹ point of `P08.12` it is 83–95 N (the momentum match, 95 N, is the more defensible), and with early transition at Re_x = 5 × 10⁵ it is 150–153 N. The layer is 25.6 mm thick at the trailing edge; a CFD mesh resolving the wall needs a first cell of about 8 μm.** The range is set by the transition point, not by the correlations.

What to improve: the transition point (a flight measurement) and the pressure gradient (a real aerofoil's pressure distribution, `P10.06`). What not to improve: the choice among log-based correlations.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (2) and (4); `turbulent_bl.png`; the four (or five) experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) giving the wing's friction drag with a range and the reason for each end, the first-cell height for CFD, and the labels "correlation" and "synthetic" where they apply.

### Adversarial self-checks

- **What if the speed doubles?** Re doubles, the Prandtl–Schlichting C_F falls by 11 % and q quadruples: turbulent drag rises about 3.6×, against 2.8× laminar. The e⁹ transition point moves forward to 29 % chord, so the transitional drag rises faster still.
- **What if viscosity goes to zero?** u_τ/U falls only logarithmically: turbulent friction does not vanish at high Re, unlike the laminar layer's 1/√Re.
- **Which assumption fails first in the real world?** Smoothness: at y⁺ = 1 being 8 μm, rivet heads, paint edges and insects of tens of micrometres sit well into the buffer layer and raise c_f.
- **What measurement could falsify the model?** A boundary-layer rake near the trailing edge: a Clauser fit giving u_τ 20 % above the correlation's 2.28 m/s would point to roughness or an adverse gradient.
- **Does the result survive a change of method?** Four correlations, three integrated local laws and the log-wake law agree within 5 % at the wing's Re.
- **Which decision changes if the uncertainty doubles?** The correlation spread (±3 %) is irrelevant to design; the transition point's uncertainty (83–153 N) decides whether laminar-flow surfaces are worth their cost.

## Check yourself

1. Define u_τ, y⁺ and u⁺, and say why they are the natural scales near the wall.
2. Name the regions of the inner layer and the y⁺ range of each.
3. What does a Clauser fit do, and what does it assume?
4. Why is turbulent skin friction so much larger than laminar, and why does it fall so slowly with Reynolds number?
5. State the power-law local and average skin-friction correlations and the relation between them.
6. Compare the wing's turbulent and laminar friction drag at cruise.
7. What does Prandtl's mixed-plate formula assume about the turbulent layer, and when does it err?
8. Which conservation laws does this lesson use, in which form, and how is the turbulent stress modelled?
