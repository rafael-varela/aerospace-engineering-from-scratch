# P10.02  Airfoil Geometry

> An airfoil is a thickness distribution wrapped around a camber line. Generate it from its equations, then check it the way you would check anyone else's coordinates.

| | |
|---|---|
| **Learning outcomes** | 1. Construct NACA 4-digit and standard 5-digit sections from their thickness and camber equations, with thickness laid perpendicular to the camber line, in Selig order, on cosine spacing, with an open or closed trailing edge. 2. Derive the thickness form's area 0.68508 t c², its leading-edge radius 1.1019 t², and the 5-digit constants r and k₁ from their defining conditions instead of copying them. 3. Implement and test a generator with geometric checks (maximum thickness and its position, camber, area, closure, symmetry) and measure the observed order of the area error on cosine and uniform spacing. 4. Use the generator to answer a design question about volume inside a wing, and judge what a section designation does and does not tell you. |
| **Prerequisites** | `P01.05` (Taylor series, behaviour of √x near zero), `P03.11` (truncation error and observed order), `P03.02` (root finding by bisection), `P02.09` (the Fourier-cosine idea behind the Glauert substitution, used once in §2.4), `P05.05` (area of a section). Skill check: you can integrate a polynomial term by term and use the shoelace formula for a polygon. |
| **Engineering question** | The reference trainer's wing is a NACA 2412 of chord 1.2 m. What is its cross-section area, and how much gross volume lies between a front spar at 15 % chord and a rear spar at 60 % chord over the 10 m span? Answers: 0.11849 m² (0.6857 t c²); 0.07245 m² per metre of span, 0.7245 m³ (724 L) over the span before ribs, skins and spar webs are subtracted. On the way: the section is 144.1 mm thick at 30 % chord and its camber is 24.0 mm at 40 %. |
| **Assumptions and limits** | The NACA analytical definitions (TR-460, TR-824): 4-digit thickness polynomial, 4-digit two-parabola camber, standard (non-reflexed) 5-digit cubic-plus-line camber. Real sections deviate (the 6-series, supercritical and low-Re sections of `P10.07` and `P10.14` are defined by tables or by design codes). Polygon geometry: the area and every later panel method see straight segments between points. The 5-digit constants are derived from thin-airfoil theory (`P10.04`), which is how they were defined. |
| **Analytical work** | The thickness polynomial integrated for the area; the leading-edge radius from the √x term; the closed-trailing-edge coefficient; the perpendicular construction x_u = x − y_t sin θ; the 5-digit condition r³ − 3r² + 6pr − 3p² = 0; k₁ proportional to the ideal lift coefficient; the h^{1.5} error of a polygon at a √x nose and how cosine spacing removes it. |
| **Computational deliverable** | `starter/main.py`: `cosine_spacing`, `naca4_thickness`, `naca4_camber`, `naca5_camber`, `naca5_r`, `naca5_ideal_cl`, `naca5_constants`, `naca_coordinates`, `polygon_area`, `thickness_area_coefficient`, `max_thickness`, `leading_edge_radius`, `area_between`, `area_convergence` (`parse_naca`, `surface_thickness_camber`, `observed_order`, `write_selig_dat` and `main` are given). Tests prove: spacing symmetric and clustered; maximum thickness t at 30 % chord for open and closed edges; trailing-edge thickness 0.021t or zero; leading-edge radius from the limit y_t²/2x; area coefficient equal to adaptive quadrature; camber lines continuous with zero slope at the maximum; the 5-digit constants put the maximum at p and give the design ideal c_l (and reproduce the 230 constants); Selig order closed and counter-clockwise; symmetric sections symmetric; area scaling with c²; area order 2.0 ± 0.1 on cosine and 1.5 ± 0.1 on uniform spacing. |
| **Practical activity** | **computational experiment**: the area convergence study on two spacings; the 5-digit constants derived and the quoted table checked against its own definition; vertical against perpendicular thickness; the spar-box volume of the trainer's wing. |
| **Equipment** | None. A computer with the environment from `P00.01`. (Optional: print the generated section at full scale as a rib template for `P10.23`.) |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `airfoil_geometry.png`, `naca2412.dat`, the derivations of §2.2–2.4, the experiment tables, and a short note on the trainer's spar-box volume and what it does not include. |
| **Assessment** | `assessments/phase-10/grading-rubric.md`; competencies 1–4. Pass: the area, the leading-edge radius and the 5-digit condition derived; all public tests pass; the observed orders on both spacings measured and explained; the 5-digit table checked, with the outlier identified. |
| **References** | `@ntrs-naca-tr-824-summary-of-airfoil-data` "Description of Airfoils": Method of Combining Mean Lines and Thickness Distributions; NACA Four-Digit Series Airfoils; NACA Five-Digit Series Airfoils; Supplementary Data I "Basic Thickness Forms", II "Data for Mean Lines", III "Airfoil Ordinates" (https://ntrs.nasa.gov/citations/19930090976) (required-free: the definitions, and tabulated ordinates to check your generator against). `@ntrs-naca-tr-460-78-related-airfoils` "Description of Airfoils" (https://ntrs.nasa.gov/citations/19930091108) (supplementary-free: the original 4-digit equations). `@leishman-intro-aerospace-flight-vehicles` "Airfoil Shapes" (https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/airfoil-geometries/): Airfoil Geometry; NACA Method of Drawing Airfoil Shapes; Symmetrical NACA Airfoils; Cambered NACA Airfoils (required-free, link only). `@nasa-glenn-beginners-guide-aeronautics` "Effects of Shape on Lift Interactive" (https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/foilshape/) (supplementary-free). Optional textbooks for deeper study: I. H. Abbott and A. E. von Doenhoff, *Theory of Wing Sections* (Dover; paid, inexpensive, often in libraries), the book form of TR-824's definitions and data. |
| **Extensions** | ★ Add the NACA 4-digit-modified and reflexed 5-digit (third digit 1) camber lines. ★ Fit a 4-digit section to measured coordinates (least squares on m, p, t) and report the residual. ★ Generate a section with a finite trailing-edge thickness for manufacturing (blend the closed form to a set thickness over the last 10 %) and see what it does to `P10.05`'s Kutta condition. |
| **Time** | unmeasured |

## 1. The engineering problem

Before an aerodynamicist can compute anything about a wing, someone has to say what shape it is. For most of the twentieth century that meant a NACA designation: four or five digits from which anyone can rebuild the section exactly. Every panel method, boundary-layer code and structural model in this phase starts from such a list of points. If the points are wrong (thickness added vertically instead of perpendicular to the camber line, a trailing edge left open when the solver assumes it closed, too few points around the nose) every later number is wrong in a way that is hard to trace.

Geometry also answers design questions directly. The reference trainer's wing (NACA 2412, chord 1.2 m, span 10 m) carries fuel and a spar inside the section. How much room is there between a front spar at 15 % chord and a rear spar at 60 %? That is an integral over the section, and it needs the surface to be right.

Success criteria, written first: the generated 12 % section is 12 % thick at 30 % chord (to the precision of the defining polynomial); the camber line peaks at the designated place with the designated height; the polygon area converges at second order to the exact integral; the 5-digit constants satisfy the conditions that define them; and the generator's output passes the checks you would apply to anyone else's coordinates.

## 2. Understand and derive

### 2.1 Thickness and camber

A NACA section is built from two curves. The **camber line** (mean line) y_c(x) is the shape of a zero-thickness airfoil; its height m and position p set the lift at zero angle and the pitching moment (`P10.04`). The **thickness distribution** y_t(x) is a symmetric half-thickness added on both sides. x runs along the chord from the leading edge, in chords.

The 4-digit thickness form is a polynomial with a √x term that gives a round nose:

**y_t = 5t (0.2969√x − 0.1260x − 0.3516x² + 0.2843x³ − 0.1015x⁴).**   (1)

The coefficients were chosen so that the maximum thickness is t at x ≈ 0.30 (the code finds 1.00029t at x = 0.2998), with a nose radius and trailing-edge angle that matched good sections of the 1930s. At x = 1 the coefficients sum to 0.0021, leaving an open trailing edge **0.021t** thick (0.00252c for 12 %). Replacing −0.1015 by **−0.1036** makes the sum zero and closes it; panel methods prefer the closed form because the Kutta condition (`P10.03`, `P10.05`) is applied at a single point.

The 4-digit camber line is two parabolas meeting at x = p with zero slope and height m:

y_c = (m/p²)(2px − x²) for x ≤ p,  y_c = m/(1 − p)² · ((1 − 2p) + 2px − x²) for x > p.   (2)

NACA 2412 means m = 0.02, p = 0.4, t = 0.12. The curvature jumps at x = p, a detail that matters for quadrature in `P10.04`.

### 2.2 Thickness perpendicular to the camber line

The thickness is laid off along the local normal to the camber line, not vertically. With θ = arctan(dy_c/dx):

**x_u = x − y_t sin θ,  y_u = y_c + y_t cos θ;  x_l = x + y_t sin θ,  y_l = y_c − y_t cos θ.**   (3)

![A strongly cambered, thick 4-digit section with its camber line. At one station x a rib is drawn along the local normal to the camber line: the upper point (x_u, y_u) lies y_t along the normal above (x, y_c) and the lower point (x_l, y_l) the same distance below. The normal is tilted by θ from the vertical, the angle the camber-line tangent makes with the horizontal; hollow circles mark where a vertical construction would put the points.](figures/thickness-normal-to-camber.svg)
*Figure 1. Thickness laid off perpendicular to the camber line, equation (3). The upper point moves forward and the lower point aft, by y_t sin θ each. Drawn for m = 0.10 and t = 0.30, far beyond the 2412, so that θ and the offsets are visible.*

The upper and lower points of one "rib" are at different x. Near the nose of a cambered section, where θ is several degrees, the difference between this and the vertical construction is up to 0.0018c on the 2412 and 0.0066c on a 6412 (Problem 6): small, but exactly where the suction peak sits.

**Worked example (Problem 1 in reverse).** For the 2412 at x = 0.25: y_t = 0.6(0.148450 − 0.031500 − 0.021975 + 0.004442 − 0.000396) = 0.059413; y_c = (0.02/0.16)(0.2 − 0.0625) = 0.017188; dy_c/dx = (0.04/0.16)(0.15) = 0.0375, θ = 2.148°. So (x_u, y_u) = (0.24777, 0.07656) and (x_l, y_l) = (0.25223, −0.04219).

### 2.3 Area and nose radius

Integrating (1) term by term (∫√x dx = ⅔x^{3/2}):

**A = ∫₀¹ 2y_t dx = 10t(2a₀/3 + a₁/2 + a₂/3 + a₃/4 + a₄/5) = 0.68508 t c²**   (4)

(0.68088 t c² closed). The often-quoted "0.685" is this number. Camber changes the area only at second order, through (3): the 2412's area is 0.6857 t c², 0.09 % more. For the trainer, A = 0.6857 × 0.12 × 1.2² = 0.1185 m².

Near the nose the √x term dominates, y_t ≈ 5ta₀√x, so y_t² ≈ 25t²a₀² x. A parabola y² = 2rx has radius r at its vertex, so the **leading-edge radius is r = (5ta₀)²/2 = 1.1019 t²**: 0.0159c for 12 %, 19 mm on the trainer.

### 2.4 The 5-digit camber line, derived

The standard 5-digit line is a cubic up to x = r and a straight line after:

y_c = (k₁/6)(x³ − 3rx² + r²(3 − r)x) for x ≤ r,  y_c = (k₁r³/6)(1 − x) for x > r.   (5)

Two numbers, r and k₁, are fixed by two conditions. **Position**: the maximum camber is at x = p, so the slope (k₁/6)(3x² − 6rx + r²(3 − r)) vanishes there:

**r³ − 3r² + 6pr − 3p² = 0.**   (6)

For p = 0.15 (the "30" in 23012) the root between p and 1 is r = 0.2027. **Size**: the first digit L sets the *ideal* (design) lift coefficient, c_l,i = 0.15L. In thin-airfoil theory (derived in `P10.04`; used here as a formula) c_l,i = πA₁ with A₁ = (2/π)∫₀^π (dy_c/dx) cos θ dθ, x = (1 − cos θ)/2. Every term of (5) carries k₁, so **c_l,i ∝ k₁** and k₁ = c_l,i / (c_l,i of the line with k₁ = 1). The integral is done by Gauss–Legendre quadrature split at the kink θ_r = arccos(1 − 2r). For the 230 line this gives k₁ = 15.920 (the commonly quoted 15.957 with r = 0.2025 gives c_l,i = 0.3000 too: the quoted pair is consistent). §4 checks the whole table.

### 2.5 Spacing and the polygon

A computer sees the surface as straight segments. Between two points on a smooth curve a chord misses an area of order h³, so a polygon with n uniform segments converges as h². But near the nose y ∝ √x is not smooth: the first segment from (0, 0) to (h, k√h) misses kh^{3/2}/6, so on **uniform spacing the area error falls only as h^{1.5}**. With **cosine spacing** x = (1 − cos β)/2 and uniform steps in β, x ≈ β²/4 near the nose and √x ≈ β/2 is smooth in β: second order is restored, and the nose gets many small panels. Cosine spacing also clusters points at the trailing edge, where `P10.05` applies the Kutta condition.

### 2.6 Where the model stops

- The NACA families are analytical shapes from the 1930s. Modern sections (the NACA 6-series, supercritical, the low-Re sections of `P10.14`) are tables or design-code outputs; read them with the same Selig-order conventions and check them with the same tools.
- A polygon is not the surface: every panel method inherits the h² geometric error, and an under-resolved nose spoils the suction peak long before it spoils the area.
- A manufactured section has a finite trailing edge, surface waviness and tolerances (Ladson's model was within 0.0002c, `P10.06`); "NACA 2412" on a drawing is a nominal shape.

## 3. Build the computational model

Open `starter/main.py`. `parse_naca`, `surface_thickness_camber`, `observed_order`, `write_selig_dat` and `main` are given.

### Step 1: spacing, thickness, camber

`cosine_spacing(n)` returns n + 1 stations. `naca4_thickness(x, t, closed_te)` is (1); `naca4_camber(x, m, p)` returns (y_c, dy_c/dx) from (2) with `np.where`; `naca5_camber(x, r, k1)` is (5).

### Step 2: the 5-digit constants

`naca5_r(p)` solves (6) with `scipy.optimize.brentq` on (p, 1). `naca5_ideal_cl(r, k1)` is the split Gauss–Legendre quadrature:

```python
th_r = np.arccos(1 - 2 * r)
for a, b in ((0.0, th_r), (th_r, np.pi)):
    tt = 0.5 * (b - a) * (xg + 1) + a
    total += 0.5 * (b - a) * np.sum(wg * naca5_camber(0.5 * (1 - np.cos(tt)), r, k1)[1] * np.cos(tt))
return np.pi * (2 / np.pi) * total
```

`naca5_constants(p, cl_design)` returns r and k₁ = cl_design / naca5_ideal_cl(r, 1).

### Step 3: coordinates

`naca_coordinates(code, n, closed_te, chord)` builds x, y_t, y_c, θ, applies (3), and concatenates upper (reversed: trailing edge to nose) and lower (nose excluded, to the trailing edge): 2n + 1 points in Selig order.

### Step 4: checks

`polygon_area` is the shoelace formula with the closing edge; `thickness_area_coefficient` is (4); `max_thickness` uses `minimize_scalar`; `leading_edge_radius` is §2.3; `area_between(code, x_front, x_rear, chord)` integrates the vertical thickness between two stations; `area_convergence` builds symmetric polygons on either spacing.

Running `main()` prints

```
thickness-form area coefficient: open TE 0.68508, closed TE 0.68088  (area = coeff * t * c^2)
NACA 0012: max thickness 0.120035 c at x = 0.2998 c (closed TE: 0.120014 at 0.2995); t/c target 0.12 -> ratio 1.00029
TE thickness (open) = 0.00252 c = 0.0210 t; LE radius = 0.015867 c = 1.1019 t^2
5-digit 3x0 series (p = 0.15): r = 0.2027, k1 = 15.920, max camber 1.839 % at x = 0.1500; quoted table (r, k1) = (0.2025, 15.957) gives ideal cl 0.3000
trainer NACA 2412, chord 1.2 m: area 0.11849 m^2 (0.6857 t c^2); max vertical thickness 144.1 mm at x = 0.300 c; max camber 24.00 mm at x = 0.402 c
spar box 15-60 % chord: 0.07245 m^2 per metre span; x 10.0 m span = 0.7245 m^3 = 724 L gross (before ribs, skins, spar webs)
NACA 23012: max camber 1.840 % at x = 0.143; area 0.08231 c^2
```

(the other four 5-digit lines, and the convergence table of §4, are printed too) and writes `airfoil_geometry.png` (four sections and the area-convergence plot) and `naca2412.dat` in Selig format.

### Verify

`python -m pytest curriculum/phase-10/lesson-02-airfoil-geometry` runs 19 tests: spacing endpoints, symmetry and clustering; maximum thickness t (to 0.05 %) at 30 % (±0.002) for three thicknesses and the closed edge; trailing-edge thicknesses; the nose radius as the limit of y_t²/2x at x = 10⁻⁸; the area coefficient against `scipy.integrate.quad`; the 4-digit camber line's continuity, maximum and slope (against `np.gradient`); for each 5-digit line, zero slope at p, ideal c_l = 0.3 to 10⁻⁸, continuity at r and k₁ doubling with the design c_l; the 230 constants within 0.3 %; Selig order, closure and counter-clockwise orientation; symmetry of a 0015 and recovered camber of a 4412; area scaling with c²; observed orders 2.0 ± 0.1 (cosine) and 1.5 ± 0.1 (uniform); the spar box inside the total.

### Compare

TR-824's Supplementary Data III tabulates ordinates for many of these sections; generate the same section with the open trailing edge and compare at the tabulated stations (the table is rounded to four or five figures). Your generator is the independent route for any coordinates you are given.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction first.

1. **Spacing and convergence.** Predict the observed order of the polygon area of a NACA 0012 on uniform spacing. *Finding:*

| n per surface | cosine error | uniform error |
|---|---|---|
| 10 | 1.324 × 10⁻³ | 2.450 × 10⁻³ |
| 40 | 8.319 × 10⁻⁵ | 2.995 × 10⁻⁴ |
| 160 | 5.201 × 10⁻⁶ | 3.702 × 10⁻⁵ |
| 320 | 1.300 × 10⁻⁶ | 1.304 × 10⁻⁵ |

Orders: cosine 1.99, 2.00, 2.00, 2.00, 2.00; uniform 1.52, 1.51, 1.51, 1.51, 1.50. At 320 panels per surface uniform spacing is ten times worse, entirely because of the nose.
2. **The 5-digit table, checked.** Predict whether every quoted (r, k₁) pair satisfies its own definition. *Finding:* derived (r, k₁) = (0.0581, 350.33), (0.1257, 51.577), (0.2027, 15.920), (0.2903, 6.624), (0.3913, 3.223), each with its camber peak exactly at p (1.084 %, 1.524 %, 1.839 %, 2.080 %, 2.263 % camber). The quoted pairs give ideal c_l 0.3084, 0.3019, 0.3000, 0.3001, 0.3000: the 230, 240 and 250 lines are consistent; the 210 pair is 2.8 % high, and the error is in k₁, not r (Problem 5).
3. **Thickness, vertical and perpendicular.** Predict whether a 12 % cambered section is exactly 12 % thick vertically. *Finding:* 0.12003 (0012), 0.12007 (2412), 0.12039 (6412), with the maximum moving forward to 0.292c on the 6412; the 23012's mean line peaks at 0.143, not at its camber line's 0.150.
4. **The trainer's spar box.** Predict the fraction of the section between 15 % and 60 % chord. *Finding:* 0.07245 of 0.11849 m², 61 %; 724 L over 10 m of span, gross.

Raw data format: CSV `n, err_cosine, err_uniform` (1); `p, r, k1, quoted_r, quoted_k1, quoted_ideal_cl` (2); `section, t_vertical, x_t, camber, x_camber` (3); `naca2412.dat` in Selig format.

## 5. Validate and improve

| check | theory | code | difference |
|---|---|---|---|
| area coefficient | (4): 0.68508 | quadrature of (1): 0.68508 | < 10⁻⁸ |
| max thickness of 0012 | t at 30 % | 1.00029t at 0.2998 | the polynomial's own rounding |
| nose radius | 1.1019 t² | lim y_t²/2x | < 10⁻³ relative |
| 230 constants | (6) and c_l,i = 0.3 | r = 0.2027, k₁ = 15.920; quoted pair gives 0.3000 | consistent |
| area convergence | order 2 (cosine), 1.5 (uniform) | 2.00, 1.51 | – |

The engineering answer: **the trainer's NACA 2412 section has an area of 0.11849 m²; between spars at 15 % and 60 % chord it holds 0.07245 m² per metre, 724 L gross over the 10 m span.** That is not a fuel capacity: ribs, skins, spar webs, the fuselage carry-through and the outboard aileron bays take most of it, which is why light trainers carry a small fraction of that figure in tanks. What the geometry does give exactly is the upper bound and the local depth (144 mm at 30 % chord) that the spar of `P05.10` must fit.

What to improve: compare your generator with TR-824's tabulated ordinates for the 2412 and 23012, and add the reflexed 5-digit lines. What not to improve: the number of points beyond the panel method's needs (`P10.05` shows how many it needs).

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (3), (4), the nose radius and (6); `airfoil_geometry.png` and `naca2412.dat`; the four experiment tables; and a short note (template in `projects/_template/report-template.md`, sections 1, 2, 5 and 7) on the trainer's spar-box volume and what it excludes.

### Adversarial self-checks

- **What if t → 0?** The section collapses to its camber line, the area to zero linearly in t and the nose radius quadratically; `P10.04` and `P10.05` use exactly that limit.
- **What if the camber is large?** At m = 0.09 the vertical thickness exceeds t by 0.7 % and the difference between constructions near the nose grows to several thousandths of chord; the perpendicular construction is the definition, so use it.
- **Which assumption fails first in the real world?** The trailing edge: a real wing has a finite one, and the open form's 0.021t is a manufacturing compromise, not a design choice.
- **What could falsify the generator?** Disagreement with TR-824's tabulated ordinates beyond their rounding; any section whose area does not converge at order 2 on cosine spacing.
- **Does the result survive a change of spacing or point count?** The area to 10⁻⁶ by n = 320 on cosine spacing; the spar-box area changes by less than 10⁻⁵ between n = 400 and n = 4000.
- **Which decision changes if the geometry is wrong by its tolerance?** A 0.0002c ordinate error changes the area by 0.03 %: irrelevant for volume, but not for the nose radius and the suction peak.

## Check yourself

1. What do the digits of NACA 2412 mean?
2. What is cosine spacing, and why is it used?
3. In which order are Selig-format coordinates listed, and in which direction does that go round the section?
4. What is the area of the 4-digit thickness form, and where does the number come from?
5. What is the leading-edge radius of a 4-digit section of thickness t, and why?
6. How are the 5-digit constants r and k₁ defined?
7. Why does the polygon area of a NACA section converge at order 2 on cosine spacing but only 1.5 on uniform spacing?
