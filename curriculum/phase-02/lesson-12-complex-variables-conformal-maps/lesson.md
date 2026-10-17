# P02.12  Complex Variables & Conformal Maps

> An analytic function is a flow. Solve the flow around a circle, bend the plane until the circle is a wing, and the flow bends with it.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the Cauchy–Riemann equations from the requirement that a complex derivative exist, show that the real and imaginary parts of an analytic function are harmonic, and explain why the complex potential of an incompressible, irrotational flow is analytic. 2. Implement and test, from scratch, numerical checks of analyticity, harmonicity and conformality, the complex potential and velocity of flow past a cylinder with circulation, the Joukowski map, the Kutta condition, the surface pressure and the Kutta–Joukowski lift, and a contour integrator checked against the residue theorem. 3. Build a Joukowski airfoil from a circle, apply the Kutta condition, and show that the lift from integrating the surface pressure equals ρUΓ to 0.2 %. 4. Decide what the potential-flow model can and cannot say about a real airfoil, and name the measurement that exposes each limit. |
| **Prerequisites** | `P02.07` (complex numbers, Euler's formula, the complex plane), `P02.08` (the Laplace transform, for the residue cross-link), `P02.11` (Laplace's equation and harmonic functions), `P01.12` (line integrals, circulation and flux, the divergence and Stokes theorems), `P01.09` (partial derivatives, the chain rule), `P00.09` (SymPy, for checking derivatives), `P00.03`, `P00.07`. Skill check: you can differentiate log z and 1/z, and evaluate e^{iθ}. |
| **Engineering question** | A circle of radius 0.276 m is mapped by z = ζ + c²/ζ (c = 0.25 m) into an airfoil of chord 1.01 m: what are its thickness and camber, what circulation does the Kutta condition demand at 4° incidence in a 50 m/s stream, what lift per span and lift coefficient follow, and does the lift obtained by integrating the surface pressure agree with ρUΓ? |
| **Assumptions and limits** | Two-dimensional, steady, incompressible (local Mach below 0.3, which for this section means U below about 65 m/s), inviscid and irrotational flow; the Kutta condition supplies the circulation that viscosity would supply in reality; the trailing edge is a cusp (zero included angle), which real sections do not have; no boundary layer, so no drag (d'Alembert's paradox), no stall, and a lift slope that does not fall off at high incidence; the lift coefficient is per unit span of an infinite wing. Numerical derivatives are central differences with step 10⁻⁶ and are accurate to about 10⁻⁹. |
| **Analytical work** | The Cauchy–Riemann equations; harmonicity of u and v; the complex potential F = φ + iψ and the complex velocity w = dF/dz = u − iv; uniform flow, doublet and vortex; the circle theorem; the Joukowski map, its critical points and the doubling of angles there; the Kutta condition Γ = 4πUR sin(α + β); the velocity transformation w_z = w_ζ/(dz/dζ) and its finite limit at the trailing edge by l'Hôpital; the Kutta–Joukowski theorem from the residue theorem (Blasius); C_l = 2π sin(α + β)(4R/chord). |
| **Computational deliverable** | `starter/main.py`: `cauchy_riemann_residual`, `complex_derivative`, `laplacian_of_real_part`, `conformality_defect`, `cylinder_potential`, `cylinder_velocity`, `kutta_circulation`, `stagnation_points`, `joukowski`, `joukowski_derivative`, `airfoil_circle`, `airfoil_coordinates`, `airfoil_geometry`, `airfoil_surface_velocity`, `pressure_coefficient`, `lift_per_span`, `lift_coefficient`, `circle_path`, `contour_integral`, `circulation_and_flux`. Tests prove: the Cauchy–Riemann residual vanishes to 10⁻⁷ for z², eᶻ, 1/z, sin z and is 2 for z̄; Re f is harmonic to 10⁻⁵; numerical derivatives match SymPy to 10⁻⁶ and the analytic cylinder velocity to 10⁻¹²; the Joukowski map preserves angles to 10⁻⁶ away from ±c and doubles them at c; the circle |ζ| = c maps to the slit [−2c, 2c]; the airfoil is closed with its trailing edge at exactly 2c; the circle is a streamline (ψ constant to 10⁻¹⁰); the Kutta circulation puts a stagnation point at the trailing edge for four incidences; the trailing-edge speed is finite with the Kutta circulation and infinite without; the pressure-integrated lift equals ρUΓ to 0.2 % and the drag is zero to 0.2 % of the lift; the lift slope lies between 2π and 1.2 × 2π per radian; contour integrals reproduce 2πi × (sum of residues) to 10⁻¹⁰ and the velocity contour integral returns −Γ + 0i. |
| **Practical activity** | **computational experiment**: sweep the circle centre to map thickness and camber against lift slope and zero-lift angle, sweep incidence for C_l(α), and compare the Joukowski lift slope with the thin-airfoil value 2π and its thickness correction. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `joukowski_airfoil.png`, the geometry and lift table, the C_l(α) table against 2π(α + β)(4R/chord), and a short note stating which assumption of the potential-flow model you would test first and how. |
| **Assessment** | `assessments/phase-02/grading-rubric.md`; competencies 1–4. Pass: the derivation of the Cauchy–Riemann equations and of the Kutta circulation is complete; all public tests pass; the hand values of Γ, L′ and C_l agree with the code to three figures; the note names a limit with a number. |
| **References** | `@mit-ocw-18-04-complex-variables` Topic 1: Complex Algebra and the Complex Plane (sections 1.1–1.6.2), Topic 2: Analytic Functions (2.1–2.5, 2.6–2.9), Topic 5: Harmonic Functions (5.1–5.4), Topic 6: Two-Dimensional Hydrodynamics and Complex Potentials (6.1–6.4, 6.5–6.7), Topic 8: Residue Theorem, Topic 10: Conformal Transformations (10.1–10.4, 10.5–10.8) (required-free; PDF notes). Topic 3: Cauchy's Theorem, Topic 4: Cauchy Integral Formula, Topic 9: Definite Integrals Using the Residue Theorem, Topic 12: Laplace Transform (supplementary-free). The Joukowski airfoil construction itself is original instruction in this lesson. |
| **Extensions** | ★ Replace the Joukowski map by the Kármán–Trefftz map z = nc[(ζ + c)ⁿ + (ζ − c)ⁿ]/[(ζ + c)ⁿ − (ζ − c)ⁿ] with n = 2 − τ/π to give the trailing edge a finite angle τ, and show the Kutta condition still works. ★ Compute the pitching moment about the quarter-chord from the C_p distribution and compare with the thin-airfoil value. ★ Apply the Prandtl–Glauert factor 1/√(1 − M²) to the C_p distribution and estimate the critical Mach number of the section. |
| **Time** | unmeasured |

## 1. The engineering problem

Before computers, the only way to know the pressure distribution on a wing section was to find a shape for which the equations could be solved exactly, and the only shapes with exact solutions were those that could be *mapped* from a circle. The Joukowski transformation does that: it takes the flow past a cylinder, which is three terms of algebra, and bends the plane so that the cylinder becomes a shape with a rounded nose and a sharp tail. Every quantity in the airfoil plane, streamline, velocity, pressure, lift, is inherited from the cylinder plane by one chain rule. The method is a century old, and it is still how airfoil codes are validated: a panel method or a Navier–Stokes solver that cannot reproduce the Joukowski pressure distribution at low Mach number is wrong.

The concrete problem: a circle of radius 0.276 m centred at (−0.025, 0.025) m, passing through ζ = c = 0.25 m, is mapped by z = ζ + c²/ζ. Find the airfoil's chord, thickness and camber. At 4° incidence in a 50 m/s stream at sea level, apply the **Kutta condition** (the flow must leave the sharp trailing edge smoothly) to fix the circulation, and from it the lift per span and C_l. Then check the lift two independent ways: by integrating the surface pressure around the airfoil, and by the **Kutta–Joukowski theorem** L′ = ρUΓ, which is a residue calculation.

The mathematics needed is the theory of **analytic functions**: functions of a complex variable that have a derivative. Their real and imaginary parts are harmonic (`P02.11`), which is exactly the property of the velocity potential and stream function of an incompressible irrotational flow; their compositions are analytic, which is why mapping works; and their contour integrals are determined by residues, which is why the lift is ρUΓ.

Success criteria, written down first: our numerical Cauchy–Riemann check distinguishes analytic from non-analytic functions with a residual ratio of 10⁷; the conformality defect of the map is below 10⁻⁶ rad away from the critical points and exactly doubles angles at them; with the Kutta circulation the complex velocity at the trailing-edge image point is zero to 10⁻⁹ U; the integrated-pressure lift equals ρUΓ to 0.2 % with zero drag to the same tolerance; the contour integral of 1/z returns 2πi to 10⁻¹⁰; and the lift slope exceeds 2π per radian by the thickness correction, within 1 %.

## 2. Understand and derive

### 2.1 Analytic functions and the Cauchy–Riemann equations

A function f(z) = u(x, y) + i v(x, y) of z = x + iy is **analytic** at a point if the limit f′(z) = lim [f(z + Δz) − f(z)]/Δz exists and is the same for every direction of approach Δz (`@mit-ocw-18-04-complex-variables` Topic 2). Approach along the real axis (Δz = h): f′ = u_x + i v_x. Approach along the imaginary axis (Δz = ih): f′ = (u_y + i v_y)/i = v_y − i u_y. Equating real and imaginary parts:

**u_x = v_y, u_y = −v_x**: the **Cauchy–Riemann equations**.   (1)

Packaged as one statement: f_x + i f_y = 0, which is what `cauchy_riemann_residual` evaluates by central differences. For f = z² it returns 5 × 10⁻¹⁰ (rounding at h = 10⁻⁶); for f = z̄ = x − iy it returns exactly 2 (f_x = 1, f_y = −i, so f_x + i f_y = 2): the conjugate is not analytic anywhere, and neither is |z|² (residual 5 at the test points). Analytic functions are rare in the space of all functions of two variables, and that is what makes them useful: two real functions tied together by (1) are one function of z, with all of calculus available.

Differentiate (1) again: u_xx = v_yx and u_yy = −v_xy, so **u_xx + u_yy = 0**, and likewise for v (Topic 5). The real and imaginary parts of an analytic function are **harmonic**: they solve Laplace's equation (`P02.11`). `laplacian_of_real_part` checks it with the five-point stencil: 4 × 10⁻⁸ for eᶻ, exactly 4 for |z|² (whose Laplacian is 4). The converse also holds: every harmonic u has a **conjugate harmonic** v, unique up to a constant, making u + iv analytic. Level curves of u and v cross at right angles, because ∇u·∇v = u_x v_x + u_y v_y = u_x v_x − v_x u_x = 0 by (1).

### 2.2 The complex potential

For a two-dimensional, incompressible, irrotational flow (`P01.12`): irrotational means a velocity potential φ exists with (u, v) = ∇φ; incompressible means ∇·(u, v) = 0, so φ is harmonic. Incompressible also means a **stream function** ψ exists with u = ψ_y, v = −ψ_x, and irrotational makes ψ harmonic too. Together, φ_x = ψ_y and φ_y = −ψ_x: the Cauchy–Riemann equations. So

**F(z) = φ + iψ is analytic, and w = dF/dz = φ_x + iψ_x = u − iv**   (2)

is the **complex velocity** (Topic 6). Streamlines are the curves ψ = constant; a solid boundary is a streamline; the speed is |w|; and the flow rate between two streamlines is the difference of their ψ. Three building blocks, each analytic away from its singularity:

- uniform stream at angle α: F = U z e^{−iα}, w = U e^{−iα};
- **doublet** of strength μ at the origin: F = μ/z;
- **vortex** of circulation Γ at the origin: F = (iΓ/2π) log z, w = iΓ/(2πz). On the circle |z| = r this gives a tangential speed Γ/(2πr); with the + sign used here, Γ > 0 is a *clockwise* circulation, the sense that produces upward lift in a stream moving to the right.

Superposition is allowed because (2) is linear. The **cylinder** of radius R comes from the stream plus a doublet of strength UR²e^{iα}: the circle |z| = R is then a streamline (ψ = 0 on it, which the test checks to 10⁻¹⁰), the normal velocity vanishes there, and at the top of the cylinder at α = 0 the speed is 2U. Add the vortex and, with the centre at ζ₀,

**F(ζ) = U[(ζ − ζ₀)e^{−iα} + R²e^{iα}/(ζ − ζ₀)] + (iΓ/2π) log(ζ − ζ₀),**   (3)

**w(ζ) = U[e^{−iα} − R²e^{iα}/(ζ − ζ₀)²] + iΓ/(2π(ζ − ζ₀)).**   (4)

The circle stays a streamline for every Γ, because the vortex adds only tangential velocity on it: potential flow past a cylinder is a one-parameter family, and nothing in the inviscid equations chooses Γ. That is the gap the Kutta condition fills.

### 2.3 Conformal maps and the Joukowski transformation

If z = g(ζ) is analytic with g′ ≠ 0, then F(ζ(z)) is analytic in z, so the flow in the ζ plane *is* a flow in the z plane: streamlines map to streamlines, and the complex velocity transforms by the chain rule,

**w_z = dF/dz = (dF/dζ)(dζ/dz) = w_ζ / g′(ζ).**   (5)

Such a map is **conformal**: it preserves the angle between any two curves, because near ζ₀ it is multiplication by the complex number g′(ζ₀), a rotation and a scaling (`@mit-ocw-18-04-complex-variables` Topic 10). `conformality_defect` measures it: 7 × 10⁻⁹ rad for the Joukowski map at an ordinary point, but at a **critical point**, where g′ = 0, the map behaves like (ζ − ζ₀)² and angles double: a right angle at ζ = c becomes a straight angle, which is why the airfoil's trailing edge is a cusp. The test demands π/2 there, from an input of π/2, to 10⁻⁴.

The **Joukowski map** is

**z = ζ + c²/ζ, dz/dζ = 1 − c²/ζ², critical points ζ = ±c.**   (6)

Three facts to check by hand. (i) z → ζ as |ζ| → ∞: the far field is undisturbed, so the free stream in the airfoil plane is still U at angle α. (ii) ζ and c²/ζ have the same image. (iii) On the circle ζ = c e^{iθ}, z = 2c cos θ: the circle of radius c maps to the slit [−2c, 2c], a flat plate of chord 4c. A circle slightly larger than c, offset so that it still passes through ζ = c, maps to a closed curve that is rounded where the circle is far from ±c and sharp at z = 2c, where it passes through the critical point: an airfoil. Offsetting the centre to the left (x₀ < 0) thickens it; offsetting upward (y₀ > 0) cambers it. With c = 0.25 m and ζ₀ = (−0.025 + 0.025i) m the radius is R = |c − ζ₀| = 0.2761 m and the code reports chord 1.0084 m (4.03c), thickness 11.9 % and camber 4.5 % of chord.

### 2.4 The Kutta condition

Map the cylinder flow (3) through (6). By (5), w_z = w_ζ/(1 − c²/ζ²). At the trailing edge ζ = c the denominator vanishes, so unless the numerator also vanishes the speed there is infinite: the flow is trying to go round the sharp edge from the lower surface to the upper. Without circulation the code finds a speed of 416 m/s one point away from the edge and ∞ at it. Real fluids cannot do this; a boundary layer separates, and the wake adjusts the circulation until the flow leaves the edge smoothly. The inviscid model encodes that as the **Kutta condition**: *the rear stagnation point of the cylinder flow must sit at ζ = c*.

Write the trailing-edge point relative to the centre as c − ζ₀ = R e^{−iβ}, with **β = atan2(y₀, c − x₀)** = 5.19° here. Setting w(ζ = c) = 0 in (4) and multiplying by e^{−iβ}:

U[e^{−i(α+β)} − e^{i(α+β)}] + iΓ/(2πR) = 0 ⟹ −2iU sin(α + β) + iΓ/(2πR) = 0,

**Γ = 4πUR sin(α + β).**   (7)

For U = 50 m/s, R = 0.2761 m, α = 4°: Γ = 4π × 50 × 0.2761 × sin 9.19° = **27.72 m²/s**. The code confirms |w_ζ(c)| = 1.5 × 10⁻¹⁴ m/s and places the two stagnation points at ζ = 0.2500 (the edge) and −0.2938 − 0.0380i (near the nose, mapped to z = −0.503 − 0.011i, just under the leading edge, which is where a lifting airfoil's front stagnation point belongs). With both numerator and denominator of (5) zero at the edge, l'Hôpital gives a finite limit, w_z = (dw_ζ/dζ)/(d²z/dζ²): 44.7 m/s here, continuous with the neighbouring points. The zero-lift angle is α₀ = −β = −5.19°, set by camber alone.

### 2.5 Pressure and lift

Bernoulli for steady incompressible flow gives p + ½ρ|w|² = p_∞ + ½ρU², so the **pressure coefficient** is

**C_p = (p − p_∞)/(½ρU²) = 1 − |w|²/U².**   (8)

C_p = 1 at a stagnation point and negative wherever the flow is faster than the free stream. On the airfoil at 4°, the maximum speed is 79.3 m/s at x = −0.48 m on the upper surface just behind the nose, C_p = −1.51 there, and C_p = 1.00 at the front stagnation point. The force per span on the airfoil is F = −∮ p **n** ds; writing it as a complex number with the outward normal n ds = −i dz on a counter-clockwise contour, F_x + iF_y = i∮ p dz = i·½ρU² ∮ C_p dz (the constant p_∞ integrates to zero round a closed curve). The lift is the component normal to the free stream, L′ = F_y cos α − F_x sin α: **1698.0 N/m** from the 2000-point trapezoidal sum. The drag component is 10⁻⁴ N/m: **d'Alembert's paradox**, zero drag in potential flow, which is a correct statement about the model and a wrong statement about air.

### 2.6 The Kutta–Joukowski theorem is a residue

The same force is given by **Blasius' formula**, F_x − iF_y = (iρ/2)∮ w² dz, obtained by substituting (8) into the pressure integral and using the fact that w dz is real on the surface (the surface is a streamline). Because w is analytic outside the airfoil, the contour can be expanded to a large circle, where w has the Laurent expansion w = U e^{−iα} + iΓ/(2πz) + O(1/z²), so w² = U²e^{−2iα} + iΓUe^{−iα}/(πz) + O(1/z²). The **residue theorem** (`@mit-ocw-18-04-complex-variables` Topic 8) says ∮ f dz = 2πi × (sum of residues inside); only the 1/z term contributes, with residue iΓUe^{−iα}/π:

F_x − iF_y = (iρ/2)(2πi)(iΓUe^{−iα}/π) = −iρUΓ e^{−iα} ⟹ **L′ = ρUΓ, D′ = 0.**   (9)

ρUΓ = 1.225 × 50 × 27.72 = **1698.0 N/m**, the same number to five figures as the pressure integral of 2.5, which used no residues at all. The lift coefficient is C_l = L′/(½ρU² chord) = 2Γ/(U chord) = **1.100**; from (7), C_l = 2π sin(α + β) × (4R/chord), and 4R/chord = 1.095 for this section against the thin-airfoil slope 2π. The lift slope over 0–8° is 6.79 per radian, 8 % above 2π; the classical thickness correction 2π(1 + 0.77 t/c) gives 6.86, so the Joukowski section's slope is explained by its thickness to 1 %.

`contour_integral` checks the residue machinery on known integrands: ∮ dz/z = 2πi to 10⁻¹², ∮ dz/(z − 3) = 0 for a pole outside, ∮ z² dz = 0 for an entire function, ∮ 3 dz/(z − ½i) = 6πi. Applied to the velocity (4) around any circle enclosing the cylinder, ∮ w dz = Γ_ccw + iQ returns −Γ (the clockwise Γ seen counter-clockwise) and Q = 0 (no source); around the mapped airfoil the same −27.723 + 0i, independent of the contour, as Cauchy's theorem requires. The residue theorem is also what inverts Laplace transforms in `P02.08`: the Bromwich integral is a contour integral closed on the left, and each pole of the transfer function contributes its residue, which is the partial-fraction term.

### 2.7 Where the model stops being valid

- **Inviscid.** No boundary layer, so no skin friction, no form drag, no separation and no stall: the lift keeps rising linearly past the 12–16° where a real section of this thickness stalls. The Kutta condition is viscosity's only appearance.
- **Incompressible.** The peak speed is 1.59U; at U = 65 m/s the local Mach number reaches 0.3 and compressibility corrections begin (extension 3).
- **Cusped trailing edge.** Real sections have a finite edge angle; the Kármán–Trefftz map (extension 1) fixes that.
- **Two-dimensional.** An infinite span. Finite wings shed trailing vorticity that reduces the lift slope and adds induced drag (Phase 10).
- **One family of shapes.** The Joukowski family has its maximum thickness far forward and a cusped tail; modern sections are designed by inverse methods, with this family as the check case.

## 3. Build the computational model

Open `starter/main.py`. `lesson_airfoil` and `main` are given.

### Step 1: analyticity by numbers

```python
def cauchy_riemann_residual(f, z, h=1e-6):
    fx = (f(z + h) - f(z - h)) / (2 * h)
    fy = (f(z + 1j * h) - f(z - 1j * h)) / (2 * h)
    return np.abs(fx + 1j * fy)                                       # (1) in one line
```

`complex_derivative` is the real-direction difference; `laplacian_of_real_part` the five-point stencil with h = 10⁻⁴ (the second difference loses precision faster, so use a larger step); `conformality_defect` compares arg[(f(z + he^{iθ₂}) − f(z))/(f(z + he^{iθ₁}) − f(z))] with θ₂ − θ₁, both reduced to (−π, π].

### Step 2: the cylinder

`cylinder_potential` and `cylinder_velocity` are (3) and (4) with s = ζ − ζ₀. `kutta_circulation` is (7). `stagnation_points` solves w = 0 as the quadratic s²e^{−iα} + iΓs/(2πU) − R²e^{iα} = 0 and adds ζ₀.

### Step 3: the map and the airfoil

`joukowski` and `joukowski_derivative` are (6). `airfoil_circle` returns (ζ₀, |c − ζ₀|) and rejects x₀ ≥ 0 or y₀ < 0. `airfoil_coordinates` walks the circle counter-clockwise *starting at the trailing edge* (angle −β), sets the first and last points to exactly c, and maps. `airfoil_geometry` splits the polygon at the leading edge (minimum x), interpolates both surfaces onto a common x grid, and returns chord, maximum thickness and maximum camber in metres.

### Step 4: velocity, pressure, lift

`airfoil_surface_velocity` applies (5) point by point and treats the trailing edge separately: if |w_ζ(c)| is below 10⁻⁹U it returns the l'Hôpital limit (dw_ζ/dζ)/(2/c), otherwise `inf`. `pressure_coefficient` is (8); `lift_per_span` is ρUΓ; `lift_coefficient` is 2Γ/(U chord).

### Step 5: contour integrals

`circle_path` returns n + 1 points with the last equal to the first *exactly* (assign it; `exp(2πi)` is not quite 1). `contour_integral` recognises a closed path, treats the points as samples of a periodic curve z(t), computes z′(t) by the spectral derivative of `P02.10` (`np.fft`), and returns Σ f(z_k) z′_k Δt; the trapezoidal rule on a periodic analytic integrand converges faster than any power of the spacing, which is why 2000 points give 10⁻¹². For an open path it falls back to the chord trapezoid. `circulation_and_flux` returns the real and imaginary parts of ∮ w dz.

Running `main()` prints

```
circle: centre -0.0250+0.0250j m, R = 0.2761 m, beta = 5.194 deg; map constant c = 0.25 m
Cauchy-Riemann residual of z^2     : 4.97e-10; Laplacian of Re: 0.00e+00
Cauchy-Riemann residual of conj(z) : 2.00e+00; Laplacian of Re: 0.00e+00
Cauchy-Riemann residual of |z|^2   : 5.00e+00; Laplacian of Re: 4.00e+00
conformality defect of the Joukowski map at zeta = 0.4 + 0.3i: 7.21e-09 rad; at the critical point zeta = c: 1.200 rad (angles double there)
contour integral of 1/z around |z| = 1: -0.000000000000+6.283185307180j (2 pi i = 6.283185307180 i)
airfoil: chord 1.0084 m, thickness 11.9 % chord, camber 4.5 % chord; trailing edge at z = 0.500+0.000j
Kutta circulation at alpha = 4.0 deg: Gamma = 27.723 m^2/s
stagnation points in the zeta plane: 0.2500-0.0000j, -0.2938-0.0380j; trailing edge zeta = 0.25
surface speed at the trailing edge (finite by l'Hopital): 44.686 m/s; max surface speed 79.28 m/s, min C_p -1.514, max C_p 0.9999
without circulation: speed at the trailing edge inf and near it 415.6 m/s
lift per span rho U Gamma = 1698.0 N/m; C_l = 1.0997; thin-airfoil 2 pi (alpha + beta) = 1.0083
lift from integrating C_p around the surface: 1698.0 N/m (Kutta-Joukowski 1698.0 N/m)
contour integral of w dz around the airfoil (ccw): -27.723 + 0.0e+00 i  (expected -Gamma = -27.723, flux 0)
```

and the C_l table 0.623, 0.862, 1.100, 1.336, 1.571 for α = 0, 2, 4, 6, 8°.

### Verify

- **Analytic cases:** Cauchy–Riemann residual below 10⁻⁷ for five analytic functions and equal to 2 for z̄; Laplacian of Re f below 10⁻⁵ and equal to 4 for |z|²; derivatives against SymPy to 10⁻⁶ and the closed-form cylinder velocity against SymPy's derivative of (3) to 10⁻¹²; ∮ dz/z = 2πi, ∮ 3 dz/(z − ½i) = 6πi, zero for a pole outside and for z² + eᶻ; the circle |ζ| = c maps onto [−2c, 2c] with zero imaginary part to 10⁻¹⁵.
- **Limiting cases:** z → ζ at infinity to 10⁻⁶; the far-field velocity is U e^{−iα} to 10⁻⁴; the cylinder top speed is 2U at α = Γ = 0; y₀ = 0 gives a symmetric airfoil with zero camber; Γ = 0 at α = β = 0; doubling −x₀ increases the thickness by more than half.
- **Structure and conservation:** ψ is constant on the circle to 10⁻¹⁰ and the normal velocity vanishes there; angles are preserved to 10⁻⁶ at four points and three angle pairs, and doubled at c; the velocity contour integral gives −Γ and zero flux on two radii, a source adds flux only, and a contour that misses the cylinder sees nothing; the airfoil polygon is closed with the trailing edge at 2c exactly.
- **The Kutta condition and the lift:** |w_ζ(c)| < 10⁻⁹U for four incidences and a stagnation point within 10⁻⁹ of c; trailing-edge speed finite and continuous with its neighbour to 5 % with Γ, infinite without; C_p max = 1 to 2 × 10⁻³ and min below −1; pressure-integrated lift = ρUΓ to 0.2 %; drag below 0.2 % of lift; lift slope in (2π, 1.2 × 2π).

### Compare

`scipy` and `numpy` offer no conformal-mapping or potential-flow routine, and SymPy is used only as the oracle for derivatives. What practice uses in place of the Joukowski solution is a **panel method** (distributed sources and vortices on the real section's surface, solved as a linear system; Phase 10) and, beyond that, Navier–Stokes solvers. Both are validated against exactly this calculation: a panel code run on the Joukowski coordinates from `airfoil_coordinates` must return C_l = 1.100 at 4° and the C_p curve of the figure, to the discretisation error of its panels. The value of the closed form is not that engineers draw Joukowski wings; it is that an exact answer exists for one wing, so that every approximate method can be measured against it.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Three sweeps; write each prediction first.

1. **Thickness sweep.** x₀/c = −0.05, −0.10, −0.15, −0.20 with y₀ = 0: tabulate t/c, chord, lift slope dC_l/dα (from α = 0 and 8°) and compare with 2π(1 + 0.77 t/c). Predict that the slope rises with thickness and that the symmetric sections have α₀ = 0.
2. **Camber sweep.** y₀/c = 0, 0.05, 0.10, 0.15 with x₀/c = −0.10: tabulate camber, β, C_l at α = 0 and the zero-lift angle α₀ = −β. Predict that α₀ moves by about −1° per 1 % camber (thin-airfoil theory gives about −1.1° per 1 % for a circular-arc mean line, which is close to the Joukowski mean line).
3. **Incidence sweep.** The lesson section at α = −6 … 12° in 2° steps: C_l, maximum surface speed, minimum C_p, and the free-stream speed at which the peak local Mach number would reach 0.3 (speeds scale with U). Predict the linear C_l(α) and that the suction peak grows faster than linearly.

Raw data format: three CSV tables with a one-line header naming the fixed parameters; `joukowski_airfoil.png`.

What you should find: lift slopes within 0.4 % of 2π(1 + 0.77 t/c) for all four thicknesses (6 % to 22 %), always slightly below it; zero-lift angles linear in camber at −1.15° per 1 % of chord; a linear C_l(α) with slope 6.79 per radian everywhere (no stall: the model does not know what stall is); and a peak-speed ratio that reaches 2 near 7° and 2.9 at 12°, where the local Mach number would reach 0.3 at a free-stream speed of only 35 m/s.

## 5. Validate and improve

| quantity | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| R, β | 0.2761 m, 5.19° | 0.2761 m, 5.194° | 0 | exact |
| chord | ≈ 4c = 1.00 m | 1.0084 m | 0.8 % | the offset circle is larger than c |
| Γ at 4° | 27.72 m²/s | 27.723 m²/s | 0 | exact, (7) |
| \|w_ζ\| at the trailing edge | 0 | 1.5 × 10⁻¹⁴ m/s | rounding | — |
| L′ = ρUΓ | 1698 N/m | 1698.0 N/m | 0 | exact, (9) |
| L′ from ∮ C_p dz | 1698 N/m | 1698.0 N/m | < 0.01 % | 2000-point trapezoid |
| D′ from ∮ C_p dz | 0 | 10⁻⁴ N/m | 10⁻⁷ of L′ | trapezoid |
| C_l at 4° | 2π sin 9.19° × 1.095 = 1.100 | 1.0997 | 0 | exact |
| lift slope | 2π(1 + 0.77 × 0.119) = 6.86 /rad | 6.79 /rad | 1 % | an empirical correction against an exact family |
| ∮ dz/z | 2πi | 2πi to 10⁻¹² | rounding | spectral trapezoid |

Every difference is named and every one is rounding or a stated approximation. The model's agreement with *itself* three ways (Kutta–Joukowski, pressure integral, contour integral of w) is complete; its agreement with a real section is a different question. Published measurements of 12 %-thick, 4 %-cambered sections at Reynolds numbers of a few million show lift slopes of about 0.10–0.11 per degree (5.7–6.3 per radian), 5–15 % below the 6.79 found here, because the boundary layer thickens toward the trailing edge and effectively de-cambers the section; they show a maximum C_l of about 1.4–1.6 at 14–16°, where this model gives 2.3–2.5 and keeps climbing; and they show a drag coefficient of about 0.006–0.008 at C_l = 0.4, where this model gives zero. The potential-flow pressure distribution is good to a few per cent over the front half of the section at low incidence, which is where it is used: to predict suction peaks, pitching moments and the pressure gradients that a boundary-layer calculation then needs as input.

Improvements worth making: the Kármán–Trefftz edge angle; a viscous correction by displacement thickness (Phase 11); the pitching-moment integral. Improvement worth *not* making: refining the 2000-point contour. The pressure integral already agrees with the exact theorem to the last printed figure, and the remaining disagreement with reality is physics, not resolution.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivation of (1), (7) and (9); `joukowski_airfoil.png`; the three sweep tables; and a one-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) giving the section's geometry, Γ, L′ and C_l at 4°, the three-way lift check, and the assumption you would test first.

### Adversarial self-checks

- **What if the speed doubles?** Γ doubles (7), L′ quadruples (9), C_l is unchanged, the peak local speed doubles to 159 m/s and the incompressible assumption fails (M = 0.47 locally). Which of the four statements would survive in a wind tunnel?
- **What if the incidence goes to −β?** Γ = 0, L′ = 0, the trailing-edge speed is finite without any circulation, and the stagnation points sit symmetrically about the centre line of the circle. Check all three.
- **Which assumption fails first in the real world?** Inviscid flow at the trailing edge: the boundary layer displaces the Kutta point and reduces the lift slope by 5–15 %. Then stall, which the model never predicts.
- **What measurement could falsify the model?** A pressure-tapped section in a wind tunnel: the measured C_p over the rear 20 % of chord will lie above the potential-flow curve (less pressure recovery), and the integrated lift will be lower than ρUΓ from the measured circulation around a far contour, by the amount of the boundary-layer de-cambering.
- **Does the result survive a change in resolution or contour?** The contour integral of w gives −27.723 on circles of radius 1.5R and 4R and on the mapped airfoil; the pressure integral with 400 points instead of 2000 changes the lift by less than 0.1 %. State the smallest number of points you would accept.
- **Which decision would change if the uncertainty doubled?** If the circle centre were known only to ±0.005 m (2 % of c) instead of ±0.0025 m, the camber would be uncertain by ±1 % of chord and the zero-lift angle by ±1°, which at cruise C_l = 0.4 is a 25 % uncertainty in the incidence the wing must be rigged at; you would measure the section, not compute it.

## Check yourself

1. State the Cauchy–Riemann equations and show in two lines that they make u harmonic.
2. What are φ and ψ for a flow, and why is F = φ + iψ analytic for an incompressible irrotational flow? What is dF/dz?
3. Write the complex potential of a uniform stream past a cylinder of radius R with clockwise circulation Γ, and name the term that produces each of the three features.
4. What does the Joukowski map do to the circle |ζ| = c, and what happens to angles at ζ = ±c?
5. State the Kutta condition and the circulation it fixes for a circle of radius R whose trailing-edge point lies at angle −β below the centre.
6. Why is the velocity at the trailing edge finite with the Kutta circulation although dz/dζ = 0 there?
7. Derive L′ = ρUΓ in one paragraph from Blasius' formula and the residue theorem, and say what the same theorem says about drag.
8. A section has 4R/chord = 1.10 and β = 3°. What is C_l at α = 5°, and what is the zero-lift angle?
