# P01.11  Vector Fields, Line & Surface Integrals

> Circulation is what a closed path *sees* of a vortex; flux is what a surface *sees* of a source. Both are integrals, and both are one line of NumPy.

| | |
|---|---|
| **Learning outcomes** | 1. Write the velocity field of a point vortex, a source, a uniform stream and their superposition, and show by hand that the vortex is irrotational everywhere except at its centre. 2. Derive the circulation of a point vortex around any enclosing path and show it is zero around any path that does not enclose it. 3. Implement and test a midpoint-rule line integral, a 2D flux integral and a parametric-surface flux integral, and measure their second-order convergence. 4. Compute the mass flow rate through an engine inlet as a surface integral and the work done by gravity between two orbits as a line integral, and check both against closed forms. 5. Decide from a convergence plot how many segments a circulation calculation needs for a stated accuracy. |
| **Prerequisites** | `P01.07` (vectors, dot product, parametric curves), `P01.08` (partial derivatives), `P01.10` (double integrals, Jacobians), `P00.03` (NumPy broadcasting), `P00.07` (pytest). Skill check: you can parametrise a circle as (R cos t, R sin t) and differentiate it; you can evaluate ∫₀^{2π} cos² t dt. |
| **Engineering question** | A wing in steady flight leaves behind a vortex with a circulation Γ; how much lift per metre of span does that circulation produce, and how do we measure Γ from the velocity field on *any* loop around the wing? |
| **Assumptions and limits** | Two-dimensional, steady, incompressible, inviscid flow for the aerodynamic examples (potential flow). The point vortex has infinite velocity at its centre: real vortices have a viscous core, modelled here by the Rankine vortex. Gravity is that of a point mass (spherical Earth, no J₂). The inlet mass flow assumes a uniform velocity across the inlet face; a real inlet has a boundary layer and distortion. The numerical integrals assume the field is smooth along the path or surface; near a kink or a singularity the convergence order drops. |
| **Analytical work** | Vortex, source and doublet fields from their potentials; circulation of a vortex around a circle by direct parametrisation and around any loop by the potential argument; flux of a source; the Kutta–Joukowski relation stated and dimensionally checked; work in a gravitational field from the potential; closed-form error of the midpoint rule on a circle; mass flow through a disc by hand. |
| **Computational deliverable** | `starter/main.py`: `vortex_velocity`, `source_velocity`, `uniform_velocity`, `cylinder_flow_velocity`, `rankine_vortex_velocity`, `gravity_field`, `gravity_potential`, `circle_path`, `rectangle_path`, `line_integral`, `circulation`, `flux_2d`, `work_along_path`, `kutta_joukowski_lift`, `point_source_3d`, `sphere_surface`, `disc_surface`, `surface_normal`, `surface_flux`, `mass_flow_rate`, `circulation_convergence`. Tests prove: circulation equals Γ on every enclosing path and 0 otherwise; flux of a source equals q; the uniform stream has neither; second-order convergence of both integrals with a closed-form error on circles; path independence of gravitational work; zero net flux of a uniform stream through a closed sphere; ρVA through a disc. |
| **Practical activity** | **computational experiment**: a convergence study of circulation against the number of path segments for circles and squares, a sweep of circulation against radius through the core of a Rankine wake vortex, and a mass-flow calculation for a boundary-layer velocity profile. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `circulation_vortex.png`, the convergence table with fitted slopes, the Rankine sweep table, and a short note on the inlet mass flow with its profile correction. |
| **Assessment** | `assessments/phase-01/grading-rubric.md`; competencies 1–4. Pass: the circulation derivation (both the enclosing and the non-enclosing case) is complete with assumptions named; all public tests pass; the convergence slope is measured and compared with the closed form; the judgment problem is answered with a measurable sign. |
| **References** | `@openstax-calculus-volume-3` 6.1 Vector Fields (required-free); `@openstax-calculus-volume-3` 6.2 Line Integrals (required-free); `@openstax-calculus-volume-3` 6.3 Conservative Vector Fields (required-free); `@openstax-calculus-volume-3` 6.6 Surface Integrals (required-free). |
| **Extensions** | ★ Add a doublet-only field and show the cylinder streamlines by integrating dx/dt = u, dy/dt = v with `solve_ivp`; ★ compute the pressure on the cylinder from Bernoulli and integrate it to recover the Kutta–Joukowski lift directly; ★ replace the point-mass gravity with the J₂ field and show that work is still path-independent. |
| **Time** | unmeasured |

## 1. The engineering problem

Every lifting wing drags the air into a slow rotation around itself. Seen from the side, the flow goes faster over the top than under the bottom; the difference is the same thing as a net swirl, and the swirl has a single number attached to it: the **circulation** Γ (m²/s), the integral of velocity around a closed loop enclosing the wing section. The Kutta–Joukowski theorem (Phase 10 derives it; here we take it as given and check its dimensions) says that the lift per unit span is L' = ρ V∞ Γ. So if we can compute Γ from a velocity field, we have the lift. The vortex a wing sheds persists for minutes behind a transport aircraft, and the circulation of that wake vortex is what sets the separation between landing aircraft.

The engineering question is therefore: given a velocity field **v**(x, y), compute the circulation Γ = ∮ **v** · d**r** around a loop, show that the answer does not depend on which loop we draw as long as it encloses the vortex, and turn Γ into a lift. Along the way we meet the second kind of integral over a field, the **flux** through a surface, which is how an engine's mass flow rate is defined: ṁ = ∬ ρ **v** · **n** dA over the inlet face. The same machinery gives the work done by a force field along a path, and the fact that gravity's work is path-independent is what makes orbital mechanics tractable in Phase 13.

Success criteria, written before we start: the circulation of a point vortex of strength Γ comes out as Γ to one part in 10⁵ on a circle, on a square and on an off-centre circle, and as 0 to the same tolerance on any loop that does not enclose it; the numerical integral converges at second order in the number of segments, with the error on a circle matching a closed form; the flux of a point source through a sphere equals its strength for any radius; the inlet mass flow at 10 km altitude and 250 m/s through a 2 m diameter fan face agrees with ρVA to 10⁻⁶; and the work done by gravity between a 400 km and an 800 km circular orbit agrees with the potential difference on two different paths.

## 2. Understand and derive

### 2.1 Vector fields

A **vector field** assigns a vector to every point: **v**(x, y) = (u(x, y), v(x, y)) for a planar velocity field, **F**(x, y, z) for a force field in space (`@openstax-calculus-volume-3` 6.1). Three fields carry this lesson, all in SI units and all written with a centre at the origin.

**Uniform stream** of speed V∞ along +x: **v** = (V∞, 0). Dimensions: m/s.

**Point vortex** of circulation Γ (m²/s), positive counter-clockwise:

**v** = (Γ/2π) (−y, x)/r², with r² = x² + y², i.e. a tangential speed v_θ = Γ/(2πr).   (1)

Dimensions: (m²/s)(m)/(m²) = m/s. Good. The speed falls as 1/r and becomes infinite at the centre; that is the model's known failure, fixed in §2.6.

**Point source** of volume flow rate q per unit depth (m²/s):

**v** = (q/2π) (x, y)/r², a radial speed v_r = q/(2πr).   (2)

Figure 1 shows the vortex and the source side by side; the arrows of the vortex go round, those of the source go out:

![Two velocity fields side by side. Left, a point vortex with positive circulation Γ: arrows go counter-clockwise around the centre with tangential speed v_θ = Γ/(2πr). Right, a point source with positive strength q: arrows point radially away from the centre with radial speed v_r = q/(2πr).](figures/vortex-and-source.svg)
*Figure 1. A point vortex and a point source. Both speeds fall as 1/r; the vortex’s velocity goes round the centre, the source’s goes straight out.*

Both are **gradient fields** away from the origin: the source is ∇φ with φ = (q/2π) ln r, and the vortex is ∇φ with φ = Γθ/(2π), where θ = atan2(y, x) is the polar angle. Check the vortex: ∂θ/∂x = −y/r², ∂θ/∂y = x/r², so ∇(Γθ/2π) = (Γ/2π)(−y, x)/r², which is (1). The catch, and it is the whole point of §2.3, is that θ is not a single-valued function on a loop around the origin: go round once and θ has increased by 2π.

**Assumption 1 (potential flow):** the flow is steady, incompressible and inviscid, so velocity fields superpose. The flow past a circular cylinder of radius a with circulation is the sum of the stream, a **doublet** (a source and a sink pushed together) of strength 2πV∞a² and a vortex:

u = V∞ [1 − a²(x² − y²)/r⁴] − Γy/(2πr²),  v = −2V∞a² x y/r⁴ + Γx/(2πr²).   (3)

On r = a the radial velocity is V∞(1 − a²/r²) cos θ = 0, so no fluid crosses the body: that is the test `test_cylinder_flow_has_no_normal_velocity_on_the_surface`. Far away the doublet and vortex terms vanish as 1/r² and 1/r and the stream is recovered.

### 2.2 The line integral

Take a curve C parametrised by **r**(t) for t from t₀ to t₁ (`P01.07`). The **line integral of a vector field along C** is

∫_C **F** · d**r** = ∫_{t₀}^{t₁} **F**(**r**(t)) · **r**′(t) dt.   (4)

(`@openstax-calculus-volume-3` 6.2.) It adds up the component of **F** along the direction of travel, weighted by distance. When **F** is a force, (4) is the work done by the field on a body moving along C (J). When **F** is a velocity, (4) around a *closed* loop is the circulation (m²/s). Reversing the direction of travel changes the sign; that is why every loop in this lesson is counter-clockwise, and the tests check it.

**Worked example: circulation of a point vortex around a circle of radius R centred on it.** Parametrise **r**(t) = (R cos t, R sin t), 0 ≤ t ≤ 2π, so **r**′(t) = (−R sin t, R cos t). From (1) at r = R, **v** = (Γ/2πR)(−sin t, cos t). The dot product is (Γ/2πR)(R sin² t + R cos² t) = Γ/2π, a constant, and

Γ_circle = ∫₀^{2π} (Γ/2π) dt = Γ.   (5)

Independent of R, which is the first hint that the circulation is a property of what is *inside* the loop, not of the loop.

**Worked example: the same loop with a source.** **v** = (q/2πR)(cos t, sin t), and **v** · **r**′ = (q/2πR)(−R cos t sin t + R sin t cos t) = 0. A source has no circulation. By the same algebra the vortex has no *flux* (§2.4). The two fields are each other's rotation by 90°, and the two integrals are each other's rotation by 90°.

### 2.3 Why any enclosing loop gives Γ, and any other loop gives zero

Now the general statement, which is the heart of the lesson. Suppose a field is a gradient, **F** = ∇φ, on a region containing the curve C from point A to point B. Then by the chain rule d/dt [φ(**r**(t))] = ∇φ · **r**′(t), so

∫_C ∇φ · d**r** = ∫_{t₀}^{t₁} d/dt [φ(**r**(t))] dt = φ(B) − φ(A).   (6)

This is the **fundamental theorem for line integrals** (`@openstax-calculus-volume-3` 6.3): the integral depends only on the endpoints. For a closed loop B = A and the integral is zero. Such a field is called **conservative**, and the function φ a **potential**.

Apply it to the vortex with φ = Γθ/2π. On any loop that does *not* enclose the origin, you can choose the branch of θ continuously all the way round (the loop never crosses a ray from the origin twice in the "wrong" way), so φ returns to its starting value and the circulation is zero. On a loop that *does* enclose the origin, θ must increase by exactly 2π in one counter-clockwise turn, so φ increases by Γ and the circulation is Γ, whatever the loop's shape. A loop that winds twice gives 2Γ. This argument needs only that the field is a gradient *away from the centre*; the centre itself is excluded, and that single excluded point is what makes the circulation non-zero. In `P01.12` the same fact reappears as "the vorticity is concentrated at the centre" and Stokes' theorem.

**Assumption 2 (the loop avoids the singularity):** a path that passes through the vortex centre has no defined circulation; a numerical path that passes *near* it needs many segments (§3, Verify).

### 2.4 Flux through a curve and through a surface

For a planar field and a closed counter-clockwise loop, the **outward flux** is

∮_C **v** · **n** ds, with **n** ds = (dy, −dx),   (7)

because rotating the tangent element (dx, dy) clockwise by 90° gives the outward normal. For the source on a circle of radius R: **v** · **n** = q/(2πR) and the flux is q, independent of R. For the vortex it is zero. In words: flux counts how much fluid leaves the loop per unit time per unit depth, and only a source inside contributes.

In three dimensions the surface takes over the role of the loop (`@openstax-calculus-volume-3` 6.6). Parametrise a surface S by **r**(s, t) over a rectangle in the (s, t) plane. The two tangent vectors **r**_s and **r**_t span a parallelogram of area |**r**_s × **r**_t| ds dt, and the cross product itself is a normal vector. The flux integral is

∬_S **F** · **n** dA = ∬ **F**(**r**(s, t)) · (**r**_s × **r**_t) ds dt.   (8)

The *order* of the cross product fixes which way is "positive". For the sphere **r**(s, t) = R(sin s cos t, sin s sin t, cos s) with s the polar angle, **r**_s × **r**_t = R² sin s (sin s cos t, sin s sin t, cos s), which points outward and has the familiar area element R² sin s ds dt. For the disc **r**(s, t) = (s cos t, s sin t, 0), **r**_s × **r**_t = (0, 0, s): the +z normal with the polar area element s ds dt.

**Worked example: the mass flow through an engine inlet.** The fan face of a large turbofan is a disc of radius R = 1 m. At 10 km in the ISA, ρ = 0.4135 kg/m³ (`P00.03`), and at a cruise speed of 250 m/s with a uniform velocity across the face, the mass flux vector is ρ**v** = (0, 0, ρV). By (8) with the disc parametrisation, ṁ = ∫₀^{2π}∫₀^R ρV s ds dt = ρVπR² = 0.4135 × 250 × π = 324.76 kg/s. Dimensions: (kg/m³)(m/s)(m²) = kg/s. Three hundred kilograms of air a second is the right order for a large engine at cruise. The code reproduces this to 10⁻⁶ relative, which only proves the disc normal is right; the physics content is the uniform-profile assumption, which problem 5 relaxes.

**Worked example: the flux of a 3D point source through a sphere.** For **v** = q **r**/(4π|**r**|³) and the sphere of radius R centred on it, **v** · **n** = q/(4πR²) and the flux is q, for every R. Through a sphere that does not contain the source the flux is zero. This is the three-dimensional version of §2.3, and in `P01.12` it becomes the divergence theorem.

### 2.5 Work done by gravity

The gravitational acceleration of a point mass with parameter μ = GM is **g** = −μ **r**/|**r**|³, and it is minus the gradient of φ = −μ/r. By (6) the work per unit mass done by gravity along any path from A to B is φ(A) − φ(B) = μ(1/r_B − 1/r_A). Moving a satellite from a 400 km circular orbit (r_A = 6771 km) to 800 km (r_B = 7171 km): with μ = 3.986 × 10¹⁴ m³/s², the work is −3.2837 × 10⁶ J/kg: gravity does negative work, the spacecraft must supply it. The number does not depend on whether you go straight out or spiral; a test checks two paths.

### 2.6 Where the model stops being valid

- **The point vortex.** Speed → ∞ as r → 0, and so does the kinetic energy. A real vortex has a viscous core of radius a inside which the fluid rotates like a solid body: the **Rankine vortex** v_θ = Γr/(2πa²) for r < a, Γ/(2πr) for r ≥ a. Circulation around a circle of radius R < a is then 2πR · ΓR/(2πa²) = Γ(R/a)²: it grows with the enclosed area until the loop contains the whole core. The loop now *does* matter, because the vorticity is spread out rather than concentrated at a point. A transport aircraft on approach (weight 687 kN, span 34 m, 70 m/s) sheds a pair of vortices with Γ ≈ 4W/(πρV∞b) = 300 m²/s and a core of a metre or two; the peak swirl at a 1.5 m core radius is Γ/(2πa) = 32 m/s. Those are the numbers that set wake-turbulence separations.
- **Potential flow.** No viscosity, so no drag and no boundary layer; the circulation of a real airfoil is fixed by the Kutta condition at the trailing edge, which potential flow cannot supply by itself. Lift from L' = ρV∞Γ is good to a few percent below stall; the model says nothing about stall.
- **Uniform inlet velocity.** A real inlet has a boundary layer on its walls and a non-uniform profile; problem 5 shows a 1/7-power profile carries 18 % less mass than the uniform assumption.
- **Point-mass gravity.** Earth's oblateness (J₂) changes the potential by about 10⁻³ at low orbit; the field is still conservative, so (6) still holds with the corrected φ.
- **Smoothness of the integrand.** The midpoint rule is second order for smooth fields. Across a kink (the Rankine core edge) or near a singularity the constant in front of h² grows, and for a profile with an infinite derivative (problem 5) the order itself drops.

## 3. Build the computational model

Open `starter/main.py`. Fields are plain functions `field(x, y) -> (u, v)` that accept arrays; paths are `(n, 2)` arrays of vertices; surfaces are functions `surface(s, t) -> (x, y, z)`. Keeping fields as callables means a superposition is just a lambda that adds two of them.

### Step 1: the fields

```python
def vortex_velocity(x, y, gamma, center=(0.0, 0.0)):
    dx = np.asarray(x, dtype=float) - center[0]
    dy = np.asarray(y, dtype=float) - center[1]
    r2 = dx**2 + dy**2
    return -gamma * dy / (2.0 * np.pi * r2), gamma * dx / (2.0 * np.pi * r2)
```

Equation (1), element by element. `source_velocity` is (2) with the signs swapped and no minus. `cylinder_flow_velocity` is (3): write the doublet terms with `r4 = r2**2` and add `vortex_velocity`. The Rankine vortex needs `np.where` for the two branches and a guard for r = 0, where the velocity is zero (not `nan`): the same "compute both branches safely, then select" idiom as `P00.03`.

### Step 2: paths

```python
def circle_path(radius, n, center=(0.0, 0.0)):
    th = 2.0 * np.pi * np.arange(n) / n
    return np.column_stack([center[0] + radius * np.cos(th), center[1] + radius * np.sin(th)])
```

The first vertex is not repeated: a closed polygon with n vertices has n segments, the last from vertex n−1 back to vertex 0. `rectangle_path` walks the four sides counter-clockwise, `n_per_side` vertices each, corners included once. Reversing a path (`points[::-1]`) must flip the sign of every integral.

### Step 3: the line integral by the midpoint rule

```python
def _segments(points, closed):
    pts = np.asarray(points, dtype=float)
    nxt = np.roll(pts, -1, axis=0) if closed else pts[1:]
    cur = pts if closed else pts[:-1]
    return 0.5 * (cur + nxt), nxt - cur

def line_integral(field, points, closed=True):
    mid, d = _segments(points, closed)
    u, v = field(mid[:, 0], mid[:, 1])
    return float(np.sum(u * d[:, 0] + v * d[:, 1]))
```

This is (4) with the parametrisation replaced by the polygon: on each segment **r**′ dt becomes the displacement Δ**r**, and the field is sampled once, at the midpoint. `np.roll(pts, -1, axis=0)` shifts every vertex up by one so that row i holds vertex i+1 and the last row holds vertex 0: that closes the loop without an `if`. `circulation` is `line_integral` with `closed=True`; `work_along_path` with `closed=False`; `flux_2d` is the same loop with (7): `u * d[:, 1] - v * d[:, 0]`.

### Step 4: the circle error has a closed form

On a circle of radius R with n segments, the midpoint of each chord sits at radius R cos(π/n), where the vortex speed is Γ/(2πR cos(π/n)); the chord has length 2R sin(π/n) and is parallel to the velocity at its midpoint. Each segment contributes Γ tan(π/n)/π and the sum is

Γ_n = Γ (n/π) tan(π/n) = Γ [1 + π²/(3n²) + …].   (9)

So the relative error is π²/(3n²): 0.8 × 10⁻³ at n = 64, 0.5 × 10⁻⁵ at n = 512. That is second-order convergence with a known constant, and `test_circulation_converges_at_second_order` checks (9) to 10⁻⁸. For a reference value you did not derive you would only see the *slope*; having the constant too is a stronger check, and it is why the circle is the first test case.

### Step 5: surfaces

```python
def surface_normal(surface, s, t, h=1e-6):
    xs1, ys1, zs1 = surface(s + h, t)
    xs0, ys0, zs0 = surface(s - h, t)
    xt1, yt1, zt1 = surface(s, t + h)
    xt0, yt0, zt0 = surface(s, t - h)
    r_s = np.stack([xs1 - xs0, ys1 - ys0, zs1 - zs0], axis=-1) / (2.0 * h)
    r_t = np.stack([xt1 - xt0, yt1 - yt0, zt1 - zt0], axis=-1) / (2.0 * h)
    return np.cross(r_s, r_t)
```

The tangents **r**_s and **r**_t are central differences in the parameters (`P01.01`), so any surface you can write as a function works without deriving its tangents by hand; the error of the central difference is ~h² ≈ 10⁻¹², far below the integration error. `np.stack(..., axis=-1)` builds an array whose last axis holds the three components, which is what `np.cross` expects.

```python
def surface_flux(field, surface, s_range, t_range, ns, nt):
    ds = (s_range[1] - s_range[0]) / ns
    dt = (t_range[1] - t_range[0]) / nt
    s = s_range[0] + ds * (np.arange(ns) + 0.5)
    t = t_range[0] + dt * (np.arange(nt) + 0.5)
    S, T = np.meshgrid(s, t, indexing="ij")
    x, y, z = surface(S, T)
    n = surface_normal(surface, S, T)
    u, v, w = field(x, y, z)
    return float(np.sum(u * n[..., 0] + v * n[..., 1] + w * n[..., 2]) * ds * dt)
```

Equation (8) with the midpoint rule in both parameters: cell centres at `(i + 0.5) * ds`, one field sample per cell, times the cell area in parameter space. The `ns × nt` grid is built by `meshgrid` with `indexing="ij"` so that `S[i, j] = s[i]`.

### Verify

Run `python -m pytest curriculum/phase-01/lesson-11-line-and-surface-integrals` from the repo root. In order of how much they would hurt if wrong:

1. **Field identities.** The vortex has no radial component and speed Γ/(2πr); the source has no tangential component; the cylinder has zero normal velocity on its surface and recovers V∞ far away; the Rankine vortex is continuous at r = a and finite at r = 0; gravity is −∇φ to 10⁻⁶ by central differences.
2. **Orientation.** Both path helpers give a positive shoelace area (counter-clockwise); reversing a path negates the circulation.
3. **The topological fact.** Circulation = Γ on three different enclosing loops to 10⁻⁵, and 0 on two non-enclosing loops. Flux of a source = q on loops of any size; a source has no circulation and a vortex no flux; a uniform stream has neither.
4. **Convergence.** Rates above 1.8 on circles and squares, and the circle error equals (9) to 10⁻⁸.
5. **Conservative field.** Gravity's work between the two orbits agrees with φ(A) − φ(B) on a straight path and on an arc to 10⁻⁶; the circulation of gravity around a loop is zero; the open line integral of ∇(x²y + sin y) along a wiggly curve equals the potential difference.
6. **Surfaces.** The sphere normal points outward with magnitude R² sin s; a source's flux through a sphere is q for three radii and centre offsets and 0 when the source is outside; a uniform stream has zero net flux through a closed sphere and ρVπR² through the disc; the flux of the position field **r** through a sphere (exactly 4πR³) converges at second order.

### Compare

There is no circulation function in SciPy, but there is a general quadrature, and comparing against it tells you what the midpoint rule costs. `scipy.integrate.quad` applied to the parametrised integrand of (4) for the vortex on a circle returns Γ to 10⁻¹⁴ with about 21 function evaluations, because it uses an adaptive high-order rule; our midpoint rule needs n ≈ 1800 segments for 10⁻⁶ by (9). The midpoint rule is not the best integrator; it is the one whose error you can derive in two lines, which is what this lesson needs. In `P03` you will replace it. For the cylinder (3) on a 4 m square, the midpoint rule reaches a relative error of 4 × 10⁻⁴ at n = 64, 2.6 × 10⁻⁵ at n = 256 and 1.6 × 10⁻⁶ at n = 1024: the same second order, with a larger constant because the field varies faster near the body.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** No measurement is involved; we are mapping the behaviour of the integrals.

**Experiment A: convergence.** Prediction first: by (9) the circle error at n = 32 is 3.2 × 10⁻³ Γ and halves twice (factor 4) for every doubling of n; the square should show the same slope with a different constant. Run `circulation_convergence(field, gamma, [8, 16, 32, 64, 128, 256, 512])` for the point vortex and for the lifting cylinder, tabulate, plot on log-log axes with a reference n⁻² line, and fit the slope over n ≥ 32 with `np.polyfit(np.log(n), np.log(err), 1)`. Record the n at which each case first reaches 10⁻⁴ relative error.

**Experiment B: the Rankine sweep.** For the approach wake vortex (Γ = 300 m²/s, a = 1.5 m), compute the circulation around circles of radius 0.25, 0.5, …, 3.0 m with n = 2000 and tabulate against Γ min(1, (R/a)²). Prediction: a parabola inside the core, a flat line outside, and a corner at R = a. Then compare the convergence *constant* inside and outside: run n = 64, 128, 256, 512 on a circle of radius 1 m (inside the core, where the field is linear in r) and of radius 3 m (outside). Prediction: outside, (9) gives +π²/(3n²); inside, the midpoint rule is exact on each chord for a linear field, so the result is the exact circulation of the inscribed *polygon*, Γ(R/a)²(n/2π) sin(2π/n), a relative error of −2π²/(3n²): twice as large and of the opposite sign.

**Experiment C: the inlet profile.** Replace the uniform mass-flux field by a boundary-layer-like profile ρV(1 − r/R)^{1/7} and compute ṁ with `surface_flux` for ns = nt = 20, 40, 80, 160, 320. Prediction: the hand integral gives 49/60 of ρVπR² (problem 5), and because the profile has an infinite slope at the wall the midpoint rule should converge slower than second order. Measure the order.

Raw data format: three plain-text tables, (n, err_circle, err_square), (R, Γ_num, Γ_pred), (n, ṁ); keep them with your lesson notes.

## 5. Validate and improve

Put the sources side by side for the lifting cylinder, V∞ = 50 m/s, a = 0.5 m, Γ = −20 m²/s (clockwise, so the lift is upward):

| quantity | hand (§2) | code (n = 400) | difference |
|---|---|---|---|
| circulation, circle r = 1 m | −20 m²/s | −20.0004 | 2 × 10⁻⁵, equation (9) |
| circulation, circle r = 3 m | −20 m²/s | −20.0004 | same: (9) is independent of R |
| circulation, square 4 m × 4 m | −20 m²/s | −20.0002 | midpoint error, smaller constant |
| circulation, circle not enclosing | 0 | −0.0000 | rounding |
| L' = −ρV∞Γ | 1225 N/m | 1225.0 | 0 |

The lift is q c c_l for a 1 m chord at c_l = 0.8, exactly, because Γ = ½ c_l c V∞ was chosen to make it so; the code adds nothing to this line except confirming the sign convention: a clockwise vortex in a +x stream lifts in +y.

Every difference must be explained, and here each has a cause you can name:

- **The circle error is the same at r = 1 m and r = 3 m.** By (9) the vortex contribution does not depend on R, and the doublet contributes exactly zero to the discrete sum on any circle centred on the cylinder by symmetry (each segment's doublet contribution is cancelled by the segment diametrically opposite). If your two circles disagree, the centre of the path is not the centre of the cylinder.
- **The square is more accurate than the circle at the same n.** The square's sides are straight, so the only error is the midpoint sampling of a smoothly varying integrand, while the circle also replaces arcs by chords. Different constants, same order.
- **The Rankine sweep at R = a.** With n = 2000 the circulation at R = 1.5 m is 299.89 against a prediction of 299.89; the kink is on the path but the integrand (v_θ at radius R cos(π/n) ≈ 0.9999988 a) is sampled a hair inside the core, where the solid-body formula still applies, so nothing breaks. Inside the core the measured relative error is −1.606 × 10⁻³ at n = 64 and −2.510 × 10⁻⁵ at n = 512, against +8.04 × 10⁻⁴ and +1.255 × 10⁻⁵ outside: both second order, constants −2π²/3 and +π²/3 as predicted. The sign tells you which branch of the field the path is sampling.
- **The inlet profile converges at order 1.14, not 2.** The integrand ρV(1 − r/R)^{1/7} r has a derivative that diverges at the wall; the midpoint rule's error on the last cell is O(h^{8/7}) instead of O(h³), and that one cell sets the order. Measured: 266.03, 265.58, 265.38, 265.29, 265.26 kg/s for n = 20 to 320 against the exact 265.22 kg/s; the successive differences shrink by a factor of 2.2 rather than 4. Fix: substitute r = R(1 − w⁷), which makes the integrand smooth in w, or use a Gauss–Jacobi rule (`P03`).

What is still unverified when everything agrees: the model matches itself. L' = ρV∞Γ with Γ from a potential-flow field says nothing about whether a real airfoil develops that Γ; that needs the Kutta condition and, ultimately, a wind-tunnel balance (`P10`). What the lesson does establish is that *if* you have a velocity field, measured or computed, Γ is a well-defined number that does not depend on the loop you draw outside the core.

Improvements worth making: replace the midpoint rule by Simpson's rule on each segment (third field sample per segment, error ∝ h⁴); let `circulation_convergence` return the fitted slope alongside the errors. Improvements worth not making: adaptive refinement near the vortex centre for paths that pass through it; such a path has no circulation to compute.

## 6. Defend your work

Submit: `main.py` with all tests passing; `circulation_vortex.png`; the three tables from part 4 with the fitted convergence slopes; your hand derivation of (5), of the non-enclosing case in §2.3 and of (9); and a one-page note (`projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) stating what the potential-flow model assumes and what the inlet calculation assumes.

### Adversarial self-checks

- **What if Γ doubles?** Every circulation doubles and the lift doubles; the convergence error (9) doubles too, so the *relative* error is unchanged. A test asserts linearity in Γ. What if V∞ → 0? The lift goes to zero while Γ stays: circulation is a property of the field, lift is a property of the field *and* the stream.
- **What if the loop shrinks onto the vortex centre?** The circulation stays Γ for a point vortex at any radius, however small (equation (9) has no R in it), but for the Rankine vortex it falls as R². A loop of radius 0.1 a sees only 1 % of Γ. That is how you tell a point vortex from a real one.
- **Which assumption fails first in the real world?** Inviscid flow: the core of a real wake vortex is viscous and diffuses, with a growing as √(νt), and the circulation outside it decays slowly by viscous and turbulent transport; the potential-flow circulation is constant forever.
- **What measurement could falsify the model?** A particle-image-velocimetry map of the velocity around a wing section gives **v**(x, y) on a grid; computing the circulation on several loops and comparing with the balance-measured lift through L' = ρV∞Γ tests Kutta–Joukowski directly. Disagreement of more than a few percent outside the boundary layer and wake means the flow is not two-dimensional or not steady.
- **Does the result survive a change in path, n or surface parametrisation?** It must: three paths to 10⁻⁵ in the tests; n only changes the error by (9); the sphere flux is the same for a sphere centred off the source as long as it still encloses it. If your sphere flux changes sign when you swap s and t in the parametrisation, that is the cross-product order, not the physics.
- **Which decision changes if the uncertainty doubled?** Wake-vortex separation standards are set from Γ with margins; a 10 % uncertainty in Γ (from weight and speed) is a 10 % uncertainty in the peak swirl a following aircraft meets. Double it and the separation would have to grow. For the inlet, ṁ uncertain by 2 % instead of 1 % moves the engine's operating point by a comparable amount; the profile correction (18 %) dwarfs both, which is why it must be modelled, not budgeted.

## Check yourself

1. What are the SI units of circulation, of 2D source strength and of 3D source strength?
2. A loop encloses a point vortex of circulation Γ and winds around it twice. What is the line integral of velocity around it?
3. The vortex field is ∇(Γθ/2π) away from the origin. Why is its circulation around an enclosing loop not zero, as the fundamental theorem for line integrals seems to say?
4. For the Rankine vortex with core radius a, what fraction of Γ does a circle of radius a/2 see?
5. The midpoint-rule circulation on a circle with n = 100 segments has what relative error for a point vortex? Does the answer depend on the radius?
6. Why does `surface_flux` with the sphere parametrisation give a positive flux for the position field **r** = (x, y, z), and what would make it negative?
7. The mass flow through the 2 m fan face at 10 km and 250 m/s is 325 kg/s with a uniform profile. Roughly what is it with the 1/7-power profile, and why is it smaller?
8. Gravity's work from a 400 km orbit to an 800 km orbit is −3.28 MJ/kg along a straight line. What is it along a spiral that goes round the Earth three times on the way?
