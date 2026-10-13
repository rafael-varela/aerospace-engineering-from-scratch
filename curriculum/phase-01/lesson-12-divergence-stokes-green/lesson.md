# P01.12  Divergence, Stokes & Green's Theorems

> What happens inside a volume is written on its surface: divergence is flux per unit volume, curl is circulation per unit area, and the three big theorems are just that sentence made exact.

| | |
|---|---|
| **Learning outcomes** | 1. Compute divergence and curl of a vector field by hand, with SymPy, and by central differences, and state what each measures physically (dilatation rate, vorticity). 2. Derive Green's theorem for a rectangle from the fundamental theorem of calculus and derive the continuity equation from the divergence theorem. 3. Implement and test the two sides of the divergence theorem on a box control volume, of Green's theorem on a rectangle and of Stokes' theorem on a hemisphere, and measure their second-order convergence. 4. Explain why the discrete flux and volume integrals balance *exactly* for some fields and only to O(h²) for others. 5. Decide, from a flux balance on a control volume, whether a computed flow field conserves mass. |
| **Prerequisites** | `P01.10` (double and triple integrals), `P01.11` (line and surface integrals, circulation, flux), `P01.08` (partial derivatives), `P00.09` (SymPy). Skill check: you can write the outward normal of each face of a box and parametrise a hemisphere. |
| **Engineering question** | Air enters a control volume around the core of a wind-tunnel contraction through its sides and leaves through its downstream face: does the volume flow in equal the volume flow out, and what single pointwise property of the velocity field guarantees it? |
| **Assumptions and limits** | Fields are smooth inside the region (the Rankine vortex's vorticity step is the deliberate exception). The flow examples are incompressible and inviscid; the divergence theorem itself has no such restriction, but the "flux in = flux out" conclusion needs steady incompressible flow. Boxes and hemispheres are the only shapes integrated numerically; the midpoint rule is second order only for smooth integrands. The uniform-sphere gravity model holds inside a body of constant density only. |
| **Analytical work** | Divergence and curl of the contraction, rotation and expansion fields; Green's theorem on a rectangle from the fundamental theorem; the divergence theorem's physical meaning via the continuity equation; Stokes' theorem applied to a rigid rotation; the discrete telescoping identity that makes the box flux exact for quadratic fields; Gauss's law for gravity inside a uniform sphere. |
| **Computational deliverable** | `starter/main.py`: `gradient_sym`, `divergence_sym`, `curl_sym`, `divergence_fd`, `curl_fd`, example fields (`contraction_flow`, `rigid_rotation`, `uniform_expansion`, `rankine_vortex_2d`, `rankine_vorticity`, `uniform_sphere_gravity`), `box_face_fluxes`, `box_flux`, `box_volume_integral`, `divergence_theorem_check`, `polygon_circulation`, `rectangle_area_integral`, `greens_theorem_check`, `surface_flux`, `rim_circulation`, `stokes_theorem_check`, `divergence_theorem_convergence`. Tests prove: curl∘grad = 0 and div∘curl = 0 symbolically; finite differences match SymPy at second order; the contraction flow has zero net flux face by face; the divergence-theorem residual converges at second order and is exactly zero for quadratic fields; Green's theorem holds inside, outside and across the Rankine core; the flux of a curl through a closed box is zero; Stokes' theorem gives the same number on the rim, the cap and the disc. |
| **Practical activity** | **computational experiment**: a convergence study of the divergence-theorem residual, a study of what happens to Green's theorem when the vorticity has a step, and a test of which fields make the discrete theorem exact. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `divergence_theorem.png`, the three convergence tables with fitted orders, and a short note explaining the exact-for-quadratics observation. |
| **Assessment** | `assessments/phase-01/grading-rubric.md`; competencies 1–4. Pass: Green's theorem derived for the rectangle and the continuity equation derived from the divergence theorem, with assumptions named; all public tests pass; the convergence orders are measured and the two anomalies (exactness, the step) are explained; the judgment problem is answered with a measurable sign. |
| **References** | `@openstax-calculus-volume-3` 6.4 Green's Theorem (required-free); `@openstax-calculus-volume-3` 6.5 Divergence and Curl (required-free); `@openstax-calculus-volume-3` 6.7 Stokes' Theorem (required-free); `@openstax-calculus-volume-3` 6.8 The Divergence Theorem (required-free). |
| **Extensions** | ★ Replace the box by a tetrahedron and show the divergence theorem still balances (this is what a finite-volume CFD code does on every cell); ★ compute the volume and centroid of an ellipsoid from surface fluxes alone (problem 7) and then of a triangulated wing surface; ★ derive the vorticity transport equation and show that a 2D inviscid flow conserves the circulation of every material loop (Kelvin's theorem). |
| **Time** | unmeasured |

## 1. The engineering problem

Every flow solver, every propulsion calculation and every structural load balance rests on one bookkeeping idea: draw a closed surface, call the inside a **control volume**, and account for what crosses the surface. Mass that enters must leave or accumulate; momentum that enters must leave or push on the walls. In `P01.11` you computed the flux of a field through a surface. This lesson supplies the theorem that turns a surface flux into a volume integral and back, so that a *pointwise* property of the field (its divergence) tells you a *global* fact about the control volume (net flux), and the mirror theorem that turns the circulation around a loop into the integral of a pointwise property (curl, the vorticity) over any surface the loop bounds.

Take the core of a wind-tunnel contraction: air converges toward the axis and accelerates downstream. The simplest model of that is the axisymmetric field **v** = (−kx/2, −ky/2, kz), with k = 10 s⁻¹. Draw a box 0.4 m × 0.4 m × 1 m around the axis between z = 1 m and z = 2 m. Air comes in through the four side faces and the upstream face and goes out through the downstream face. The engineering question is: does the volume flow balance, which pointwise property guarantees that it does, and how accurately does a numerical flux balance confirm it?

Success criteria, written before we start: the six face fluxes of the contraction field, computed by the midpoint rule, sum to zero to 10⁻¹² of the largest face flux; the two sides of the divergence theorem for a smooth field with non-constant divergence agree with a residual that falls by a factor of four per halving of the cell size; Green's theorem on the Rankine wake vortex gives the same circulation from the loop and from the vorticity inside it, to 10⁻¹⁰ inside the core and to 1 % across the core edge with the reason for the 1 % stated; Stokes' theorem gives the same number for the rim circulation, the flux of curl through a hemisphere and the flux through the flat disc, to 10⁻³ at 40 cells and to 10⁻⁴ at 80.

## 2. Understand and derive

### 2.1 Divergence: flux per unit volume

For a field **F** = (F_x, F_y, F_z), the **divergence** is

∇ · **F** = ∂F_x/∂x + ∂F_y/∂y + ∂F_z/∂z.   (1)

(`@openstax-calculus-volume-3` 6.5.) Its meaning comes from a small box of sides Δx, Δy, Δz centred at a point. The outward flux through the two x faces is [F_x(x + Δx/2) − F_x(x − Δx/2)] Δy Δz ≈ (∂F_x/∂x) Δx Δy Δz, and likewise for y and z; the total outward flux divided by the volume Δx Δy Δz is (1). So divergence is **outward flux per unit volume**: positive where the field "spreads out", zero where what comes in goes out. For a velocity field the dimensions are 1/s, and ∇ · **v** is the fractional rate of growth of the volume of a small fluid parcel, its **dilatation rate**.

Three fields you should be able to differentiate in your head:

| field | ∇ · **v** | ∇ × **v** | what it is |
|---|---|---|---|
| contraction (−kx/2, −ky/2, kz) | −k/2 − k/2 + k = 0 | (0, 0, 0) | incompressible, irrotational acceleration toward +z |
| rigid rotation (−ωy, ωx, 0) | 0 | (0, 0, 2ω) | solid-body spin about z |
| uniform expansion k(x, y, z) | 3k | (0, 0, 0) | every parcel's volume grows as e^{3kt} |

**Assumption 1 (incompressible flow):** density is constant along a parcel's path, so its volume cannot change, so ∇ · **v** = 0 everywhere. That single equation is why the contraction's volume flow balances; §2.3 makes it exact.

### 2.2 Curl: circulation per unit area

The **curl** is

∇ × **F** = (∂F_z/∂y − ∂F_y/∂z, ∂F_x/∂z − ∂F_z/∂x, ∂F_y/∂x − ∂F_x/∂y).   (2)

Take a small rectangle Δx × Δy in a plane z = const, traversed counter-clockwise. The circulation is F_x(bottom) Δx + F_y(right) Δy − F_x(top) Δx − F_y(left) Δy ≈ (∂F_y/∂x − ∂F_x/∂y) Δx Δy. So the z component of the curl is **circulation per unit area** of a loop in the xy plane, and the other components are the same for loops in the other planes. For a velocity field the curl is the **vorticity** ω, and the rigid-rotation row of the table shows the normalisation: a fluid spinning at angular velocity ω has vorticity 2ω. In `P01.11` the point vortex had circulation Γ on every enclosing loop and zero on every other; now we can say why: its vorticity is zero everywhere except at the centre, where all of Γ is concentrated. The Rankine vortex spreads it uniformly over the core: ω_z = Γ/(πa²) for r < a, and 0 outside.

Two identities follow from equality of mixed partials (`P01.08`) and are worth memorising because the tests check them: ∇ × (∇φ) = 0 for any scalar φ (gradient fields are irrotational, which is why the conservative fields of `P01.11` had zero circulation), and ∇ · (∇ × **F**) = 0 for any field **F** (curl fields are solenoidal: they have no sources).

### 2.3 The divergence theorem

For a region V with closed boundary surface S and outward normal **n** (`@openstax-calculus-volume-3` 6.8),

∯_S **F** · **n** dA = ∭_V ∇ · **F** dV.   (3)

Why it is true, in one paragraph: tile V with small boxes. By §2.1 the outward flux of each box equals its divergence times its volume. Add up all the boxes: every interior face is shared by two boxes whose outward normals are opposite, so the fluxes through interior faces cancel in pairs and only the outer surface survives on the left, while the right side becomes the volume integral. The theorem is the statement that this cancellation is exact in the limit.

**Worked example: the contraction control volume.** Box [−0.2, 0.2] × [−0.2, 0.2] × [1, 2] m, k = 10 s⁻¹. On the face x = +0.2, u = −kx/2 = −1 m/s, pointing inward (−x direction) while the outward normal is +x; the face area is 0.4 × 1 = 0.4 m², so the outward flux is −0.4 m³/s. Same on x = −0.2 (u = +1 m/s inward against normal −x) and on both y faces: four side faces at −0.4 m³/s each. On z = 1, w = 10 m/s against the outward normal −z: −1.6 m³/s. On z = 2, w = 20 m/s along +z: +3.2 m³/s. Sum: −1.6 − 1.6 + 3.2 = 0. The right side of (3) is ∭ 0 dV = 0. In words: 1.6 m³/s enters through the sides, 1.6 m³/s enters upstream, 3.2 m³/s leaves downstream. Dimensions: (m/s)(m²) = m³/s. ✓

**Worked example: uniform expansion.** **v** = k**r** with k = 2 s⁻¹ in the same box: ∇ · **v** = 3k = 6 s⁻¹, so the right side is 6 × 0.16 = 0.96 m³/s. The left side: each face contributes k × (coordinate) × (area); the four side faces give 2 × 0.2 × 0.4 × 4 = 0.64, and the two z faces give 2 × (2 − 1) × 0.16 = 0.32, total 0.96. ✓ Fluid leaves at 0.96 m³/s: the box's contents are not conserved; the field has sources everywhere.

**Why engineers care: the continuity equation.** Let ρ(**r**, t) be a fluid's density and **v** its velocity. Mass inside a fixed control volume V is ∭ ρ dV, and it changes only because mass crosses S: d/dt ∭_V ρ dV = −∯_S ρ**v** · **n** dA. Apply (3) to the right side and move the time derivative inside the (fixed) volume: ∭_V [∂ρ/∂t + ∇ · (ρ**v**)] dV = 0. Because this holds for *every* control volume, however small, the integrand must vanish:

∂ρ/∂t + ∇ · (ρ**v**) = 0.   (4)

Steady flow with ρ constant reduces (4) to ∇ · **v** = 0, which is Assumption 1. The divergence theorem is what lets you pass from a statement about control volumes to a differential equation and back, and (4) is the first of the conservation laws of Phase 05.

### 2.4 Green's theorem

In two dimensions, for a region D with counter-clockwise boundary C (`@openstax-calculus-volume-3` 6.4),

∮_C (P dx + Q dy) = ∬_D (∂Q/∂x − ∂P/∂y) dA.   (5)

The integrand on the right is the z component of the curl of (P, Q, 0), so (5) says *circulation around C equals the vorticity integrated over D*. There is also a flux form, ∮ **F** · **n** ds = ∬ ∇ · **F** dA, which is the divergence theorem in the plane.

**Derivation for a rectangle** [a, b] × [c, d]. Take the Q term. By the fundamental theorem of calculus in x, ∫_a^b ∂Q/∂x dx = Q(b, y) − Q(a, y) for each y; integrate over y: ∬ ∂Q/∂x dA = ∫_c^d Q(b, y) dy − ∫_c^d Q(a, y) dy. The first term is the line integral of Q dy up the right side (y increasing, counter-clockwise), the second is the line integral down the left side (dy negative). On the top and bottom dy = 0. So ∬ ∂Q/∂x dA = ∮ Q dy. The P term works the same way with the roles of x and y swapped and a sign from the orientation: ∬ ∂P/∂y dA = ∫_a^b [P(x, d) − P(x, c)] dx = −∮ P dx, because the bottom is traversed with x increasing and the top with x decreasing. Subtract: (5). A general region is a union of rectangles in the limit, with interior edges cancelling exactly as in §2.3.

**Worked example: a wake vortex.** The Rankine vortex with Γ = 300 m²/s and a = 1.5 m (`P01.11`). A 1 m square centred on the axis lies inside the core, where ω_z = Γ/(πa²) = 42.44 s⁻¹: the right side of (5) is 42.44 × 1 = 42.44 m²/s, and the circulation around the square must be the same. A 6 m square encloses the whole core: the right side is ω_z × πa² = Γ = 300 m²/s. A 1 m square at (2, 2)–(3, 3) m is outside the core, where ω_z = 0: circulation zero, although the velocity on that loop is 10 m/s. Circulation is not "how fast the fluid goes round the loop"; it is how much vorticity the loop encloses.

### 2.5 Stokes' theorem

For an oriented surface S with boundary curve C (`@openstax-calculus-volume-3` 6.7), with the orientation of C and the normal of S related by the right-hand rule,

∮_C **F** · d**r** = ∬_S (∇ × **F**) · **n** dA.   (6)

Green's theorem is the flat case. The surprising content is on the right: *any* surface spanning C gives the same number, a hemisphere or the flat disc it shares its rim with. The reason: two surfaces with the same rim together form a closed surface, and the flux of a curl through a closed surface is zero by (3) and ∇ · (∇ × **F**) = 0.

**Worked example: rigid rotation.** **v** = (−ωy, ωx, 0), ω = 0.5 rad/s, around the circle of radius R = 2 m in the plane z = 0. On the circle the speed is ωR, tangential, so the circulation is ωR × 2πR = 2πωR² = 12.566 m²/s. The curl is (0, 0, 2ω) = (0, 0, 1) s⁻¹; its flux through the flat disc is 2ω × πR² = 12.566 m²/s. Through the upper hemisphere of radius R, **n** = **r**/R and (∇ × **v**) · **n** = 2ω z/R, so the flux is 2ω/R ∬ z dA = 2ω/R × πR³ = 2πωR² again (∬_hemisphere z dA = πR³, a standard integral you should verify in spherical coordinates). ✓

In fluid mechanics (6) is the bridge between the circulation a wing section develops and the vorticity in its boundary layer and wake; with the vorticity transport equation it gives Kelvin's theorem, that the circulation of a material loop is conserved in an inviscid barotropic flow, which is why a starting vortex of −Γ must be left behind when a wing develops +Γ.

### 2.6 Where the model stops being valid

- **Smoothness.** Both theorems assume the field is continuously differentiable inside the region. The Rankine vortex has a *step* in vorticity at r = a; the theorems still hold (the step is integrable), but the midpoint rule's accuracy for the area integral drops from second to first order across the step. Problem 4 measures it.
- **Singularities.** The point vortex has infinite vorticity at one point; Green's theorem over a region containing it needs the vorticity to be treated as a delta function with integral Γ. Numerically, never put a cell centre on a singularity.
- **Orientation.** Reversing the normal flips the sign of every flux; reversing the loop flips every circulation. The right-hand rule linking C and **n** in (6) is a convention, and a mismatch is the most common Stokes error.
- **Control-volume conclusions.** "Flux in equals flux out" requires steady, incompressible flow. Unsteady: mass accumulates, ∂ρ/∂t ≠ 0 in (4). Compressible: ρ**v**, not **v**, is divergence-free, so volume flow does not balance even when mass flow does.
- **Uniform-sphere gravity.** ∇ · **g** = −4πGρ holds pointwise for any density distribution (Gauss's law), but the linear field **g** = −(4π/3)Gρ**r** is only correct inside a *uniform* sphere; the real Earth's density rises from 2700 kg/m³ at the surface to 13 000 kg/m³ in the core.

## 3. Build the computational model

Open `starter/main.py`. The shared symbols come from `coords()`; symbolic fields are lists of three expressions, numerical fields are callables as in `P01.11`.

### Step 1: symbolic operators

```python
X, Y, Z = sp.symbols("x y z", real=True)

def divergence_sym(field):
    return sp.diff(field[0], X) + sp.diff(field[1], Y) + sp.diff(field[2], Z)

def curl_sym(field):
    fx, fy, fz = field
    return [sp.diff(fz, Y) - sp.diff(fy, Z), sp.diff(fx, Z) - sp.diff(fz, X), sp.diff(fy, X) - sp.diff(fx, Y)]
```

Equations (1) and (2) literally. `gradient_sym` is the three partials. `lambdify_field` (given) turns a symbolic field into a NumPy callable and broadcasts constant components to the input shape, because `sp.lambdify` of a constant returns a Python scalar, not an array, and the integrators need arrays.

### Step 2: finite-difference operators

```python
def divergence_fd(field, x, y, z, h=1e-4):
    dudx = (field(x + h, y, z)[0] - field(x - h, y, z)[0]) / (2.0 * h)
    dvdy = (field(x, y + h, z)[1] - field(x, y - h, z)[1]) / (2.0 * h)
    dwdz = (field(x, y, z + h)[2] - field(x, y, z - h)[2]) / (2.0 * h)
    return dudx + dvdy + dwdz
```

Central differences in each coordinate (`P01.01`, `P00.09`): error ∝ h², and the test measures that order on a field whose derivatives are not accidentally exact. `curl_fd` needs all nine partials; evaluate the field at the six shifted points once each and assemble (2). Both exist so that you can take the divergence of a field you only have as *data* (a CFD result, a PIV map), where SymPy cannot help.

### Step 3: the box flux, face by face

```python
def box_face_fluxes(field, bounds, n):
    (x0, x1), (y0, y1), (z0, z1) = bounds
    mids = [lo + (hi - lo) * (np.arange(n) + 0.5) / n for lo, hi in bounds]
    cell = [(hi - lo) / n for lo, hi in bounds]
    out = {}
    Ym, Zm = np.meshgrid(mids[1], mids[2], indexing="ij")
    da = cell[1] * cell[2]
    out["x0"] = float(-np.sum(field(np.full_like(Ym, x0), Ym, Zm)[0]) * da)
    out["x1"] = float(np.sum(field(np.full_like(Ym, x1), Ym, Zm)[0]) * da)
    ...  # y and z faces likewise
    return out
```

On the face x = x₀ the outward normal is −**x̂**, so the flux is −∬ u dA; on x = x₁ it is +∬ u dA. Midpoint rule: cell centres at `(i + 0.5) / n`, one sample per cell, times the cell area. Returning the six fluxes separately, not just their sum, is deliberate: the sum being zero tells you the flow balances, the parts tell you *where* the flow goes, which is what an engineer reads off a control volume. `box_flux` sums the dictionary.

### Step 4: the volume integral and the theorem

```python
def box_volume_integral(f, bounds, n):
    mids = [lo + (hi - lo) * (np.arange(n) + 0.5) / n for lo, hi in bounds]
    dv = np.prod([(hi - lo) / n for lo, hi in bounds])
    Xm, Ym, Zm = np.meshgrid(*mids, indexing="ij")
    return float(np.sum(f(Xm, Ym, Zm)) * dv)

def divergence_theorem_check(field, div, bounds, n):
    flux = box_flux(field, bounds, n)
    vol = box_volume_integral(div, bounds, n)
    return flux, vol, abs(flux - vol)
```

`div` is passed in separately so you can hand it the *exact* divergence from SymPy and test the theorem, not your differencing. For a field with non-constant divergence the residual should fall as n⁻².

### Step 5: a telescoping surprise

Run `divergence_theorem_check` on **F** = (x²y, yz², z sin x) with n = 2, 3, 4. The residual is zero to rounding at every n, even though neither side is the exact integral. Reason: in the volume sum, the x-derivative of F_x = x²y is 2xy, which is linear in x, and the midpoint rule integrates a linear function exactly over each column of cells: Σ_i 2x_i y_j Δx = (x₁² − x₀²) y_j. That is exactly what the two x faces contribute at the same (y_j, z_k) sample. The discrete flux and the discrete volume integral *telescope* into each other whenever each component depends on its own coordinate at most quadratically, so the discrete theorem is an identity, not an approximation. Finite-volume CFD codes are built on this: they define the discrete divergence *as* the face-flux sum, so mass is conserved to rounding on every cell whatever the accuracy of the fluxes themselves. The test `test_discrete_divergence_theorem_is_exact_...` pins this down, and `test_divergence_theorem_residual_converges_at_second_order` uses (y eˣ, z sin y, x cos z), where no component is quadratic in its own coordinate, to see a real residual.

### Step 6: Green's theorem

`polygon_circulation` is the midpoint-rule line integral from `P01.11`; `rectangle_area_integral` is the 2D midpoint rule; `greens_theorem_check` returns both sides and their difference. For the rigid rotation (−y, x) the circulation is exactly twice the area for *any* n, by the same telescoping argument in two dimensions.

### Step 7: Stokes' theorem on a hemisphere

```python
def stokes_theorem_check(field, curl, radius, n):
    circ = rim_circulation(field, radius, 4 * n)
    cap = surface_flux(curl, hemisphere_surface(radius), (0.0, np.pi / 2), (0.0, 2 * np.pi), n, n)
    disc = surface_flux(curl, disc_surface(radius), (0.0, radius), (0.0, 2 * np.pi), n, n)
    return circ, cap, disc
```

`surface_flux` is the parametric-surface integrator of `P01.11` with the normal **r**_s × **r**_t by central differences. The hemisphere's parametrisation puts +z at the pole, so its normal is outward, matching the counter-clockwise rim by the right-hand rule. `rim_circulation` is a 3D line integral around the circle in z = 0.

### Verify

Run `python -m pytest curriculum/phase-01/lesson-12-divergence-stokes-green` from the repo root. In order of how much they would hurt if wrong:

1. **Symbolic identities.** div and curl of (y eˣ, z sin y, x cos z) against hand results; ∇ × ∇φ = 0; ∇ · (∇ × **F**) = 0; ∇²(r²) = 6. A wrong sign in (2) fails the first and the third.
2. **Finite differences.** Second-order convergence against SymPy for both operators on fields chosen so that no term is accidentally exact; the table in §2.1 (contraction: div 0, curl 0; rotation: curl 2ω; expansion: div 3k; uniform-sphere gravity: div −4πGρ).
3. **The control volume.** Each face flux of the contraction has the sign it should (in through five faces, out through one), the top face carries kz × area, and the sum is zero to 10⁻¹².
4. **The divergence theorem.** Exact for the uniform expansion (linear field), second order for the smooth field, exact for the quadratic field, zero for the flux of a curl through the closed box, linear in the field.
5. **Green's theorem.** Rigid rotation: 2 × area exactly; a general smooth field: second order; the Rankine vortex: exact inside the core (uniform vorticity, linear field), zero outside, Γ across the core edge with the area integral only first-order accurate.
6. **Stokes' theorem.** The hemisphere normal points outward (flux of **r** is 2πR³); rigid rotation gives 2πωR² on rim, cap and disc; a general field agrees on all three to 10⁻⁴ with the cap converging at roughly second order; the rim circulation of a gradient field is zero.

### Compare

SymPy ships these operators: `sympy.vector` has `divergence` and `curl` for fields built on a `CoordSys3D`. Compare: `from sympy.vector import CoordSys3D, divergence, curl; N = CoordSys3D("N"); divergence(N.y*sp.exp(N.x)*N.i + ...)` gives the same expressions as `divergence_sym` with the symbols renamed. The library handles curvilinear coordinate systems (cylindrical, spherical) correctly, which our Cartesian-only version does not; that is the reason to use it once you need ∇ · in cylindrical coordinates for an axisymmetric nozzle. For the integrals there is no library function to compare with; `scipy.integrate.tplquad` on the volume integral of the smooth field gives 0.2384032, and `box_volume_integral` at n = 64 gives 0.2384036: a difference of 4 × 10⁻⁷, which is the midpoint rule's own error at that n (it falls by four when n doubles). That is the right check to run once (problem 4).

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** No measurement is involved; we are studying the behaviour of the discrete theorems.

**Experiment A: convergence of the divergence-theorem residual.** Prediction first: both the face fluxes and the volume integral are midpoint sums with error ∝ h² = (L/n)², so the residual should fall by a factor of 4 per doubling of n. Run `divergence_theorem_convergence` for **F** = (y eˣ, z sin y, x cos z) on the contraction box with n = 2, 4, 8, 16, 32, 64; tabulate; fit the slope with `np.polyfit(np.log(n), np.log(res), 1)`. Then swap in **F** = (x²y, yz², z sin x) and record the residuals: they should sit at rounding (10⁻¹⁶). Then try (x³, 0, 0), cubic in its own coordinate, and see the residual come back.

**Experiment B: a step in the vorticity.** Green's theorem for the Rankine vortex on the 6 m square with n = 25, 50, 100, 200, 400, 800. Prediction: the circulation (a smooth integrand along the loop, which lies outside the core) converges at second order; the area integral of the stepped vorticity converges only at first order, because the cells straddling the circle r = a are counted either fully in or fully out. Tabulate both errors against Γ = 300 and fit both slopes; expect the area slope to be noisy.

**Experiment C: which surface?** Stokes' theorem for **F** = (y eˣ, z sin y, x cos z) with R = 1.3 m and n = 10, 20, 40, 80: tabulate rim, cap and disc. Prediction: all three converge to the same number at second order; the disc has the smallest error constant because its integrand is smooth in s all the way to the rim, while the cap's area element R² sin s and the rim's chord-versus-arc geometry each add a larger constant.

Raw data format: three plain-text tables (n, residual), (n, err_circ, err_area), (n, rim, cap, disc), kept with your lesson notes.

## 5. Validate and improve

Put the sources side by side for the contraction control volume and the two theorems:

| quantity | hand (§2) | code | difference |
|---|---|---|---|
| side-face fluxes, contraction | −0.4 m³/s each | −0.40000 | 0 (linear field, midpoint exact) |
| z = 1 and z = 2 faces | −1.6, +3.2 m³/s | −1.60000, +3.20000 | 0 |
| net flux | 0 | 0.0 × 10⁰ | rounding |
| expansion, flux and volume integral | 0.96 m³/s | 0.960000, 0.960000 | 1 × 10⁻¹⁶ |
| Green, 1 m square inside the core | 42.44 m²/s | 42.4413, 42.4413 | 7 × 10⁻¹⁵ |
| Green, 6 m square around the core, n = 200 | 300 m²/s | 300.0008, 300.2299 | 8 × 10⁻⁴ (loop), 0.23 (area) |
| Stokes, rigid rotation, n = 40 | 12.5664 m²/s | 12.5631, 12.5696, 12.5664 | 3 × 10⁻³, 3 × 10⁻³, 0 |

Every difference must be explained, and here each has a cause you can name:

- **The exact rows.** The contraction, expansion and rigid-rotation fields are linear, and the midpoint rule integrates linear functions exactly, so the only error is rounding. If your contraction net flux is 10⁻³ rather than 10⁻¹⁶, a face sign is wrong, not the resolution.
- **The area integral across the core edge is off by 0.08 %.** With n = 200 the cells are 3 cm; the circle r = 1.5 m cuts through about 2πa/0.03 ≈ 300 cells, each counted fully inside or fully outside the core, with errors of either sign that only partly cancel. The errors for n = 50, 100, 200, 400, 800 are 4.2, 1.9, 0.23, 0.12, 0.03 m²/s: falling roughly as 1/n on average, but erratically, because the result depends on exactly how the grid happens to straddle the circle. The loop integral, which never touches the step, is at 8 × 10⁻⁴ and falls cleanly by four per doubling: the second-order midpoint error of `P01.11`.
- **The rim circulation of the rotation is 3 × 10⁻⁴ low at 160 segments** by exactly the (n/2π) sin(2π/n) polygon-area factor of `P01.11` (the field is linear, so the polygon circulation is exact for the *polygon*, which has 1 − 2π²/(3 · 160²) = 0.99974 of the circle's area). The hemisphere flux is 2.6 × 10⁻⁴ high because the midpoint rule in the polar angle s over-weights the sin s factor near the rim. The disc is exact: its integrand 2ω s is linear in s.
- **The sum of residuals across tests is not an uncertainty.** Each residual is a deterministic discretisation error with a known order; the numbers above are not noisy and will reproduce to every digit on another machine.

What is still unverified when everything agrees: the theorems have been checked on boxes, rectangles and hemispheres. A finite-volume solver applies them on skewed hexahedra and tetrahedra, where the face normals are not axis-aligned and the midpoint of a face is not its centroid; the telescoping identity still holds (it depends only on each interior face being shared with opposite normals), but the second-order accuracy does not survive badly skewed cells. That is a Phase 10 topic.

Improvements worth making: let `box_face_fluxes` accept a different n per direction; add Simpson weights to `rectangle_area_integral` and watch the Rankine area integral *not* improve (the step defeats any polynomial rule). Improvements worth not making: refining the grid near the core edge by hand; the right tool is to integrate the vorticity in polar coordinates where the step is a grid line.

## 6. Defend your work

Submit: `main.py` with all tests passing; `divergence_theorem.png`; the three tables of part 4 with fitted slopes; your derivations of Green's theorem for the rectangle and of the continuity equation; and a one-page note (`projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) stating what "the flow balances" assumes and explaining the exactness observation in Step 5.

### Adversarial self-checks

- **What if k doubles?** Every contraction face flux doubles and the net stays zero; the expansion flux doubles. Linearity in the field is a test. What if the box is moved off the axis? The side fluxes are no longer equal and opposite in pairs, but their sum with the z faces is still zero: the balance does not depend on symmetry, only on ∇ · **v** = 0.
- **What if the box shrinks to a point?** Flux/volume → ∇ · **v** at that point: that is the definition. Numerically, at fixed n, halving every side of the box halves h (residual ÷ 4) and divides the volume by 8, so the residual of the smooth field should fall by about 2⁵ = 32: measured 2.48 × 10⁻⁵ → 6.50 × 10⁻⁷ → 1.83 × 10⁻⁸ at n = 8, factors of 38 and 36 (the fifth-power scaling plus a slowly varying derivative).
- **Which assumption fails first in the real world?** Steadiness: in a wind tunnel starting up, ∂ρ/∂t ≠ 0 and the flux balance has a storage term; next, incompressibility above about Mach 0.3, where ∇ · **v** ≠ 0 even in steady flow while ∇ · (ρ**v**) = 0 still holds.
- **What measurement could falsify the model?** A hot-wire or pitot traverse of the six faces of a real control volume in a tunnel: the measured volume flows must balance to within the probe uncertainty (typically 1–2 %). A persistent imbalance means a leak, an unsteady flow or a probe calibration error, and the face-by-face breakdown says which face to suspect.
- **Does the result survive a change in n, in the surface, in the orientation?** The theorems must, and the tests check two surfaces for Stokes and six values of n for Gauss; the *sign* must flip if you reverse the normal or the loop, and if it does not, your integrator is ignoring orientation.
- **Which decision changes if the uncertainty doubled?** A 2 % instead of 1 % flux imbalance in a tunnel calibration moves the inferred test-section speed by 1 % and dynamic pressure by 2 %, which is the whole error budget of a drag measurement; the balance check is what decides whether a run is usable.

## Check yourself

1. What are the dimensions of the divergence and of the curl of a velocity field, and what does each measure?
2. A fluid rotates as a solid body at 3 rad/s. What is its vorticity?
3. State the divergence theorem in words, naming which side is a surface integral and which a volume integral, and say which way the normal points.
4. Why is the flux of ∇ × **F** through any closed surface zero? Give both a one-line symbolic reason and a one-line Stokes reason.
5. The circulation around a 1 m square outside the Rankine core is zero although the fluid moves at 10 m/s along the loop. Explain in one sentence.
6. Which of these fields make the discrete divergence theorem on a box exact at every n: (x, y, z); (x², y², z²); (x³, 0, 0); (y eˣ, 0, 0)?
7. Steady compressible flow through a nozzle: which field is divergence-free, **v** or ρ**v**, and which quantity is therefore the same at every cross-section?
8. Stokes' theorem gives the same number for the hemisphere and the disc with the same rim. What single identity guarantees it, and which earlier theorem do you combine it with?
