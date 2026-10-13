# P01.10  Multiple Integrals

> Mass, centre of mass, inertia tensor: three integrals over the same lump of metal, and every flight-dynamics equation you will ever write needs all three.

| | |
|---|---|
| **Learning outcomes** | 1. Derive, as iterated integrals, the mass, centre of mass and inertia tensor of a tapered swept wing box; reduce them to one-dimensional closed forms by integrating over the cross-section first; derive the Jacobian of the map from the unit cube to the box and the parallel-axis theorem for the full tensor. 2. Implement and test a mass-properties module that evaluates the same integrals three ways (composite midpoint rule on the unit cube, SymPy term-by-term, `scipy.integrate.tplquad`) and agree. 3. Measure, by computational experiment, the convergence order of the midpoint rule for the mass and for $I_{xx}$, and explain why the second is not the clean second order of the first. 4. Decide what level of geometric detail the inertia tensor needs for a flight-dynamics model, from how the tensor's entries respond to taper and to the smeared-density assumption. |
| **Prerequisites** | `P01.03`–`P01.04` (integration, numerical quadrature and its order), `P01.07` (vectors, cross product, moments), `P01.08` (functions of several variables), `P00.03` (`np.meshgrid`, broadcasting), `P00.09` (SymPy). Skill check: you can integrate a polynomial in one variable, and you know what $\int_0^1\int_0^1 xy\,\mathrm{d}x\,\mathrm{d}y$ means. |
| **Engineering question** | What are the mass, centre of mass and inertia tensor of a tapered, swept wing box, how accurately does a simple grid-based integration recover them, and which of the geometric details (taper, sweep, depth) actually move the numbers that a flight-dynamics model needs? |
| **Assumptions and limits** | The wing box is a solid of uniform *effective* density $\rho_{\text{eff}}$: the skins, spars, ribs and fuel are smeared into one homogeneous block occupying the box between two chord fractions and of depth $\tau c$. Linear taper, straight leading-edge sweep, no twist, no dihedral; the chord is horizontal. The smeared density gets the mass and CG about right if $\rho_{\text{eff}}$ is calibrated to the real wing's mass, but it misplaces the inertia whenever the real mass sits at the skins rather than filling the box; the error is bounded in §2.9. The model says nothing about the stiffness, the loads or the structure's layout; it is the first input to Phase 09, not an answer about Phase 06. |
| **Analytical work** | Double integrals over rectangles and general regions as iterated integrals (Fubini); triple integrals for $m$, $\mathbf r_{cg}$ and $I_{ij} = \rho\int(r^2\delta_{ij} - r_ir_j)\,\mathrm{d}V$; the change of variables from the unit cube to the box with Jacobian $(f_1 - f_0)\tau Lc(\eta)^2$; closed forms $m = \rho(f_1 - f_0)\tau L(c_r^2 + c_rc_t + c_t^2)/3$ and $y_{cg} = L(c_r^2 + 2c_rc_t + 3c_t^2)/[4(c_r^2 + c_rc_t + c_t^2)]$; the polar double integral for a ring; the parallel-axis theorem $I_O = I_{cg} + m(|\mathbf d|^2E - \mathbf d\mathbf d^{\mathsf T})$; the midpoint rule's $h^2$ error and why it cancels for one integrand. |
| **Computational deliverable** | `starter/main.py`: `chord`, `map_to_box`, `midpoint_3d`, `mass_properties_midpoint`, `mass_properties_analytic`, `mass_closed_form`, `spanwise_cg_closed_form`, `parallel_axis`, `inertia_about_cg`, `ring_inertia_polar`, `ring_inertia_closed_form`, `convergence_study`, `tplquad_mass_and_ixx`. Tests prove: the map hits the box corners and its Jacobian is the analytic one; the midpoint rule integrates constants and the volume; mass by three methods agrees to $10^{-10}$; the midpoint mass error falls by 4 per doubling; the CG lies inside the box and on the symmetry plane; the inertia tensor is symmetric, positive definite and obeys the triangle inequalities; midpoint and `tplquad` tensors match SymPy; the parallel-axis shift round-trips and the CG minimises every diagonal entry; $I \propto \rho L^5$; the polar ring matches its closed form and the thin-ring limit $mr^2$. |
| **Practical activity** | **computational experiment**: convergence of the midpoint mass and $I_{xx}$ with $n$ from 2 to 128, Richardson extrapolation, and a sensitivity sweep of the tensor against taper ratio and sweep. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `wing_box_mass.png`, the hand calculation of $m$, $y_{cg}$ and $I_{xx}$ with its comparison to the code, the convergence table with both observed orders and an explanation of the sign change, and the taper-sweep sensitivity table with one paragraph. |
| **Assessment** | `assessments/phase-01/grading-rubric.md`; competencies 1–4. Pass: the Jacobian and both closed forms derived; all public tests pass; the hand values of $m$, $y_{cg}$ and $I_{xx}$ agree with the code to four figures; the convergence experiment reports the order of the mass and explains with a formula why $I_{xx}$ behaves differently. |
| **References** | `@openstax-calculus-volume-3` §5.1 "Double Integrals over Rectangular Regions" (required-free); §5.2 "Double Integrals over General Regions" (required-free); §5.3 "Double Integrals in Polar Coordinates" (required-free); §5.4 "Triple Integrals" (required-free); §5.6 "Calculating Centers of Mass and Moments of Inertia" (required-free); §5.7 "Change of Variables in Multiple Integrals" (required-free); §5.5 "Triple Integrals in Cylindrical and Spherical Coordinates" (supplementary-free, needed for Problem 3). |
| **Extensions** | ★ Replace the solid box by two skins of thickness $t$ at $z = \pm\tfrac12\tau c$ and two spar webs, with the same total mass, and compare the inertia tensors. ★ Add linear twist $\theta(y)$ to the map and recompute; show the $I_{xz}$ product of inertia is no longer zero. |
| **Time** | unmeasured |

## 1. The engineering problem

Every equation of motion for an aircraft or a spacecraft contains the mass, the position of the centre of mass and the inertia tensor. The mass sets the acceleration for a given force; the centre of mass is where that force appears to act and the point about which the vehicle rotates; the inertia tensor sets the angular acceleration for a given moment and, through its off-diagonal entries, couples roll into yaw. Phase 09 will write those equations; this lesson computes their inputs for one component, a wing box, from its geometry.

A wing box is a tapered beam between two spars: at each spanwise station $y$ it occupies the chord fractions $f_0c(y)$ to $f_1c(y)$ and the depth $\pm\tfrac12\tau c(y)$, and its leading edge is swept back by an angle $\Lambda$. We smear its skins, spars, ribs and fuel into a solid of uniform effective density. The mass is then the volume integral of the density, the centre of mass is the mass-weighted mean position, and each inertia entry is a mass-weighted mean of a quadratic in the coordinates. All three are **multiple integrals** over a region whose boundaries depend on the coordinates, which is the one thing that makes them harder than the single integrals of `P01.03`.

Three independent routes will be built: closed-form integration (SymPy, exact because the integrands are polynomials), a composite midpoint rule on a grid after a change of variables that straightens the box into a unit cube, and the library's adaptive `tplquad` in the original coordinates. Success criteria, written before we start: the three routes agree on the mass to $10^{-10}$ and on the inertia tensor to $10^{-6}$; the midpoint mass error falls by a factor of 4 per doubling of the grid; the hand calculation of $m$, $y_{cg}$ and $I_{xx}$ agrees with the code to four figures; and the parallel-axis theorem round-trips the tensor between the root and the centre of mass to machine precision.

## 2. Understand and derive

### 2.1 Double integrals and Fubini's theorem

The integral of $f(x, y)$ over a rectangle $R = [a, b] \times [c, d]$ is the limit of Riemann sums over a grid of cells, $\iint_R f\,\mathrm{d}A = \lim\sum f(x_i^*, y_j^*)\Delta x\Delta y$, and **Fubini's theorem** (`@openstax-calculus-volume-3` §5.1) says it equals either iterated integral:

$$\iint_R f\,\mathrm{d}A = \int_c^d\left[\int_a^b f(x, y)\,\mathrm{d}x\right]\mathrm{d}y = \int_a^b\left[\int_c^d f(x, y)\,\mathrm{d}y\right]\mathrm{d}x.$$   (1)

The inner integral is an ordinary one-variable integral with the outer variable held fixed, exactly the "hold everything else still" of `P01.08`. Over a **general region** (§5.2) bounded by curves, the inner limits depend on the outer variable: a region $a \le y \le b$, $g_1(y) \le x \le g_2(y)$ gives $\int_a^b\int_{g_1(y)}^{g_2(y)}f\,\mathrm{d}x\,\mathrm{d}y$. The order matters for the *limits*, not for the value, and the practical rule is: integrate first over the direction in which the region is described by two simple bounding functions.

### 2.2 Triple integrals for mass properties

For a solid $E$ of density $\rho(\mathbf r)$ (§5.4, §5.6):

$$m = \iiint_E\rho\,\mathrm{d}V, \qquad \mathbf r_{cg} = \frac1m\iiint_E\rho\,\mathbf r\,\mathrm{d}V, \qquad I_{ij} = \iiint_E\rho\,(|\mathbf r|^2\delta_{ij} - r_ir_j)\,\mathrm{d}V.$$   (2)

The diagonal entries are the familiar moments of inertia, $I_{xx} = \int\rho(y^2 + z^2)\,\mathrm{d}V$: the mass-weighted mean of the squared distance from the $x$-axis. The off-diagonal **products of inertia** $I_{xy} = -\int\rho xy\,\mathrm{d}V$ are zero whenever the body is symmetric about a coordinate plane, and non-zero for a swept wing about its root, which is why a swept aircraft's roll and yaw are coupled. The tensor is symmetric and positive definite, and its diagonal entries obey $I_{xx} + I_{yy} \ge I_{zz}$ (and permutations), because $I_{xx} + I_{yy} - I_{zz} = 2\int\rho z^2\,\mathrm{d}V \ge 0$. Those two facts are tests.

**Geometry.** Axes at the root leading edge: $x$ aft along the chord, $y$ spanwise, $z$ up. Chord $c(y) = c_r + (c_t - c_r)y/L$; leading edge $x_{le}(y) = y\tan\Lambda$. The box is

$$E = \{0 \le y \le L,\ x_{le}(y) + f_0c(y) \le x \le x_{le}(y) + f_1c(y),\ -\tfrac12\tau c(y) \le z \le \tfrac12\tau c(y)\}.$$

Our numbers: $L = 5$ m, $c_r = 1.6$ m, $c_t = 0.8$ m, $f_0 = 0.15$, $f_1 = 0.65$, $\tau = 0.10$, $\Lambda = 10°$, $\rho = 180$ kg/m³ (a solid box of that density has the mass of a real light-aircraft wing of this size with fuel).

### 2.3 Integrate over the cross-section first: closed forms

The box's cross-section at station $y$ is a rectangle of width $(f_1 - f_0)c$ and height $\tau c$, so its area is $A(y) = (f_1 - f_0)\tau c(y)^2$. Integrating over $x$ and $z$ first turns every triple integral into a single integral in $y$:

$$m = \rho(f_1 - f_0)\tau\int_0^Lc(y)^2\,\mathrm{d}y = \rho(f_1 - f_0)\tau L\,\frac{c_r^2 + c_rc_t + c_t^2}{3}.$$   (3)

(The last step: $c$ is linear in $y$, so $\int_0^Lc^2\,\mathrm{d}y = L(c_r^2 + c_rc_t + c_t^2)/3$, the same algebra as the volume of a frustum.) Hand values: $c_r^2 + c_rc_t + c_t^2 = 2.56 + 1.28 + 0.64 = 4.48$ m²; $\int c^2 = 5 \times 4.48/3 = 7.4667$ m³; volume $0.5 \times 0.1 \times 7.4667 = 0.37333$ m³; $m = 180 \times 0.37333 = 67.20$ kg. Similarly

$$y_{cg} = \frac{\int_0^Lyc^2\,\mathrm{d}y}{\int_0^Lc^2\,\mathrm{d}y} = L\,\frac{c_r^2 + 2c_rc_t + 3c_t^2}{4(c_r^2 + c_rc_t + c_t^2)} = 5 \times \frac{2.56 + 2.56 + 1.92}{4 \times 4.48} = 1.9643\ \mathrm{m},$$   (4)

which is $0.393L$: inboard of mid-span because the root is wider. For an untapered box (4) gives $L/2$; for a triangle ($c_t = 0$) it gives $L/4$ ✓. The chordwise CG needs $\int c^3$: $x_{cg} = [\tan\Lambda\int yc^2 + \tfrac12(f_0 + f_1)\int c^3]/\int c^2 = [0.17633 \times 14.667 + 0.4 \times 9.600]/7.4667 = 0.8606$ m. And $I_{xx}$ splits into a spanwise part and a depth part:

$$I_{xx} = \rho(f_1 - f_0)\tau\left[\int_0^Ly^2c^2\,\mathrm{d}y + \frac{\tau^2}{12}\int_0^Lc^4\,\mathrm{d}y\right] = 9\,[42.667 + 0.000833 \times 12.698] = 384.0 + 0.095 = 384.10\ \mathrm{kg\,m^2}.$$   (5)

The depth part is 0.025 % of the total: for the rolling inertia of a wing, the box's thickness is irrelevant and only the spanwise mass distribution matters. That is the first design conclusion of the lesson, and it came from a hand calculation.

### 2.4 Change of variables: straighten the box into a cube

A grid-based rule wants a rectangular domain. Define the map from the unit cube $(\xi, \eta, \zeta) \in [0, 1]^3$ (§5.7):

$$y = \eta L, \qquad x = y\tan\Lambda + [f_0 + (f_1 - f_0)\xi]\,c(y), \qquad z = (\zeta - \tfrac12)\tau c(y).$$   (6)

The **Jacobian determinant** is the volume-scaling factor of the map:

$$J = \det\frac{\partial(x, y, z)}{\partial(\xi, \eta, \zeta)} = \det\begin{pmatrix}(f_1 - f_0)c & \partial x/\partial\eta & 0\\ 0 & L & 0\\ 0 & \partial z/\partial\eta & \tau c\end{pmatrix} = (f_1 - f_0)\tau L\,c(\eta)^2.$$   (7)

The two $\partial/\partial\eta$ entries (shear from sweep and taper) sit off the diagonal in a column with zeros elsewhere and do not enter the determinant: shearing a slab does not change its volume. Then $\iiint_Ef\,\mathrm{d}V = \int_0^1\int_0^1\int_0^1f(x(\xi, \eta, \zeta), \ldots)\,J\,\mathrm{d}\xi\,\mathrm{d}\eta\,\mathrm{d}\zeta$, and $J$ is just the cross-section area times $L$, as (3) already said. The map's four corners at $\eta = 0$ and $\eta = 1$ and its Jacobian against (7) are the first test.

### 2.5 The composite midpoint rule in three dimensions

Divide the cube into $n^3$ cells of side $h = 1/n$ and evaluate the integrand at each cell centre:

$$\int_{[0,1]^3}g\,\mathrm{d}\xi\,\mathrm{d}\eta\,\mathrm{d}\zeta \approx \frac{1}{n^3}\sum_{i,j,k}g(\xi_i, \eta_j, \zeta_k), \qquad \xi_i = \frac{i + \tfrac12}{n}.$$   (8)

The one-dimensional midpoint rule has error $\tfrac{h^2}{24}[g'(1) - g'(0)] + O(h^4)$ (`P01.04`); in three dimensions the leading error is the sum of that expression in each direction. For the mass the integrand is $g = \rho J = \rho(f_1 - f_0)\tau Lc(\eta)^2$, which depends on $\eta$ only, and $g'(1) - g'(0) = 2\rho(f_1 - f_0)\tau L(c_t - c_r)^2 \cdot(\ldots) \ne 0$, so the error is a clean $h^2$: the mass errors at $n = 2, 4, 8, 16, 32$ are $-0.600, -0.150, -0.0375, -0.00938, -0.00234$ kg, dividing by exactly 4 each time, and always negative because $c^2$ is convex and the midpoint rule underestimates a convex integrand. For $I_{xx}$ the integrand is $\rho(y^2 + z^2)J$, and §4 shows that the $y^2$ part's leading $h^2$ term happens to vanish for this taper ratio, so the two parts converge at different orders and the total error changes sign at $n \approx 16$. A rule that converges "at second order" does so for each term separately, and the observed order of a sum can be anything until the slowest term dominates.

### 2.6 Polar coordinates: a fuselage frame

A circular ring (a fuselage frame, $r_{in} \le r \le r_{out}$, length $\ell$ along its axis) is a general region in $(x, y)$ but a rectangle in polar coordinates $(r, \theta)$ (§5.3). The area element is $\mathrm{d}A = r\,\mathrm{d}r\,\mathrm{d}\theta$, the Jacobian of $(r, \theta) \mapsto (r\cos\theta, r\sin\theta)$:

$$m = \rho\ell\int_0^{2\pi}\int_{r_{in}}^{r_{out}}r\,\mathrm{d}r\,\mathrm{d}\theta = \rho\ell\pi(r_{out}^2 - r_{in}^2), \qquad I_{\text{axis}} = \rho\ell\int_0^{2\pi}\int_{r_{in}}^{r_{out}}r^2\cdot r\,\mathrm{d}r\,\mathrm{d}\theta = \tfrac12m(r_{in}^2 + r_{out}^2).$$   (9)

Forgetting the $r$ in the area element is the classic error and gives a mass proportional to $r_{out} - r_{in}$ instead of $r_{out}^2 - r_{in}^2$. An aluminium frame of radius 0.95–1.00 m and 50 mm length has $m = 2700 \times 0.05 \times \pi(1 - 0.9025) = 41.35$ kg and $I = \tfrac12 \times 41.35 \times 1.9025 = 39.34$ kg m², within 0.07 % of the thin-ring value $mr^2$ at the mean radius 0.975 m (the ring is thin: 5 % of its radius). The midpoint rule on the $(r, \theta)$ rectangle reproduces both to $10^{-6}$ with $n = 64$, because the integrands $r$ and $r^3$ are smooth and periodic in $\theta$.

### 2.7 The parallel-axis theorem for the tensor

With $\mathbf r = \mathbf r_{cg} + \mathbf r'$ in (2) and $\int\rho\mathbf r'\,\mathrm{d}V = \mathbf 0$ by definition of the CG, the cross terms vanish and

$$I_O = I_{cg} + m\left(|\mathbf d|^2E - \mathbf d\mathbf d^{\mathsf T}\right), \qquad \mathbf d = \mathbf r_{cg} - \mathbf r_O,$$   (10)

$E$ the identity. The added matrix is the inertia tensor of a point mass $m$ at $\mathbf d$; it is positive semi-definite, so every diagonal entry of $I$ is smallest about the CG, which is a test. Our box about the root leading edge has $I_{xx} = 384.10$; about its CG at $(0.861, 1.964, 0)$ the shift removes $m(y_{cg}^2 + z_{cg}^2) = 67.2 \times 3.858 = 259.3$ and leaves $124.81$ kg m². A uniform rod of the same mass and length would have $mL^2/12 = 140.0$; the tapered box is 11 % lower because its mass sits toward the root. The principal moments about the CG, the eigenvalues of $I_{cg}$, are $2.45, 126.4, 128.7$ kg m²: a slender body with one small axis (the span) and two nearly equal large ones.

### 2.8 Which route when

- **Closed form** (3)–(5): exact, instant, and the only route that shows *which* geometric parameter matters. Available only because the integrand is a polynomial over a region described by linear functions; a real wing with a curved planform or a cambered section loses it.
- **Midpoint on the mapped cube**: $n^3$ function evaluations, second-order, trivial to vectorise (`np.meshgrid`), and extends to any integrand (a density that varies with position, Problem 6) and to any region you can map to a cube. $n = 32$ ($32\,768$ points) gives the mass to $3 \times 10^{-5}$ relative in a few milliseconds.
- **`scipy.integrate.tplquad`**: adaptive Gauss–Kronrod in each of three nested directions; handles $y$-dependent limits directly, needs no map, returns an error estimate, and is the oracle the tests compare against. It is slow (hundreds of thousands of integrand calls for $10^{-10}$ relative) and the integrand must be a scalar function, not vectorised.

### 2.9 Where the model stops being valid

- **Smeared density.** A real wing box carries most of its structural mass in the skins at $z = \pm\tfrac12\tau c$. For $I_{xx}$ that is irrelevant (the depth part was 0.025 %); for $I_{yy}$ about the CG, which is dominated by chordwise and depth terms, concentrating the mass at the skins instead of filling the box changes the $z^2$ contribution by a factor of up to 3 ($\int z^2$ for two skins is $\tfrac14(\tau c)^2$ per unit mass against $\tfrac{1}{12}(\tau c)^2$ for a solid). Extension 1 quantifies it.
- **Fuel.** Fuel is a large fraction of a wing's mass and moves; the tensor changes during the flight and Phase 09 needs it at several fuel states.
- **The rest of the aircraft.** One half-wing's tensor is a sum term; the aircraft's tensor is the sum over all components *after* each is shifted to the common CG by (10). The products of inertia of the two half-wings cancel by symmetry; the $I_{xx}$ values add.

## 3. Build the computational model

Open `starter/main.py`. The geometry is a frozen `WingBox` dataclass; `_check` validates it.

### Step 1: the map and its Jacobian

```python
def map_to_box(xi, eta, zeta, box=DEFAULT_BOX):
    y = eta * box.L
    c = box.c_root + (box.c_tip - box.c_root) * eta
    x = y * np.tan(box.sweep) + (box.f0 + (box.f1 - box.f0) * xi) * c                      # (6)
    z = (zeta - 0.5) * box.tau * c
    jac = (box.f1 - box.f0) * box.tau * box.L * c**2                                        # (7)
    return x, y, z, jac
```

### Step 2: the midpoint rule on the cube

```python
def midpoint_3d(f, n, box=DEFAULT_BOX):
    m = (np.arange(n) + 0.5) / n
    xi, eta, zeta = np.meshgrid(m, m, m, indexing="ij")
    x, y, z, jac = map_to_box(xi, eta, zeta, box)
    return float(np.sum(f(x, y, z) * jac) / n**3)                                           # (8)
```

`mass_properties_midpoint` calls it ten times (mass, three first moments, six tensor entries) with lambdas such as `lambda x, y, z: y**2 + z**2`.

### Step 3: the exact route

```python
def mass_properties_analytic(box=DEFAULT_BOX):
    xi, eta, zeta = sp.symbols("xi eta zeta", real=True)
    ...  # build x, y, z, jac as SymPy expressions of (xi, eta, zeta) exactly as in map_to_box
    def integ(expr):
        poly = sp.Poly(sp.expand(expr * jac), xi, eta, zeta)
        return float(sum(coef / ((a + 1) * (b + 1) * (c + 1)) for (a, b, c), coef in poly.terms()))
    ...  # m = rho * integ(1); cg = rho * [integ(x), integ(y), integ(z)] / m; I from (2)
```

Each monomial $\xi^a\eta^b\zeta^c$ integrates over the cube to $1/[(a + 1)(b + 1)(c + 1)]$, so expanding the polynomial and summing its terms *is* the exact integral, with no call to `sp.integrate` and in a few milliseconds.

### Step 4: closed forms, shifts, the ring, the oracle

`mass_closed_form` and `spanwise_cg_closed_form` are (3) and (4); `parallel_axis` and `inertia_about_cg` are (10) and its inverse; `ring_inertia_polar` is the midpoint rule on the $(r, \theta)$ rectangle with `dA = dr * dtheta` and the integrands `R * dA` and `R**3 * dA`; `tplquad_mass_and_ixx` hands `tplquad` the physical limits with $z$ innermost (`lambda z, x, y`), then $x$ with limits depending on $y$, then $y$ from 0 to $L$.

### Verify

Run `python -m pytest curriculum/phase-01/lesson-10-multiple-integrals` from the repo root. The tests, in order of how much they would hurt if wrong:

1. **The map.** `chord` is linear and raises outside $[0, L]$; the eight cube corners land on the box corners; the Jacobian equals (7) and equals the cross-section area times $L$.
2. **The rule.** `midpoint_3d` of 1 gives the volume for every $n$ to $10^{-3}$ and of a constant times the volume; $n = 0$ raises.
3. **Mass three ways.** Closed form, SymPy and `tplquad` agree to $10^{-10}$; midpoint at $n = 32$ to $10^{-4}$; the error divides by $3.9$–$4.1$ per doubling from $n = 4$ to 64.
4. **CG.** Inside the box; $z_{cg} = 0$; $y_{cg}$ equals (4); for an untapered box $y_{cg} = L/2$; for zero sweep $x_{cg} = \tfrac12(f_0 + f_1)\bar c$ with $\bar c$ the $c^2$-weighted mean chord.
5. **Tensor.** Symmetric, positive definite, triangle inequalities, $I_{xz} = I_{yz} = 0$ by the $z$ symmetry, $I_{xy} < 0$ because $x, y > 0$ throughout the box, $I_{zz} > I_{xx} > I_{yy}$ for a slender span; midpoint at $n = 32$ and `tplquad` match SymPy to $2 \times 10^{-4}$ and $10^{-8}$ relative.
6. **Parallel axis.** Shifting to the CG and back reproduces the tensor to $10^{-12}$; every diagonal entry about the CG is below that about any of twenty random points.
7. **Scaling.** Doubling $\rho$ doubles $I$; scaling all lengths by $k$ multiplies $m$ by $k^3$ and $I$ by $k^5$.
8. **Ring.** Polar midpoint equals (9) to $10^{-6}$ at $n = 64$; a thin ring tends to $mr^2$.
9. **Interface.** `main()` writes the figure.

### Compare

`main()` prints:

```
  mass: closed form 67.2000 kg, SymPy 67.2000 kg
  centre of mass (m): [0.86064 1.96429 0.     ]; closed-form y_cg = 1.96429 m (0.3929 L)
  inertia tensor about the root leading edge (kg m^2):
  [[ 384.0952 -127.6136   -0.    ]
   [-127.6136   53.825    -0.    ]
   [  -0.       -0.      437.7297]]
  inertia tensor about the centre of mass (kg m^2):
  [[124.8095 -14.0088  -0.    ]
   [-14.0088   4.0496  -0.    ]
   [ -0.      -0.     128.6686]]
  principal moments about the CG: [  2.4458 126.4133 128.6686]
  midpoint n =  2: mass 66.600000 kg (err -6.00e-01), y_cg 2.094595, I_xx 385.38058 (rel err 3.35e-03)
  ...
  midpoint n = 32: mass 67.197656 kg (err -2.34e-03), y_cg 1.964790, I_xx 384.09514 (rel err 2.36e-07)
  observed orders for the mass: [2. 2. 2. 2.]
  scipy tplquad: mass 67.200000 kg, I_xx 384.095232 kg m^2 (SymPy 384.095232)
  fuselage frame ...: polar midpoint m = 41.3512 kg, I = 39.33534 kg m^2; closed form 41.3512 kg, 39.33534 kg m^2
```

SymPy and `tplquad` agree on $I_{xx}$ to all nine printed figures, and both agree with the hand value 384.10 of (5). `tplquad` is the established tool: it is what you reach for on a region with curved limits and no polynomial structure, and it took a third of a second here against a few milliseconds for the term-by-term sum. The midpoint rule at $n = 32$ is within $2 \times 10^{-7}$ on $I_{xx}$ and $3 \times 10^{-5}$ on the mass; the inertia converges faster than the mass for a reason §4 asks you to find.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.**

**Experiment A: two convergence curves.** Tabulate the signed errors of the midpoint mass and of $I_{xx}$ for $n = 2, 4, \ldots, 128$ against the SymPy values. Prediction, written before running: the mass error is negative and divides by 4 per doubling (order 2); the $I_{xx}$ error is positive at small $n$, crosses zero near $n = 16$, and is negative and second order at large $n$. Then split $I_{xx}$ into its $\rho\int y^2\,\mathrm{d}V$ and $\rho\int z^2\,\mathrm{d}V$ parts with two extra `midpoint_3d` calls and tabulate each part's error. Compute the midpoint error coefficient $g'(1) - g'(0)$ for the mapped $y^2$ integrand $g(\eta) = \eta^2c(\eta)^2$ by hand and show it is zero *for this taper ratio* $c_t = c_r/2$; then repeat the experiment with $c_t = 0.6$ m and show the sign change disappears. Finally apply Richardson extrapolation $(4I_{2n} - I_n)/3$ to the $n = 32, 64$ pair and report the gain.

**Experiment B: what moves the tensor.** Sweep the taper ratio $c_t/c_r$ from 0.3 to 1.0 at fixed mass (adjust $\rho$ so that $m = 67.2$ kg) and the sweep angle from 0° to 30°, and tabulate $y_{cg}$, $I_{xx}$, $I_{zz}$ and $I_{xy}$ about the root. Prediction: taper moves $I_{xx}$ by tens of per cent at fixed mass (it changes where the mass is along the span); sweep leaves $I_{xx}$ and $y_{cg}$ exactly unchanged, raises $I_{zz}$ about the root by tens of per cent at 30°, and multiplies the product of inertia $|I_{xy}|$ several times over (it is not zero even unswept, because the box lies at $x > 0$). Decide, from the table, which two of the three parameters a flight-dynamics model must get right and which it can approximate.

Raw data format: the two error tables of Experiment A with the observed orders per row, and the parameter table of Experiment B, in your lesson notes.

## 5. Validate and improve

| quantity | hand (§2) | SymPy | midpoint $n = 32$ | `tplquad` | difference |
|---|---|---|---|---|---|
| $m$ (kg) | 67.20 | 67.2000 | 67.1977 | 67.2000 | $3 \times 10^{-5}$ (midpoint) |
| $y_{cg}$ (m) | 1.9643 | 1.96429 | 1.96479 | | $2.5 \times 10^{-4}$ |
| $x_{cg}$ (m) | 0.8606 | 0.86064 | | | $< 10^{-4}$ |
| $I_{xx}$ (kg m²) | 384.10 | 384.0952 | 384.0951 | 384.0952 | $2 \times 10^{-7}$ |
| $I_{xx,cg}$ (kg m²) | 124.8 | 124.8095 | | rod $mL^2/12 = 140.0$ | 11 % (taper, expected) |
| ring $m$, $I$ | 41.35, 39.34 | | polar $n = 64$: 41.3512, 39.33534 | closed 41.3512, 39.33534 | $< 10^{-6}$ |

Every residual has a cause. The midpoint mass is low by $3 \times 10^{-5}$ because $c^2$ is convex and the rule's $h^2$ error is $-2.34 \times 10^{-3}$ kg at $n = 32$; $y_{cg}$ is high by $2.5 \times 10^{-4}$ because the first moment and the mass have errors of different size and the ratio inherits the difference. $I_{xx}$ at $n = 32$ is closer than the mass because its dominant $y^2$ part has no $h^2$ error term at this taper ratio (Experiment A). The 11 % between $I_{xx,cg}$ and the uniform rod is not an error but the taper: a uniform rod is the wrong comparison and the right one, a rod with linear mass density matching $c^2$, gives 124.8 exactly.

What is still unverified when everything agrees: the smeared-density assumption and the geometry itself. Three integrators agreeing to nine figures on a solid box says nothing about whether a real wing of this planform has $I_{xx} = 384$ kg m² about its root; the extension with skins instead of a solid, and a mass breakdown of a real wing, would each move the number by more than every numerical error in the table combined. The honest statement is "$I_{xx} = 384$ kg m² for the solid-box model, exact to $10^{-7}$; model error unknown but bounded by the skin-versus-solid comparison for the small entries and by the mass calibration for the large ones".

## 6. Defend your work

Submit: the completed `starter/main.py` with all tests passing; `wing_box_mass.png`; the hand calculation of $m$, $y_{cg}$, $x_{cg}$ and $I_{xx}$ with the comparison table; the two experiment tables with one paragraph each; and a one-page note (sections 1, 2, 5, 6 and 7 of `projects/_template/report-template.md`) stating what the smeared-density model assumes and which tensor entries you would trust in a flight-dynamics model.

### Adversarial self-checks

- **What if the taper ratio goes to zero?** A triangular planform: (3) gives $m = \rho(f_1 - f_0)\tau Lc_r^2/3$ and (4) gives $y_{cg} = L/4$. Check with `WingBox(c_tip=1e-9)` (the code requires $c_t > 0$); the midpoint rule still converges at second order because $c^2$ is still a polynomial.
- **What if the box is untapered and unswept?** Then $I_{xy} = -m\,x_{cg}y_{cg}$ exactly (the CG shift of a body symmetric about its own centre), $y_{cg} = L/2$ and $I_{xx,cg} = m(L^2 + h^2)/12$ with $h = \tau c$. Problem 1 uses this as a hand check.
- **Which assumption fails first in the real world?** Smeared density, always: a real box is skins and spars, hollow or fuel-filled, and $\rho_{\text{eff}}$ differs between the structure and the fuel. The *mass* survives calibration; $I_{yy}$ and $I_{zz}$ about the CG do not.
- **What measurement could falsify the model?** A bifilar pendulum swing test of the wing gives $I_{xx}$ about the CG directly to about 2 %; a CG measurement on two scales gives $y_{cg}$ to a few millimetres. If the measured $y_{cg}$ is outboard of (4), the mass is not distributed as $c^2$ and a spanwise density $\rho(y)$ (Problem 6) is needed.
- **Does the result survive a change in integrator?** Midpoint, SymPy and `tplquad` must agree to $10^{-6}$ on every tensor entry; if they do not, the usual cause is a Jacobian error (a missing $c^2$ or $r$), which shows up as a mass error first.
- **Which decision changes if the uncertainty doubled?** If the effective density were known only to ±20 % instead of ±10 %, every entry scales with $\rho$ and the *ratios* $I_{xx} : I_{yy} : I_{zz}$ and the CG are unchanged. The decision that changes is whether to weigh the wing before the swing test: at ±10 % the model's $I_{xx}$ is as good as a rough swing test, at ±20 % the swing test wins and the model is only a check on its plausibility.

## Check yourself

1. State Fubini's theorem for a rectangle and say what changes for a general region.
2. Write $I_{xx}$, $I_{xy}$ and the full tensor as integrals; which entries are zero for a body symmetric about the $xy$-plane, and why?
3. For a linearly tapered box, why does the mass depend on $c_r^2 + c_rc_t + c_t^2$ and not on $(c_r + c_t)^2$?
4. The map (6) has off-diagonal entries from sweep and taper. Why do they not appear in the Jacobian?
5. In polar coordinates, where does the factor $r$ in $\mathrm{d}A = r\,\mathrm{d}r\,\mathrm{d}\theta$ come from, and what mass would you get for a ring if you forgot it?
6. State the parallel-axis theorem for the full tensor and show that the diagonal entries are smallest about the CG.
7. The midpoint mass error divides by 4 per doubling; the $I_{xx}$ error changes sign. Give the one-line reason.
8. The depth term in (5) was 0.025 % of $I_{xx}$. Is the same true for $I_{yy}$ about the CG? Compute the ratio from the printed tensor and the geometry.
