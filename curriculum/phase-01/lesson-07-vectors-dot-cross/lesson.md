# P01.07  Vectors, Dot & Cross Products

> The dot product asks "how much along?"; the cross product asks "how much around?". Every force, moment and orbit is one of those two questions.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the component formulas for the dot and cross products from their geometric definitions, and prove the identities (Lagrange, triple products, transmissibility of a moment) that the tests check. 2. Implement and test a vector library from components, vectorized over stacks of vectors, and compare it with NumPy's `dot`, `cross` and `linalg.norm`. 3. Compute, from one position–velocity pair, a satellite's angular momentum, orbital plane, inclination, semi-latus rectum and flight-path angle, and the root bending and torsion moments of a swept wing from its lift distribution. 4. Decide when the arccos formula for an angle is untrustworthy and when a strip model of a distributed load is converged. |
| **Prerequisites** | `P01.06` (polar conic, $h = r^2\dot\theta$), `P01.01` (product rule), `P00.03` (NumPy broadcasting); Bridge B1 (trigonometry, $\sin$, $\cos$, $\operatorname{atan2}$). Skill check: you can expand a $3 \times 3$ determinant and you know that $\sin^2 + \cos^2 = 1$. |
| **Engineering question** | From a single tracked position and velocity of a satellite, which plane is it orbiting in and at what inclination; and from the lift distribution on a swept wing, what bending and torsion moments must the wing root carry? |
| **Assumptions and limits** | Right-handed Cartesian axes fixed in an inertial frame (satellite) or in the aircraft (wing); forces act at points (strip loads concentrated at strip centres); the wing is rigid, so loads do not move when it bends; gravity is the only force on the satellite, so $\mathbf r \times \mathbf v$ is constant. Fails when the wing twists under load (the lift distribution changes), when the strips are too coarse near the tip where the elliptic distribution has a vertical tangent, and for the satellite whenever a non-central force (drag, $J_2$, thrust) acts. |
| **Analytical work** | Dot and cross products from geometry to components; the Lagrange identity; the scalar triple product as a determinant; moment of a force and the principle of transmissibility; $\mathrm{d}(\mathbf r \times \mathbf v)/\mathrm{d}t = 0$ for a central force; inclination from $\mathbf h$; elliptic lift root bending moment $2Ls/(3\pi)$ by hand. |
| **Computational deliverable** | `starter/main.py`: `dot`, `cross`, `norm`, `unit`, `angle_between`, `projection`, `rejection`, `scalar_triple`, `moment`, `resultant`, `specific_angular_momentum`, `orbital_plane`, `semi_latus_rectum`, `flight_path_angle`, `tangential_normal_components`, `point_line_distance`, `plane_from_points`, `point_plane_distance`, `elliptic_lift_distribution`, `wing_root_loads`. Tests prove: agreement with NumPy on single and stacked vectors; symmetry, bilinearity, antisymmetry, perpendicularity, right-handedness; the Lagrange and Jacobi identities and the BAC–CAB rule; the triple product equals a determinant; projection + rejection = original; transmissibility and Varignon's theorem; the angular momentum of the `P01.06` conic is the same vector at nine true anomalies; the root loads match the closed form; a circular orbit has purely normal gravity with $v^2/a_n = r$. |
| **Practical activity** | **computational experiment**: measure the convergence of the strip model's root moment with strip count for the elliptic distribution and explain the observed order; measure where `arccos` loses the angle between nearly parallel vectors. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `vector_mechanics.png`, the hand calculation of $\mathbf h$ and the inclination with its comparison to the code, the strip-convergence table and the arccos table. |
| **Assessment** | `assessments/phase-01/grading-rubric.md`; competencies 1–4. Pass: the three proofs in Problem 2 complete; all public tests pass; the hand calculation of $\mathbf h$, $i$ and the root moment agrees with the code to four figures; the convergence experiment reports an observed order. |
| **References** | `@openstax-calculus-volume-3` §2.1 "Vectors in the Plane" and §2.2 "Vectors in Three Dimensions" (required-free); §2.3 "The Dot Product" (required-free); §2.4 "The Cross Product" (required-free); §2.5 "Equations of Lines and Planes in Space" (required-free); §3.1 "Vector-Valued Functions and Space Curves", §3.2 "Calculus of Vector-Valued Functions" and §3.4 "Motion in Space" (required-free); §2.6 "Quadric Surfaces" and §3.3 "Arc Length and Curvature" (supplementary-free). |
| **Extensions** | ★ From $\mathbf r$, $\mathbf v$ and $\mathbf h$ compute the eccentricity vector $\mathbf e = \mathbf v \times \mathbf h/\mu - \hat{\mathbf r}$ and show its magnitude matches $e$ from vis-viva and the direction points to periapsis. ★ Replace the elliptic distribution by a lifting-line result (Phase 10) and see how the root torsion changes with twist. |
| **Time** | unmeasured |

## 1. The engineering problem

A ground station tracks a satellite for a few minutes and hands the flight dynamics team one position vector $\mathbf r$ and one velocity vector $\mathbf v$ in Earth-centred axes with $z$ along the spin axis. From those six numbers the team must say which plane the orbit lies in, how steeply it is tilted to the equator (the **inclination**), whether the satellite is climbing or descending at that instant, and the size of the orbit. All of it follows from one vector, $\mathbf h = \mathbf r \times \mathbf v$, and the two products this lesson is about.

The second half of the problem comes from the structures office. A wing of 16 m span carries 15.5 kN of lift in level flight, distributed along the span, with the quarter-chord line swept back 25°. The root fitting must carry the total lift, the **bending moment** that tries to fold the wing upward, and the **torsion** that the sweep introduces because the outboard lift acts behind the root. Each strip of wing contributes a force at a point; the moment of each about the root is a cross product; the sum is the design load. Get the sign of the torsion wrong and the spar is designed to twist the wrong way.

Success criteria, written before we start: the library agrees with NumPy's `dot`, `cross` and `norm` to $10^{-14}$ relative on twenty random vectors; every identity in §2 holds to the same precision; the angular momentum computed at nine points of the `P01.06` conic is one vector to $10^{-12}$; the strip model reproduces the closed-form root bending moment of an elliptic distribution to $5 \times 10^{-4}$ with 400 strips; and the hand calculation of $\mathbf h$ and the inclination agrees with the code to four figures.

## 2. Understand and derive

### 2.1 Vectors, components and length

A **vector** has a magnitude and a direction; a position, a velocity, a force, a moment. In right-handed axes with unit vectors $\hat{\mathbf i}, \hat{\mathbf j}, \hat{\mathbf k}$, a vector is three numbers, $\mathbf a = a_1\hat{\mathbf i} + a_2\hat{\mathbf j} + a_3\hat{\mathbf k}$, written $(a_1, a_2, a_3)$ (`@openstax-calculus-volume-3` §2.2). Its length is $|\mathbf a| = \sqrt{a_1^2 + a_2^2 + a_3^2}$ by Pythagoras applied twice, and $\hat{\mathbf a} = \mathbf a/|\mathbf a|$ is the **unit vector** in its direction (undefined for the zero vector; our `unit` raises). Addition and scaling are component-wise. Units: all three components carry the vector's unit (m, m/s, N), and a unit vector is dimensionless.

```
        z (k)                      right-handed: curl the fingers of the right hand
        │                          from x toward y; the thumb points along z.
        │     a = (a1, a2, a3)     i x j = k,  j x k = i,  k x i = j
        │    /
        │   /  |a| = sqrt(a1² + a2² + a3²)
        │  /
        │ /
        └──────────── y (j)
       /
      /
     x (i)
```

### 2.2 The dot product: "how much along?"

**Definition (geometric):** $\mathbf a \cdot \mathbf b = |\mathbf a||\mathbf b|\cos\theta$, with $\theta \in [0, \pi]$ the angle between the vectors (`@openstax-calculus-volume-3` §2.3). It is a scalar. In components, expand $(a_1\hat{\mathbf i} + a_2\hat{\mathbf j} + a_3\hat{\mathbf k}) \cdot (b_1\hat{\mathbf i} + \dots)$ using $\hat{\mathbf i} \cdot \hat{\mathbf i} = 1$, $\hat{\mathbf i} \cdot \hat{\mathbf j} = 0$ (unit length, mutually perpendicular):

$$\mathbf a \cdot \mathbf b = a_1 b_1 + a_2 b_2 + a_3 b_3.$$   (1)

Consequences, each a test: $\mathbf a \cdot \mathbf b = \mathbf b \cdot \mathbf a$; $(2\mathbf a + \mathbf c) \cdot \mathbf b = 2\,\mathbf a \cdot \mathbf b + \mathbf c \cdot \mathbf b$ (bilinear); $\mathbf a \cdot \mathbf a = |\mathbf a|^2$; $\mathbf a \cdot \mathbf b = 0$ exactly when the vectors are perpendicular. The **projection** of $\mathbf a$ onto $\mathbf b$ is the part of $\mathbf a$ along $\hat{\mathbf b}$,

$$\operatorname{proj}_{\mathbf b}\mathbf a = (\mathbf a \cdot \hat{\mathbf b})\,\hat{\mathbf b},$$   (2)

and the **rejection** $\mathbf a - \operatorname{proj}_{\mathbf b}\mathbf a$ is the part perpendicular to $\mathbf b$. Work is a dot product: a force $\mathbf F$ moving its point of application by $\mathbf d$ does $W = \mathbf F \cdot \mathbf d$ joules; only the component along the motion counts.

### 2.3 The cross product: "how much around?"

**Definition (geometric):** $\mathbf a \times \mathbf b$ is the vector perpendicular to both $\mathbf a$ and $\mathbf b$, of length $|\mathbf a||\mathbf b|\sin\theta$ (the area of the parallelogram they span), pointing the way a right-handed screw advances when turned from $\mathbf a$ toward $\mathbf b$ (`@openstax-calculus-volume-3` §2.4). From $\hat{\mathbf i} \times \hat{\mathbf j} = \hat{\mathbf k}$, $\hat{\mathbf j} \times \hat{\mathbf k} = \hat{\mathbf i}$, $\hat{\mathbf k} \times \hat{\mathbf i} = \hat{\mathbf j}$, $\hat{\mathbf i} \times \hat{\mathbf i} = \mathbf 0$ and antisymmetry, expanding the product of two component sums gives

$$\mathbf a \times \mathbf b = (a_2 b_3 - a_3 b_2,\; a_3 b_1 - a_1 b_3,\; a_1 b_2 - a_2 b_1),$$   (3)

which is the formal determinant with rows $(\hat{\mathbf i}, \hat{\mathbf j}, \hat{\mathbf k})$, $(a_1, a_2, a_3)$, $(b_1, b_2, b_3)$. Each component is a cyclic shift of the previous one (1→2→3→1), which is how to remember it and how the code writes it. Properties, each a test: $\mathbf a \times \mathbf b = -\mathbf b \times \mathbf a$; $\mathbf a \times \mathbf a = \mathbf 0$; $(\mathbf a \times \mathbf b) \cdot \mathbf a = 0$.

**Lagrange's identity** connects the two products:

$$|\mathbf a \times \mathbf b|^2 + (\mathbf a \cdot \mathbf b)^2 = |\mathbf a|^2|\mathbf b|^2,$$   (4)

which is just $\sin^2\theta + \cos^2\theta = 1$ multiplied by $|\mathbf a|^2|\mathbf b|^2$; Problem 2 proves it directly from the components, which is the proof that (3) really has length $|\mathbf a||\mathbf b|\sin\theta$. It gives a numerically robust way to find the angle:

$$\theta = \operatorname{atan2}\big(|\mathbf a \times \mathbf b|,\ \mathbf a \cdot \mathbf b\big).$$   (5)

Why not $\arccos(\mathbf a \cdot \mathbf b/|\mathbf a||\mathbf b|)$? Near $\theta = 0$ the cosine is $1 - \theta^2/2$; for $\theta = 10^{-8}$ that is $1 - 5 \times 10^{-17}$, which rounds to exactly 1 in double precision, and `arccos` returns 0. The sine, $\theta$ itself, is represented perfectly. Part 4 measures this.

### 2.4 Triple products

The **scalar triple product** $\mathbf a \cdot (\mathbf b \times \mathbf c)$ is the signed volume of the parallelepiped on $\mathbf a, \mathbf b, \mathbf c$ (base area $|\mathbf b \times \mathbf c|$ times the height $\mathbf a \cdot \hat{\mathbf n}$), and in components it is the $3 \times 3$ determinant with rows $\mathbf a, \mathbf b, \mathbf c$. Hence it is unchanged by a cyclic permutation, flips sign under a swap, and vanishes when the three vectors are coplanar. The **vector triple product** obeys the BAC–CAB rule $\mathbf a \times (\mathbf b \times \mathbf c) = \mathbf b(\mathbf a \cdot \mathbf c) - \mathbf c(\mathbf a \cdot \mathbf b)$, and summing its cyclic versions gives the Jacobi identity $\mathbf a \times (\mathbf b \times \mathbf c) + \mathbf b \times (\mathbf c \times \mathbf a) + \mathbf c \times (\mathbf a \times \mathbf b) = \mathbf 0$. The tests check both on random vectors; you will meet BAC–CAB again when the inertia tensor appears in Phase 04.

### 2.5 Moments of forces

The **moment** (torque) of a force $\mathbf F$ applied at position $\mathbf r$ about the origin is

$$\mathbf M = \mathbf r \times \mathbf F \quad [\mathrm{N\,m}].$$   (6)

Its magnitude is $|\mathbf F|$ times the perpendicular distance from the origin to the line of action (the lever arm), because $|\mathbf r||\mathbf F|\sin\theta = |\mathbf F| \cdot (|\mathbf r|\sin\theta)$. Sliding the force along its own line, $\mathbf r \to \mathbf r + \lambda\mathbf F$, changes nothing since $\mathbf F \times \mathbf F = \mathbf 0$: the **principle of transmissibility**. Moments add (Varignon's theorem: the moment of a resultant equals the sum of the moments), and moving the reference point from $A$ to $B$ changes the total moment by $\mathbf M_B = \mathbf M_A - \mathbf r_{BA} \times \mathbf F_{\text{total}}$. Both are tests of `resultant`.

**The wing.** Axes: $x$ aft, $y$ along the right wing, $z$ up. A strip at station $y$ of width $\mathrm{d}y$ carries lift $l(y)\,\mathrm{d}y$ upward at the quarter-chord point $(y\tan\Lambda, y, 0)$ for sweep $\Lambda$. Its moment about the root is $(y\tan\Lambda, y, 0) \times (0, 0, l\,\mathrm{d}y) = (y\,l\,\mathrm{d}y,\ -y\tan\Lambda\,l\,\mathrm{d}y,\ 0)$. So

$$M_x = \int_0^s y\,l(y)\,\mathrm{d}y \ \text{(root bending)}, \qquad M_y = -\tan\Lambda\,M_x \ \text{(root torsion)},$$   (7)

with $s$ the semi-span. The torsion is negative: positive $M_y$ would pitch the nose up about the $y$ axis, and the aft-acting lift of a swept wing instead twists the root leading-edge *down*. For an **elliptic distribution** $l(y) = (4L/\pi b)\sqrt{1 - (2y/b)^2}$ with $b = 2s$ (chosen so $\int_{-s}^{s} l\,\mathrm{d}y = L$; this is the minimum-induced-drag distribution you will derive in Phase 10), $\int_0^s y\sqrt{1 - y^2/s^2}\,\mathrm{d}y = s^2/3$ (substitute $u = 1 - y^2/s^2$), so

$$M_x = \frac{4L}{2\pi s}\cdot\frac{s^2}{3} = \frac{2Ls}{3\pi} = \frac{L}{2}\cdot\frac{4s}{3\pi}.$$   (8)

The half-wing lift $L/2$ acts at $4s/(3\pi) = 0.4244\,s$ from the root, the centroid of a quarter-ellipse. For $L = 15.5$ kN and $s = 8$ m: $M_x = 7750 \times 3.395 = 26\,314$ N m and $M_y = -\tan 25° \times 26\,314 = -12\,270$ N m.

### 2.6 Angular momentum and the orbital plane

A satellite at $\mathbf r$ with velocity $\mathbf v$ has **specific angular momentum** (per unit mass)

$$\mathbf h = \mathbf r \times \mathbf v \quad [\mathrm{m^2/s}].$$   (9)

Differentiate with the product rule for cross products (`@openstax-calculus-volume-3` §3.2): $\dot{\mathbf h} = \dot{\mathbf r} \times \mathbf v + \mathbf r \times \dot{\mathbf v} = \mathbf v \times \mathbf v + \mathbf r \times \mathbf a = \mathbf r \times \mathbf a$. **Assumption (central force):** gravity points along $-\mathbf r$, so $\mathbf r \times \mathbf a = \mathbf 0$ and $\mathbf h$ is constant. This is the fact `P01.06` imported: $|\mathbf h| = r^2\dot\theta$ (the transverse velocity is $r\dot\theta$ and $|\mathbf r \times \mathbf v| = r \cdot r\dot\theta$), and $h^2 = \mu p$ gives the semi-latus rectum $p = |\mathbf h|^2/\mu$ from one state vector. Because $\mathbf h$ is perpendicular to both $\mathbf r$ and $\mathbf v$ at every instant and never changes, the motion stays in the plane through the origin with normal $\hat{\mathbf h}$: the **orbital plane**. Its tilt from the equator is the **inclination**,

$$\cos i = \hat{\mathbf h} \cdot \hat{\mathbf k} = h_z/|\mathbf h|, \qquad 0 \le i \le \pi,$$   (10)

with $i < 90°$ for prograde orbits (eastward) and $i > 90°$ for retrograde. The **flight-path angle** of `P01.06` is the angle between $\mathbf v$ and the local horizontal: $\tan\gamma = (\mathbf v \cdot \hat{\mathbf r})/|\mathbf v - (\mathbf v \cdot \hat{\mathbf r})\hat{\mathbf r}|$, projection over rejection.

**Worked example.** $\mathbf r = (3500, -5000, 3000)$ km, $\mathbf v = (6.02, 4.67, 0.73)$ km/s. By (3):
$h_1 = (-5000)(0.73) - (3000)(4.67) = -3650 - 14\,010 = -17\,660$;
$h_2 = (3000)(6.02) - (3500)(0.73) = 18\,060 - 2555 = 15\,505$;
$h_3 = (3500)(4.67) - (-5000)(6.02) = 16\,345 + 30\,100 = 46\,445$ km²/s.
$|\mathbf h| = \sqrt{17\,660^2 + 15\,505^2 + 46\,445^2} = \sqrt{2.7094 \times 10^9} = 52\,052$ km²/s $= 5.2052 \times 10^{10}$ m²/s. Inclination: $\cos i = 46\,445/52\,052 = 0.89228$, $i = 26.84°$. Semi-latus rectum: $p = (5.2052 \times 10^{10})^2/3.986 \times 10^{14} = 6797.3$ km, against $|\mathbf r| = 6800.7$ km: nearly circular. $\mathbf r \cdot \mathbf v = 21\,070 - 23\,350 + 2190 = -90$ km²/s, so the radial speed is $-90/6800.7 = -0.0132$ km/s and $\gamma = -0.10°$: descending, barely. Sanity check: $|\mathbf v| = 7654$ m/s against circular speed $\sqrt{\mu/r} = 7656$ m/s at this radius. This is a 420 km orbit at 26.8°, nothing exotic.

### 2.7 Tangential and normal acceleration

Along any path the acceleration splits into a part along $\hat{\mathbf v}$ that changes the speed and a part perpendicular to it that turns the path (`@openstax-calculus-volume-3` §3.4):

$$a_t = \mathbf a \cdot \hat{\mathbf v}, \qquad a_n = |\mathbf a \times \hat{\mathbf v}| = \frac{v^2}{R},$$   (11)

with $R$ the radius of curvature. Again projection (dot) and rejection (cross). For a circular orbit gravity is purely normal and $v^2/a_n$ returns the orbit radius; for the example state, $a_t = +0.015$ m/s² (it is descending, so gravity has a small component along the velocity and the satellite is speeding up, consistent with $\gamma < 0$) and $v^2/a_n = 6797$ km, equal to $p$ to four figures because $R = p$ at the apsides of any conic and this orbit is nearly circular. For an aircraft, $a_n/g$ is the load factor's horizontal part.

### 2.8 Lines and planes

A line through $\mathbf p_0$ with direction $\mathbf d$ is $\mathbf p_0 + t\mathbf d$; the distance from a point $\mathbf q$ to it is $|(\mathbf q - \mathbf p_0) \times \hat{\mathbf d}|$, the rejection of $\mathbf q - \mathbf p_0$ from $\mathbf d$ (`@openstax-calculus-volume-3` §2.5). A plane through three points has normal $\mathbf n = (\mathbf p_2 - \mathbf p_1) \times (\mathbf p_3 - \mathbf p_1)$ and equation $\hat{\mathbf n} \cdot \mathbf x = \hat{\mathbf n} \cdot \mathbf p_1 = d$; the signed distance of a point is $\hat{\mathbf n} \cdot \mathbf q - d$. Three tracked positions therefore give the orbital plane with no velocity at all (Problem 6).

### 2.9 Where the model stops being valid

- **Non-central forces.** Drag, $J_2$ and thrust all have $\mathbf r \times \mathbf a \ne \mathbf 0$; then $\mathbf h$ drifts and the "orbital plane" is an instantaneous one. Two state vectors a day apart give two slightly different planes.
- **Rigid wing.** A real swept wing bends and twists under load, moving the quarter-chord line and changing $l(y)$. The root loads here are first-iteration numbers for the structures loop.
- **Point loads.** Concentrating a strip's lift at its centre is exact for a uniform strip and second-order accurate otherwise; near the tip the elliptic $l(y)$ has infinite slope and the error decays more slowly (part 4).
- **Axes.** Everything assumes a right-handed, orthonormal basis. A left-handed convention flips every cross product; mixed units (km and m/s) multiply into km·m/s nonsense. Check both before trusting a sign.

## 3. Build the computational model

Open `starter/main.py`. Every function takes a single vector (shape `(3,)`) or a stack (shape `(N, 3)`) and returns the matching shape; the private helpers `_as_vectors`, `_pair` and `_out` are provided and handle that bookkeeping so the mathematics stays visible.

### Step 1: the two products from components

```python
def dot(a, b):
    A, B, single = _pair(a, b)                      # both (N, 3), broadcast
    return _out(A[:, 0] * B[:, 0] + A[:, 1] * B[:, 1] + A[:, 2] * B[:, 2], single)

def cross(a, b):
    A, B, single = _pair(a, b)
    out = np.column_stack([A[:, 1] * B[:, 2] - A[:, 2] * B[:, 1],      # cyclic: 2,3
                           A[:, 2] * B[:, 0] - A[:, 0] * B[:, 2],      #         3,1
                           A[:, 0] * B[:, 1] - A[:, 1] * B[:, 0]])     #         1,2
    return _out(out, single)
```

Write the components out, as in (1) and (3); the whole point is that you can compare them with `np.dot` and `np.cross`, which you will use from Phase 02 on.

### Step 2: length, direction, angle

```python
def norm(a):
    A = _as_vectors(a)
    return _out(np.sqrt(np.sum(A * A, axis=1)), np.asarray(a).ndim == 1)

def unit(a):
    ...  # raise ValueError on a zero vector; divide by n[:, None] to broadcast over rows

def angle_between(a, b):
    A, B, single = _pair(a, b)
    ...  # raise if either is zero
    c = _as_vectors(cross(A, B))
    return _out(np.arctan2(np.sqrt(np.sum(c * c, axis=1)), np.sum(A * B, axis=1)), single)   # (5)
```

`projection` is (2) with `unit(onto)`; `rejection` subtracts it; `scalar_triple(a, b, c)` is `dot(a, cross(b, c))`.

### Step 3: mechanics

```python
def moment(r, force):
    return cross(r, force)                                             # (6)

def resultant(points, forces, about=(0.0, 0.0, 0.0)):
    P, F = _as_vectors(points), _as_vectors(forces)
    arm = P - np.asarray(about, dtype=float)
    return F.sum(axis=0), _as_vectors(cross(arm, F)).sum(axis=0)

def orbital_plane(r, v):
    h = specific_angular_momentum(r, v)                                # cross(r, v)
    n = h / norm(h)                                                    # raise if |h| == 0
    return n, float(np.arctan2(np.hypot(n[0], n[1]), n[2]))            # (10), via atan2

def tangential_normal_components(v, a):
    vhat = unit(v)
    return float(dot(a, vhat)), float(norm(cross(a, vhat)))            # (11)
```

`flight_path_angle` is `atan2(dot(v, rhat), norm(rejection(v, r)))`; `point_line_distance`, `plane_from_points` and `point_plane_distance` are §2.8 line by line.

### Step 4: the wing

```python
def elliptic_lift_distribution(y, lift_total, semi_span):
    b = 2.0 * semi_span
    return 4.0 * lift_total / (np.pi * b) * np.sqrt(1.0 - (2.0 * y / b) ** 2)

def wing_root_loads(lift_total, semi_span, sweep_rad, n_strips=50):
    edges = np.linspace(0.0, semi_span, n_strips + 1)
    y = 0.5 * (edges[:-1] + edges[1:])                                 # strip centres
    lift = elliptic_lift_distribution(y, lift_total, semi_span) * np.diff(edges)
    points = np.column_stack([y * np.tan(sweep_rad), y, np.zeros_like(y)])
    forces = np.column_stack([np.zeros_like(y), np.zeros_like(y), lift])
    return resultant(points, forces)
```

No special wing code: strips are just points and forces, and `resultant` does the rest. That is the sign of a vector library doing its job.

### Verify

Run `python -m pytest curriculum/phase-01/lesson-07-vectors-dot-cross` from the repo root. The tests, in order of how much they would hurt if wrong:

1. **Against NumPy.** `dot`, `cross`, `norm` equal `np.dot`, `np.cross`, `np.linalg.norm` on twenty random vectors, singly and stacked, with the right output shapes.
2. **Identities.** Symmetry, bilinearity, antisymmetry, $\mathbf a \times \mathbf a = \mathbf 0$, perpendicularity, right-handedness of $\hat{\mathbf i}, \hat{\mathbf j}, \hat{\mathbf k}$, Lagrange (4), the triple product as a determinant with its cyclic and swap behaviour, BAC–CAB and Jacobi, projection + rejection = original with Pythagoras.
3. **Robustness.** The angle between $(1, 0, 0)$ and $(1, 10^{-9}, 0)$ is $10^{-9}$ to $10^{-6}$ relative; zero vectors and wrong shapes raise `ValueError`.
4. **Mechanics.** Transmissibility; $|\mathbf M| = |\mathbf F| \times$ lever arm; Varignon; the reference-point shift rule; the root loads match (8) and the torsion vanishes without sweep; $\mathbf h$ is the same vector at nine true anomalies of the `P01.06` conic and gives back $p$; equatorial and polar orbits give $i = 0$ and $90°$; $\gamma$ is invariant under a rotation of the axes; a circular orbit's gravity is purely normal with $v^2/a_n = r$.
5. **Interface.** `main()` writes the figure.

### Compare

NumPy is the established tool here. `main()` prints `max |cross - np.cross| = 0.0e+00` and `|dot - np.dot| = 0.0e+00`: identical, because both evaluate the same six multiplications in the same order. `np.cross` is what you will call from now on; the reason to have written your own once is that you now know exactly what it computes, in which order, and with which sign convention (`np.cross` also takes an `axis` argument and handles 2-vectors by returning a scalar, two conveniences that have produced real bugs). `np.linalg.norm` differs in one way worth knowing: for very large or very small components it scales before squaring to avoid overflow, which `sqrt(a·a)` does not; at $10^{200}$ our `norm` returns `inf`. Not a concern in SI mechanics, but the kind of difference a comparison exists to find.

The wing output of `main()`:

```
Right half-wing root loads (elliptic lift, 50 strips)
  F = [   0.     0.  7752.4] N  (should be 0, 0, 7750)
  M = [ 26334.  -12279.8      0. ] N m: bending M_x = 26334.0, torsion M_y = -12279.8
  analytic bending (L/2)(4 b/2 / (3 pi)) = 26313.6 N m; torsion = -tan(sweep) * bending = -12270.2 N m
```

Fifty strips overshoot the lift by 0.03 % and the bending by 0.08 %. Both errors have the same cause, examined in part 4.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.**

**Experiment A: strip convergence.** Prediction, written before running: the midpoint rule is second-order for a smooth integrand, so doubling the strip count should divide the root-moment error by 4. But $l(y) \propto \sqrt{1 - y^2/s^2}$ has an infinite slope at the tip, which usually degrades a midpoint rule to order 1.5. Run `wing_root_loads(15.5e3, 8.0, np.radians(25), n)` for $n = 5, 10, 20, \dots, 800$, tabulate $M_x - 2Ls/(3\pi)$, and compute the observed order $\log_2(\text{err}_n/\text{err}_{2n})$. Then decide: how many strips for 0.1 % on the bending moment? For 0.01 %?

**Experiment B: where arccos fails.** For $\theta = 10^{-2}, 10^{-4}, 10^{-6}, 10^{-8}, 10^{-9}$ build $\mathbf b = (\cos\theta, \sin\theta, 0)$ and compute the angle to $\hat{\mathbf i}$ both by `np.arccos(dot/(|a||b|))` and by `angle_between`. Prediction: arccos is fine to about $10^{-6}$, loses its fourth significant figure at $10^{-7}$ (where $\theta^2/2 \approx 5 \times 10^{-15}$, only twenty machine epsilons below 1) and returns exactly 0 at $10^{-8}$ and below; `angle_between` is exact at every size. Why would a star tracker, which measures the angle between two star directions to 1 arcsecond ($5 \times 10^{-6}$ rad), care?

Raw data format: two plain-text tables, $(n, F_z, M_x, \text{error}, \text{order})$ and $(\theta, \theta_{\arccos}, \theta_{\operatorname{atan2}})$, kept in your lesson notes.

## 5. Validate and improve

Hand calculation (§2.6), code and an independent route side by side:

| quantity | hand | code | independent check | difference |
|---|---|---|---|---|
| $\mathbf h$ (km²/s) | $(-17\,660, 15\,505, 46\,445)$ | same | `np.cross` identical | 0 |
| $|\mathbf h|$ (m²/s) | $5.2052 \times 10^{10}$ | $5.2052 \times 10^{10}$ | | 0 |
| $i$ | 26.84° | 26.839° | $\arccos(h_z/|\mathbf h|)$ | 0 |
| $p$ (km) | 6797.3 | 6797.3 | $v^2/a_n = 6797.4$ km | 0.001 % |
| $\gamma$ | $-0.10°$ | $-0.099°$ | sign of $\mathbf r \cdot \mathbf v$ | 0 |
| root $M_x$ (N m) | 26 314 (closed form) | 26 334 (50 strips) | 26 314.4 (400 strips) | 0.08 % → 0.003 % |
| root $M_y$ (N m) | $-12\,270$ | $-12\,280$ | $-\tan\Lambda\,M_x$ exactly | same ratio |

The strip error is the only non-trivial residual, and it has a cause you can quantify: in Experiment A the observed order is close to 1.5, not 2, because of the tip singularity, and the lift error ($+0.03$ % at 50 strips) is always positive because the midpoint of the last strip over-estimates a concave $\sqrt{\cdot}$ profile. Two fixes exist. Cluster the strips toward the tip with a cosine spacing $y = s\sin(\pi k/2n)$, which is what Phase 10's lifting-line code does; or integrate $l(y)$ exactly over each strip, which is possible here ($\int\sqrt{1 - u^2}\,\mathrm{d}u$ is elementary) but not for a measured distribution. Try the cosine spacing and report the new order.

What is still unverified when everything agrees: the physics inputs. The elliptic distribution is a model, not a measurement, and the quarter-chord line is where thin-airfoil theory puts the lift, not where a pressure survey would. The satellite's inclination is exact *given* $\mathbf r$ and $\mathbf v$, but a tracking solution carries uncertainty; propagate $\delta\mathbf v = (0, 0, 10)$ m/s through `orbital_plane` and the inclination moves by 0.014°, so a radar fix good to 1 m/s pins $i$ to about 0.0014° and a precise orbit solution good to 0.01 m/s to $10^{-5}$ degrees. That, not the arithmetic, is the precision of the answer. Agreement between hand and code proves the code implements the vectors; the uncertainty analysis is what makes the number usable.

## 6. Defend your work

Submit: the completed `starter/main.py` with all tests passing; `vector_mechanics.png`; the hand calculation of $\mathbf h$, $i$, $p$, $\gamma$ and the root moments with the comparison table; the two experiment tables with one paragraph each; and a one-page note (sections 1, 2, 5, 6 and 7 of `projects/_template/report-template.md`) stating what each model assumes and the uncertainty of the inclination.

### Adversarial self-checks

- **What if the lift doubles?** Every moment doubles and the lever arm $4s/(3\pi)$ does not move: moments are linear in force. If instead the span doubles at fixed lift, $M_x$ doubles too (lever arm $\propto s$) while the torsion/bending ratio $\tan\Lambda$ is unchanged.
- **What if $\mathbf v \to 0$? What if $\mathbf v \parallel \mathbf r$?** $\mathbf h \to \mathbf 0$ and the orbital plane is undefined; the code must raise, not return a normal made of `nan`. Physically that is a body falling straight down, which has no orbit.
- **Which assumption fails first in the real world?** For the wing, rigidity: a swept wing's bending produces washout that unloads the tip and moves the lift inboard, reducing $M_x$ by a few percent at limit load. For the satellite, the central force: $J_2$ turns the plane by a fraction of a degree per day.
- **What measurement could falsify the model?** Strain gauges at the root give $M_x$ and $M_y$ directly; a measured $M_y/M_x$ ratio different from $-\tan\Lambda$ means the lift does not act at the quarter chord. For the satellite, a second state vector a day later gives a second $\mathbf h$; the angle between them is the perturbation, zero under the model.
- **Does the result survive a change in axes?** Rotate $\mathbf r$ and $\mathbf v$ by the same orthogonal matrix: $|\mathbf h|$, $p$, $\gamma$ and $a_t, a_n$ must not change, and $i$ changes only because the $z$ axis moved. The flight-path test does exactly this with a random rotation.
- **Which decision changes if the uncertainty doubled?** A velocity error of 2 m/s instead of 1 m/s gives 0.0028° instead of 0.0014° on the inclination, which at 6800 km is a cross-track position uncertainty of 330 m instead of 170 m. A conjunction screen with a 1 km miss-distance threshold is unaffected; one with a 200 m threshold now cannot decide whether to manoeuvre. The decision changes only when the uncertainty crosses the threshold it is compared with.

## Check yourself

1. State the dot and cross products of $(1, 2, 3)$ and $(0, 1, 0)$ and the angle between them to the nearest degree.
2. Why is $(\mathbf a \times \mathbf b) \cdot \mathbf a$ always zero? Give the one-line geometric reason.
3. The moment of a force about a point is unchanged when the force slides along its line of action. Which product property is that?
4. A satellite has $h_z < 0$. Is it prograde or retrograde, and what can you say about $i$?
5. Write the relation between $|\mathbf h|$, $\mu$ and the semi-latus rectum $p$ of `P01.06`.
6. In a coordinated level turn the acceleration is purely normal. What is $a_t$, and which product gives $a_n$ from $\mathbf a$ and $\mathbf v$?
7. For an elliptic lift distribution, at what fraction of the semi-span does the half-wing lift act? Is the root bending moment larger or smaller than for a uniform distribution of the same total lift, and by what factor?
8. Why does `angle_between` use `atan2(|a × b|, a · b)` rather than `arccos`? At roughly what angle does the difference appear in double precision?
