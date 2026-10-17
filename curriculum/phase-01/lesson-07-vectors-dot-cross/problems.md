# P01.07  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and the axis convention you used.

Several problems use your personal seed **s**: the day of the month you were born (1–31). Use $\mu = 3.986\,004\,418 \times 10^{14}$ m³/s², $R_E = 6378.137$ km, $g = 9.80665$ m/s².

## 1. (Warm-up: one state vector)

A tracking solution gives, in Earth-centred axes with $z$ along the spin axis, $\mathbf r = (6000,\ 1000 + 100s,\ 2000)$ km and $\mathbf v = (-3.0,\ 7.2 - 0.03s,\ 2.5)$ km/s.

Find by hand, then check with your code: (a) $|\mathbf r|$ and $|\mathbf v|$; (b) $\mathbf h = \mathbf r \times \mathbf v$ and $|\mathbf h|$; (c) the unit normal of the orbital plane and the inclination; (d) the radial and transverse speeds (projection and rejection of $\mathbf v$ on $\mathbf r$) and the flight-path angle; (e) the semi-latus rectum $p = |\mathbf h|^2/\mu$; (f) the tangential and normal components of the gravitational acceleration $-\mu\mathbf r/r^3$ and the radius of curvature $v^2/a_n$. Answer form: numbers with units, four significant figures; the component arithmetic of (b) written out.

## 2. (Derivation: three identities)

(a) Prove Lagrange's identity $|\mathbf a \times \mathbf b|^2 + (\mathbf a \cdot \mathbf b)^2 = |\mathbf a|^2|\mathbf b|^2$ directly from the component formulas (1) and (3), and conclude that $|\mathbf a \times \mathbf b| = |\mathbf a||\mathbf b|\sin\theta$.
(b) Prove that $\mathbf a \cdot (\mathbf b \times \mathbf c)$ equals the determinant with rows $\mathbf a, \mathbf b, \mathbf c$, and hence that it is unchanged by a cyclic permutation and changes sign under a swap of any two vectors.
(c) For a body with velocity $\mathbf v = \dot{\mathbf r}$ and acceleration $\mathbf a = \dot{\mathbf v}$, prove the product rule $\mathrm{d}(\mathbf r \times \mathbf v)/\mathrm{d}t = \mathbf v \times \mathbf v + \mathbf r \times \mathbf a$ from the components, and conclude that $\mathbf h$ is constant for any central force $\mathbf a = f(r)\,\mathbf r$. State where you use that $f$ is a scalar.

Answer form: three proofs, each step justified.

## 3. (Unfamiliar scenario: the wind triangle, in 3D)

An aircraft flies at a true airspeed of 120 kt (1 kt $= 0.514444$ m/s) and must track due east (090°) over the ground. The wind blows *from* 030° at 30 kt. Use axes $x$ east, $y$ north, $z$ up; a bearing $\beta$ is the direction $(\sin\beta, \cos\beta, 0)$.

(a) Write the wind velocity vector (remember: "from 030°" means it blows toward 210°).
(b) The ground velocity is the vector sum of the air velocity (airspeed along the heading) and the wind. Find the heading that makes the ground velocity point east, the resulting ground speed in knots, and the crab angle between heading and track. Use `angle_between` for the crab angle and `projection`/`rejection` of the wind on the track direction to explain the two numbers.
(c) Now the aircraft also climbs at 3 m/s with the same 120 kt true airspeed. Recompute the heading and the ground velocity vector, and find the climb angle relative to the ground and the angle between the air velocity and the ground velocity.
(d) In one sentence: which product (dot or cross) answers "how long to reach a waypoint 50 km east", and which answers "how far off track am I after 10 minutes if I fly the no-wind heading"?

Answer form: three vectors, four angles, one speed, one sentence.

## 4. (Numerical, uses your code: lift distributions)

For the lesson's wing ($L = 15.5$ kN, $s = 8$ m, $\Lambda = 25°$):

(a) Replace the elliptic distribution by a triangular one, $l(y) = (L/s)(1 - y/s)$ (which also integrates to $L/2$ over the half-span), and by a uniform one, $l(y) = L/(2s)$. Derive the closed-form root bending moment for each and compute all three with `resultant` using 400 strips. Report the ratios of elliptic to triangular and elliptic to uniform bending.
(b) Measure the observed convergence order of `wing_root_loads` for the elliptic case with $n = 5, 10, \dots, 640$ strips, then repeat with cosine-spaced strip edges $y_k = s\sin(\pi k/2n)$ and report the new order.
(c) A structures engineer asks for the root loads to 0.1 %. How many uniform strips do you need, and how many cosine-spaced ones?

Answer form: three closed forms, two ratios, two observed orders with the tables, two strip counts.

## 5. (Judgment: where the models stop)

For each situation name the assumption that fails, estimate the size of the effect with one number, and say what you would measure to detect it:

(a) the wing of Problem 4 is a sailplane wing with a tip deflection of 1.5 m at limit load;
(b) the satellite of Problem 1 (compute its perigee altitude from $a = 1/(2/r - v^2/\mu)$ and $e = \sqrt{1 - p/a}$ before you answer);
(c) the inclination of Problem 1 is quoted to six decimal places from a radar fix with a velocity uncertainty of 1 m/s;
(d) a colleague computes the root torsion as $M_y = +\tan\Lambda\,M_x$ and argues that the sign is "a convention".

Answer form: four short paragraphs, each with a number.

## 6. (Hard: the orbital plane from three position fixes)

Three position fixes of a satellite, about ten minutes apart, are $\mathbf p_1 = (3500.0, -5000.0, 3000.0)$ km, $\mathbf p_2 = (6123.3, -1168.6, 2718.4)$ km and $\mathbf p_3 = (5881.4, 3209.7, 1164.8)$ km.

(a) With `plane_from_points`, find the unit normal and the offset $d$ of the plane through the three fixes; explain what a non-zero $d$ would mean physically and why it should be close to zero.
(b) Find the inclination and compare with the lesson's state vector, which belongs to the same orbit.
(c) Find the angle swept between consecutive fixes with `angle_between`, and estimate the orbital period if the orbit is circular.
(d) A ground station sits at latitude 40° N, longitude 20° E on a spherical Earth of radius $R_E$ (at the instant of the fixes, Earth-centred axes coincide with the station's Earth-fixed frame). Find its signed distance from the orbital plane and say whether the satellite can ever pass directly overhead.
(e) Estimate the sensitivity of the normal to a 1 km error in $\mathbf p_2$ by perturbing it along each axis in turn, and compare with the sensitivity of `orbital_plane` to a 1 m/s error in $\mathbf v$. Which method of finding the plane would you trust more, and why?

Answer form: one unit vector and offset, two angles, one period, one distance with a yes/no, a short comparison.
