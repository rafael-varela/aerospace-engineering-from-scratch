# P01.07  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for $s = 15$; numerical results to four significant figures. Axes are right-handed; $\mu = 3.986\,004\,418 \times 10^{14}$ m³/s².

## 1. One state vector

**Assumptions.** Inertial Earth-centred axes, $z$ along the spin axis; two-body gravity.

**Derivation** ($s = 15$: $\mathbf r = (6000, 2500, 2000)$ km, $\mathbf v = (-3.0, 6.75, 2.5)$ km/s).
(a) $|\mathbf r| = \sqrt{36 + 6.25 + 4} \times 10^3 = 6800.7$ km; $|\mathbf v| = \sqrt{9 + 45.5625 + 6.25} = 7.798$ km/s.
(b) $h_1 = r_2v_3 - r_3v_2 = 2500 \times 2.5 - 2000 \times 6.75 = 6250 - 13\,500 = -7250$; $h_2 = r_3v_1 - r_1v_3 = 2000 \times (-3.0) - 6000 \times 2.5 = -6000 - 15\,000 = -21\,000$; $h_3 = r_1v_2 - r_2v_1 = 6000 \times 6.75 - 2500 \times (-3.0) = 40\,500 + 7500 = 48\,000$ km²/s. $|\mathbf h| = \sqrt{7250^2 + 21\,000^2 + 48\,000^2} = 52\,892$ km²/s $= 5.2892 \times 10^{10}$ m²/s.
(c) $\hat{\mathbf n} = (-0.1371, -0.3970, 0.9075)$; $i = \arccos 0.9075 = 24.84°$ (prograde).
(d) $\mathbf r \cdot \mathbf v = -18\,000 + 16\,875 + 5000 = 3875$ km²/s; $v_r = 3875/6800.7 = 0.5698$ km/s (projection of $\mathbf v$ on $\hat{\mathbf r}$); $v_t = \sqrt{7.798^2 - 0.5698^2} = 7.777$ km/s (rejection; check $|\mathbf r|v_t = 6800.7 \times 7.777 = 52\,890$ km²/s $= |\mathbf h|$ ✓); $\gamma = \operatorname{atan2}(0.5698, 7.777) = 4.190°$, climbing.
(e) $p = (5.2892 \times 10^{10})^2/3.986 \times 10^{14} = 7018.5$ km.
(f) $|\mathbf a| = \mu/r^2 = 3.986 \times 10^{14}/(6.8007 \times 10^6)^2 = 8.618$ m/s²; $a_t = \mathbf a \cdot \hat{\mathbf v} = -|\mathbf a|\sin\gamma = -0.6297$ m/s² (slowing while climbing); $a_n = |\mathbf a|\cos\gamma = 8.595$ m/s²; $R = v^2/a_n = 7.798^2 \times 10^6/8.595 = 7075$ km.

**Result.** $|\mathbf r| = 6801$ km, $|\mathbf v| = 7.798$ km/s, $\mathbf h = (-7250, -21\,000, 48\,000)$ km²/s, $|\mathbf h| = 5.289 \times 10^{10}$ m²/s, $i = 24.84°$, $v_r = 0.5698$ km/s, $v_t = 7.777$ km/s, $\gamma = 4.190°$, $p = 7018$ km, $a_t = -0.6297$ m/s², $a_n = 8.595$ m/s², $R = 7075$ km. The code agrees to every figure.

**Check.** Units: km × km/s = km²/s for $\mathbf h$ ✓; $p = h^2/\mu$ in km with $\mu$ in km³/s² ✓. $R \ne p$ here because the satellite is not at an apsis ($\gamma \ne 0$). From vis-viva $a = 1/(2/r - v^2/\mu) = 7066$ km and $e = \sqrt{1 - p/a} = 0.0821$: perigee altitude 108 km, apogee 1268 km. Problem 5(b) uses this.

**Common mistakes.** Mixing km in $\mathbf r$ with m/s in $\mathbf v$; computing $i$ from $\arctan$ of a ratio (loses the quadrant); taking $a_t$ as $|\mathbf a|\cos\gamma$ (that is $a_n$).

## 2. Three identities

**Derivation.** (a) $|\mathbf a \times \mathbf b|^2 = (a_2b_3 - a_3b_2)^2 + (a_3b_1 - a_1b_3)^2 + (a_1b_2 - a_2b_1)^2 = a_2^2b_3^2 + a_3^2b_2^2 + a_3^2b_1^2 + a_1^2b_3^2 + a_1^2b_2^2 + a_2^2b_1^2 - 2(a_2b_2a_3b_3 + a_1b_1a_3b_3 + a_1b_1a_2b_2)$. And $(\mathbf a \cdot \mathbf b)^2 = a_1^2b_1^2 + a_2^2b_2^2 + a_3^2b_3^2 + 2(a_1b_1a_2b_2 + a_1b_1a_3b_3 + a_2b_2a_3b_3)$. Adding, the mixed terms cancel and the nine squares $a_i^2b_j^2$ ($i, j = 1, 2, 3$) remain, which is $(a_1^2 + a_2^2 + a_3^2)(b_1^2 + b_2^2 + b_3^2)$. Then $|\mathbf a \times \mathbf b|^2 = |\mathbf a|^2|\mathbf b|^2 - |\mathbf a|^2|\mathbf b|^2\cos^2\theta = |\mathbf a|^2|\mathbf b|^2\sin^2\theta$, and since $0 \le \theta \le \pi$, $\sin\theta \ge 0$ and $|\mathbf a \times \mathbf b| = |\mathbf a||\mathbf b|\sin\theta$.

(b) $\mathbf b \times \mathbf c = (b_2c_3 - b_3c_2, b_3c_1 - b_1c_3, b_1c_2 - b_2c_1)$, so $\mathbf a \cdot (\mathbf b \times \mathbf c) = a_1(b_2c_3 - b_3c_2) - a_2(b_1c_3 - b_3c_1) + a_3(b_1c_2 - b_2c_1)$, which is the cofactor expansion along the first row of $\det[\mathbf a; \mathbf b; \mathbf c]$. Swapping two rows of a determinant changes its sign (each term of the expansion picks up one transposition); a cyclic permutation $(\mathbf a, \mathbf b, \mathbf c) \to (\mathbf b, \mathbf c, \mathbf a)$ is two swaps, so the sign is unchanged.

(c) First component: $\frac{\mathrm{d}}{\mathrm{d}t}(r_2v_3 - r_3v_2) = \dot r_2v_3 + r_2\dot v_3 - \dot r_3v_2 - r_3\dot v_2 = (v_2v_3 - v_3v_2) + (r_2a_3 - r_3a_2)$: the first bracket is the first component of $\mathbf v \times \mathbf v$ (zero), the second that of $\mathbf r \times \mathbf a$; the other components follow cyclically. For $\mathbf a = f(r)\,\mathbf r$ with $f$ a scalar, $\mathbf r \times \mathbf a = f(r)\,(\mathbf r \times \mathbf r) = \mathbf 0$; pulling $f$ out uses that each component of the cross product is bilinear in the components, which holds only for a scalar factor (a matrix factor would not pull out). Hence $\dot{\mathbf h} = \mathbf 0$.

**Check.** (a) for $\mathbf a = \hat{\mathbf i}$, $\mathbf b = \hat{\mathbf j}$: $1 + 0 = 1 \cdot 1$ ✓. (b) $\hat{\mathbf i} \cdot (\hat{\mathbf j} \times \hat{\mathbf k}) = 1 = \det I$ ✓.

**Common mistakes.** In (a), forgetting the factor 2 on the mixed terms; in (c), writing $\mathbf a \times \mathbf r$ instead of $\mathbf r \times \mathbf a$ (same conclusion, wrong sign in the middle).

## 3. The wind triangle, in 3D

**Assumptions.** Flat Earth, steady wind, axes $x$ east, $y$ north, $z$ up. TAS $V = 120 \times 0.514444 = 61.733$ m/s; wind speed $15.433$ m/s.

**Derivation.** (a) Blowing toward 210°: $\mathbf w = 15.433(\sin 210°, \cos 210°, 0) = (-7.717, -13.366, 0)$ m/s: a headwind component of 7.72 m/s (against east) and a crosswind of 13.37 m/s from the north.
(b) Air velocity $V(\sin\psi, \cos\psi, 0)$. North component of ground velocity $= V\cos\psi - 13.366 = 0 \Rightarrow \cos\psi = 0.21652$, $\psi = 77.50°$. Ground speed $= V\sin\psi - 7.717 = 60.269 - 7.717 = 52.55$ m/s $= 102.2$ kt. Crab angle: `angle_between(air, ground)` $= 12.50°$ ($= 90° - 77.50°$). The rejection of $\mathbf w$ from $\hat{\mathbf x}$ is $(0, -13.366, 0)$, the crosswind the heading cancels ($\sin 12.50° \times 61.733 = 13.37$ ✓); the projection is $(-7.717, 0, 0)$, the headwind that reduces ground speed ($60.27 - 7.72$) — the along-track part of the air velocity is $V\cos 12.50°$, not $V$.
(c) Horizontal air speed $V_h = \sqrt{61.733^2 - 3^2} = 61.660$ m/s; $\cos\psi = 13.366/61.660$, $\psi = 77.48°$; ground velocity $(61.660\sin\psi - 7.717, 0, 3) = (52.48, 0, 3.00)$ m/s, $|\mathbf v_g| = 52.56$ m/s; climb angle $\arcsin(3/52.56) = 3.272°$ over the ground (against $\arcsin(3/61.73) = 2.786°$ through the air); angle between air and ground velocity 12.51°.
(d) Time to a waypoint: project the ground velocity onto the track direction (dot); off-track distance: the rejection, $|\mathbf d \times \hat{\mathbf t}|$ (cross).

**Result.** $\mathbf w = (-7.717, -13.366, 0)$ m/s; heading 077.5°, ground speed 102.2 kt, crab 12.50°; with climb: heading 077.5°, $\mathbf v_g = (52.48, 0, 3.00)$ m/s, climb angle 3.272°, air–ground angle 12.51°.

**Check.** No wind: $\psi = 90°$, ground speed 120 kt ✓. Pure headwind of 30 kt: ground speed 90 kt ✓ (set the crosswind to zero in the formula). Units: knots converted once, at the start.

**Common mistakes.** Using "from 030°" as the direction the wind blows toward (heading comes out 102.5°, wrong side); subtracting the full wind speed from the TAS for the ground speed; forgetting that the climb takes part of the TAS out of the horizontal plane (small here: 0.07 m/s).

## 4. Lift distributions

**Derivation.** (a) $M_x = \int_0^s y\,l(y)\,\mathrm{d}y$. Triangular: $(L/s)\int_0^s (y - y^2/s)\,\mathrm{d}y = (L/s)(s^2/2 - s^2/3) = Ls/6 = (L/2)(s/3)$. Uniform: $(L/2s)(s^2/2) = Ls/4 = (L/2)(s/2)$. Elliptic: $(L/2)(4s/3\pi)$ (lesson). With $L/2 = 7750$ N, $s = 8$ m: triangular 20 667 N m, uniform 31 000 N m, elliptic 26 314 N m. `resultant` with 400 strips: 20 666.7, 31 000.0 (both exact to rounding: the midpoint rule is exact for linear integrands), 26 314.4. Ratios: elliptic/triangular $= 8/(\pi \cdot 2) \cdot 3/... = 4s/(3\pi) \div s/3 = 4/\pi = 1.273$; elliptic/uniform $= 8/(3\pi) = 0.8488$.
(b) Uniform strips, elliptic case, error in $M_x$: 698, 241, 82.9, 20.4 (at $n = 5, 10, 20, 50$) ... 0.305 N m at $n = 800$; observed orders 1.54, 1.54, 1.53, 1.52, 1.52, 1.51, 1.51 for $n = 5 \to 640$: order 1.5, as predicted for a $\sqrt{\cdot}$ end-point singularity. Cosine spacing: relative errors $1.80 \times 10^{-2}$, $4.82 \times 10^{-3}$, $1.24 \times 10^{-3}$, $3.16 \times 10^{-4}$, $7.97 \times 10^{-5}$, $2.00 \times 10^{-5}$ for $n = 5, 10, 20, 40, 80, 160$: ratios 3.9, 3.9, 3.9, 4.0, 4.0, order 2. The substitution $y = s\sin\phi$ turns $\sqrt{1 - y^2/s^2}\,\mathrm{d}y$ into $\cos^2\phi\,\mathrm{d}\phi$, smooth.
(c) 0.1 % on $M_x$: uniform strips need $n = 50$ (0.078 %); cosine-spaced need $n = 80$ (0.0080 %) or, interpolating, about 25 for 0.1 % — $n = 40$ gives 0.032 %. Note the cosine grid is *worse* at small $n$ (1.8 % at $n = 5$ against 2.7 % uniform is comparable) and wins only through its order; with $n = 400$ it reaches $3 \times 10^{-6}$ against $3 \times 10^{-5}$.

**Result.** Closed forms $Ls/6$, $Ls/4$, $2Ls/3\pi$; ratios 1.273 and 0.8488; orders 1.5 (uniform) and 2.0 (cosine); for 0.1 %: 50 uniform strips or 40 cosine strips.

**Check.** The three lever arms $s/3 < 4s/3\pi < s/2$ order as the distributions concentrate lift outboard ✓. The lift itself is exact for triangular and uniform (midpoint exact for polynomials of degree $\le 1$) ✓.

**Common mistakes.** Normalising the triangular distribution to $L$ instead of $L/2$ over the half-span; reporting the order from two points only; comparing cosine and uniform grids at the same $n$ without stating the error level required.

## 5. Where the models stop

(a) **Rigid wing.** A 1.5 m tip deflection on an 8 m semi-span gives a bending slope at the tip of about $2 \times 1.5/8 = 0.375$ rad; on a 25° swept wing, bending alone rotates the streamwise sections by roughly slope $\times\sin\Lambda = 0.16$ rad $\approx 9°$ of washout at the tip, far more than the few degrees of twist a wing is built with. The tip unloads, the lift moves inboard, and $M_x$ falls by tens of percent relative to the rigid elliptic estimate. Measure: tip twist with a laser pointer and a target, or root strain against lift.
(b) **No drag.** Problem 1's orbit has perigee at 108 km altitude, where the density is of order $10^{-7}$ kg/m³ — ten thousand times the 300 km value of `P01.06`. The satellite would lose tens of metres per second per pass and re-enter within a few revolutions; $\mathbf h$ is not even approximately constant. Measure: the apogee radius from successive fixes drops orbit to orbit.
(c) **Measurement precision.** A 1 m/s error along $z$ shifts $i$ by 0.0048° (0.0007° along $x$, 0.0021° along $y$); the sixth decimal place ($10^{-6}$ degrees) is $10^4$ times below the uncertainty. Quote $24.84 \pm 0.005°$. Measure: repeat the fix; the scatter tells you the honest precision.
(d) **Axis convention.** With $x$ aft, $y$ right, $z$ up and lift at $x = y\tan\Lambda > 0$, $\mathbf r \times \mathbf F = (y\tan\Lambda, y, 0) \times (0, 0, l) = (yl, -yl\tan\Lambda, 0)$: $M_y$ is negative, nose-down torsion, and the *sign* is forced by the right-hand rule once the axes are declared. A "convention" could relabel $+M_y$ as nose-down, but then $\tan\Lambda$ must appear with the opposite sign in that convention's formula; the physical statement "the root is twisted leading-edge down by $12.3$ kN m" has no sign freedom, and the spar cap in tension is the lower rear one either way. Measure: strain gauges on the front and rear spar caps.

**Common mistakes.** Calling the 108 km perigee "low but fine"; treating six decimals as free; arguing about signs without first writing down the axes.

## 6. The orbital plane from three fixes

**Assumptions.** Unperturbed two-body orbit; fixes exact to the 0.1 km given.

**Derivation.** (a) `plane_from_points` gives $\hat{\mathbf n} = (-0.33927, 0.29787, 0.89228)$ and $d = 0.034$ km. $d$ is the distance of the fitted plane from the Earth's centre; a two-body orbit's plane contains the attracting centre, so $d$ should vanish. The 34 m comes from rounding the fixes to 0.1 km, as the sensitivity in (e) confirms.
(b) $i = \arccos 0.89228 = 26.839°$; the lesson's state vector gives 26.839°: same orbit.
(c) `angle_between` gives 40.00° between fixes 1 and 2 and 40.00° between 2 and 3; $|\mathbf p| = 6800.7$ km for all three, consistent with a circular orbit. At 10 min $= 600$ s per 40°, $T = 600 \times 9 = 5400$ s $= 90$ min; the circular period at this radius is $2\pi\sqrt{r^3/\mu} = 93.0$ min. The 3 % mismatch means the fixes are not exactly 10 min apart (they are 10.33 min apart if the orbit is circular); the problem's "about ten minutes" is the loose input.
(d) Station $= 6378.137(\cos 40°\cos 20°, \cos 40°\sin 20°, \sin 40°) = (4591.3, 1671.1, 4099.8)$ km; signed distance $\hat{\mathbf n} \cdot \mathbf s - d = 2598$ km (on the side the normal points to). Overhead passes require the station to lie in the orbital plane at some instant; as the Earth rotates, the station's latitude is fixed at 40° but the plane reaches at most latitude $i = 26.8°$. The satellite never passes overhead; its maximum elevation from that station is well below the zenith.
(e) Perturbing $\mathbf p_2$ by 1 km along $x$, $y$, $z$ turns the normal by 0.0122°, 0.0107°, 0.0321°. Perturbing $\mathbf v$ by 1 m/s along each axis in `orbital_plane` turns it by 0.0025°, 0.0022°, 0.0067°. The three-point method is about five times more sensitive to a 1 km position error than the state-vector method is to a 1 m/s velocity error; which is "more trustworthy" depends on which the tracking system delivers better. Radar gives range to tens of metres but range-rate to centimetres per second, so for radar the state vector wins by a wide margin; for angles-only optical fixes with kilometre-level errors, spreading the fixes further apart (longer lever arms) is the remedy, since the sensitivity scales as (position error)/(fix spacing).

**Result.** $\hat{\mathbf n} = (-0.3393, 0.2979, 0.8923)$, $d = 0.03$ km; $i = 26.84°$; 40.00° per fix; $T \approx 93$ min (circular, from $|\mathbf p|$); station 2598 km from the plane, no overhead passes; sensitivities 0.01–0.03° per km against 0.002–0.007° per m/s.

**Check.** $\hat{\mathbf n} \cdot \mathbf p_k = d$ for all three fixes to the metre ✓. $|\hat{\mathbf n}| = 1$ ✓. The normal from fixes equals $\hat{\mathbf h}$ from the lesson to five figures ✓.

**Common mistakes.** Ordering the points so the normal points to $-z$ and reporting $i = 153°$ (retrograde) for a prograde orbit — check the sign of $\hat n_3$ against the sense of motion; using the geodetic latitude on a spherical Earth (fine here, stated); concluding from $d \ne 0$ that the orbit is perturbed instead of first checking the input precision.

## Check-yourself answers

1. $\mathbf a \cdot \mathbf b = 2$; $\mathbf a \times \mathbf b = (-3, 0, 1)$; $\cos\theta = 2/\sqrt{14}$, $\theta = 58°$.
2. $\mathbf a \times \mathbf b$ is by construction perpendicular to $\mathbf a$, and perpendicular vectors have zero dot product.
3. $\mathbf F \times \mathbf F = \mathbf 0$: adding $\lambda\mathbf F$ to $\mathbf r$ adds $\lambda\mathbf F \times \mathbf F$ to the moment.
4. Retrograde; $i > 90°$ (because $\cos i = h_z/|\mathbf h| < 0$).
5. $|\mathbf h|^2 = \mu p$, i.e. $p = |\mathbf h|^2/\mu$.
6. $a_t = 0$ (speed constant); $a_n = |\mathbf a \times \hat{\mathbf v}|$, the cross product.
7. At $4/(3\pi) = 0.4244$ of the semi-span; smaller than uniform (which acts at $0.5s$) by the factor $8/(3\pi) = 0.849$.
8. Near $\theta = 0$ the cosine is $1 - \theta^2/2$, which rounds to 1 once $\theta^2/2$ is below machine epsilon, so `arccos` returns 0; the sine (from the cross product) is $\theta$ itself. The difference shows below about $10^{-7}$ rad and `arccos` returns exactly 0 at $10^{-8}$.
