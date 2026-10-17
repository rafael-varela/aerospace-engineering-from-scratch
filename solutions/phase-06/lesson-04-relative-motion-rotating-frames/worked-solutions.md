# P06.04  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe script that calls it, given to the figures shown.

## 1. The fictitious terms on an airliner (s = 15: λ = 35° N)

**Assumptions.** Spherical Earth, Ω = 7.2921159 × 10⁻⁵ rad/s, the aircraft at the surface radius.

(a) The fictitious Coriolis acceleration is −2Ω × v with Ω = Ω(0, cos λ, sin λ) in ENU.
Flying north, v = (0, 250, 0): **(+0.02091, 0, 0) m/s²**, all of it east, i.e. sideways to the right; size 2ΩV sin λ.
Flying east, v = (250, 0, 0): **(0, −0.02091, +0.02987) m/s²**: 0.02091 m/s² to the south (to the right again) and 0.02987 m/s² upward (2ΩV cos λ, the Eötvös term: eastbound aircraft weigh slightly less).

(b) Centrifugal −Ω × (Ω × r): **(0, −0.01592, +0.02273) m/s²**, toward the equator and outward.

(c) atan(0.02091/9.80665) = **0.122°** of bank. No pilot would notice it, and an autopilot holding the track absorbs it. An inertial navigator that ignored it would integrate 0.0209 m/s² twice: ½at² = 135 km of error after an hour. That is why an inertial navigation system must model the Earth's rotation.

**Common mistakes.** Using +2Ω × v (the term in a_P) instead of −2Ω × v (what the ground observer sees); taking Ω along local up instead of along the Earth's axis.

## 2. Derivations

(a) A unit vector e fixed in B rotates rigidly with B: in time dt it turns through ω dt, so de = ω dt × e, i.e. (de/dt)_I = ω × e. For A = ΣAᵢeᵢ: (dA/dt)_I = ΣȦᵢeᵢ + ΣAᵢ ω × eᵢ = (dA/dt)_B + ω × A.

(b) R_P = R_O + r. By (1), v_P = v_O + v_rel + ω × r with v_rel = (dr/dt)_B. Differentiate again in I. (dv_rel/dt)_I = a_rel + **ω × v_rel** (first half). (d(ω × r)/dt)_I = α × r + ω × (dr/dt)_I = α × r + ω × (v_rel + ω × r), which contains the second **ω × v_rel**. Collecting: a_P = a_O + α × r + ω × (ω × r) + 2ω × v_rel + a_rel.

(c) Dot (3) with v: v · v̇ = v · f − 2v · (Ω × v) − v · (Ω × (Ω × r)). The middle term is zero (Ω × v ⟂ v). f = −∇(−μ/r) and −Ω × (Ω × r) = ∇(½|Ω × r|²), so v · f = −d(−μ/r)/dt and the last term is d(½|Ω × r|²)/dt. Hence d/dt(½v² − μ/r − ½|Ω × r|²) = 0.

(d) r̈₁ = −2Ω × (v₀ + gt) with r₁(0) = ṙ₁(0) = 0. Integrate: ṙ₁ = −2Ω × (v₀t + gt²/2), r₁ = −Ω × (v₀t² + gt³/3). At the equator Ω = (0, Ω, 0); for v₀ = (0, v_n, v_u) and g = (0, 0, −g), Ω × v₀ = (Ωv_u, 0, 0) and Ω × g = (−Ωg, 0, 0), so the east component is r₁,E = −Ωv_u t² + Ωg t³/3. At T = 2v_u/g: −ΩT²(v_u − gT/3) = −ΩT²(v_u − 2v_u/3) = **−ΩT²v_u/3**. Its vertical part is zero (Ω × anything in the y–z plane has only an x component), so there is no time correction.

(e) Write r = Q s with Q = R_z(−Ω_z t), so Q̇ = −Ω_z (e_z ×) Q, and e_z × commutes with Q. Then ṙ = Q(ṡ − Ω_z e_z × s) and r̈ = Q(s̈ − 2Ω_z e_z × ṡ − Ω_z² s) for horizontal s. The right side of (7), −ω₀²r − 2Ω_z e_z × ṙ, is Q(−ω₀² s − 2Ω_z e_z × ṡ − 2Ω_z² s). Equating, the e_z × ṡ terms cancel and **s̈ = −(ω₀² + Ω_z²) s**. Release from rest on the floor: r(0) = (a, 0) and ṙ(0) = 0 ⇒ s(0) = (a, 0), ṡ(0) = Ω_z e_z × s(0) = (0, Ω_z a). The solution s = (a cos Wt, (Ω_z a/W) sin Wt) is a fixed ellipse; r = R_z(−Ω_z t)s turns it at −Ω_z: clockwise seen from above in the north. The tests confirm (8) satisfies (7) to 10⁻⁶ by finite differences.

**Common mistakes.** Writing a single ω × v_rel; using the full Ω instead of Ω sin λ for the pendulum.

## 3. Dropping a spanner in a spinning station

**Assumptions.** Rigid station spinning at constant ω about its axis; no air; the spanner released at rest relative to the station at radius r = 98 m.

(a) ω²R = g ⇒ **ω = √(g/R) = 0.3132 rad/s = 2.99 rpm** (period 20.06 s). Walking spinward at 1.5 m/s: (ωR + v)²/R = **1.098 g**; anti-spinward: **0.906 g**. (The Coriolis term 2ωv is the main part of the change.)

(b) Inertial frame: the spanner starts at (r, 0) with velocity (0, ωr) = (0, 30.69 m/s) and flies straight. It reaches the floor where r² + (ωrt)² = R²: **t = √(R² − r²)/(ωr) = 0.6484 s**. It has then turned through atan(ωt) = 0.20033 rad about the axis, while the floor point below the release has turned through ωt = 0.20306 rad. It lands **0.272 m anti-spinward** (behind) of the point directly below the release.

(c) The rotating-frame integration (f = 0, only the fictitious terms) gives a fall time of 0.648427 s and an arc offset of −0.27238 m, agreeing with (b) to 10⁻¹² m.

(d) On Earth a 2 m drop takes √(2h/g) = 0.6387 s; in the station 0.6484 s, 1.5 % longer. The apparent gravity at 98 m is only 0.98 g, and it is not a uniform field: the spanner is really flying a straight line while the floor curves up to meet it, and part of its apparent fall is sideways (the Coriolis term).

**Common mistakes.** Giving the spanner the floor's speed ωR; expecting it to land directly below.

## 4. Your own long-range shot (s = 15: λ = 35° N, azimuth 150°)

(a) Flight **146.933 s**, range **103.329 km**; rotating and inertial impacts **2.6 × 10⁻⁵ m** apart.

(b) Coriolis shift: down-range **+287.8 m**, cross-range **+892.8 m** (to the right); first order: +292.8 m and +868.6 m. Agreement to 2–3 %, as in the lesson.

(c) Against a non-rotating Earth: +680.1 m down-range, +977.4 m cross-range. The non-Coriolis part, **392.4 m down-range and 84.6 m cross-range**, comes from the centrifugal field: at 35° N it pulls toward the equator and outward, and firing south-south-east, "toward the equator" is mostly along the line of fire.

(d) J changes by **−1.2 × 10⁻¹¹** relative between launch and impact.

## 5. When does Earth's rotation matter?

(a) Flat formula, 45°, north from 49° N: **20.7 m at 330 m/s (11.1 km range)**, **575 m at 1000 m/s (102 km)**, **2354 m at 1600 m/s (261 km)**. It passes 20 m at about **11 km** of vacuum range.

(b) The airliner measures its position continuously (satellite navigation, radio aids) and its autopilot holds the track by banking about 0.1°; the Coriolis term is absorbed in the feedback loop. An inertial navigator measures only acceleration and integrates it twice; any acceleration it fails to model becomes a position error growing as t².

(c) No air. Drag shortens the flight and the time of flight; since the deflection grows as T² (or faster), the real deflection is smaller than the vacuum value for the same muzzle speed, though large in proportion to the (shorter) range.

(d) Take-off distance: no (seconds of flight and nothing to integrate; the effect is 10⁻³ g). 10 km mortar: marginal (tens of metres at most, comparable with other errors; corrected in precise fire). Ballistic missile: yes, essential (kilometres). 10-hour inertial navigation: yes, essential (the error would be thousands of kilometres). Weather forecast: yes, the Coriolis term is what organises large-scale winds.

## Check-yourself answers

1. (dA/dt)_I = (dA/dt)_B + ω × A: the left side is the rate seen in the inertial frame, (dA/dt)_B the rate of A's components in the rotating frame.
2. Origin a_O, Euler α × r, centripetal ω × (ω × r), Coriolis 2ω × v_rel, relative a_rel.
3. It is perpendicular to the relative velocity, so its power v · a is zero.
4. To the left of its track (southern hemisphere: left, whatever the heading).
5. Because Ω has a horizontal component (north) at the equator: the rising shell is pushed west and the falling shell east, with a net westward shift of ΩT²v_up/3.
6. Because the gun moves with the ground at Ω × r₀; v₀ is measured relative to the ground.
7. 2π/(Ω sin 30°) = 2 × 23.93 h = 47.9 h.
8. Firing east, the Coriolis acceleration −2Ω × v has an upward component 2ΩV cos λ, so the shell stays up longer and lands further away (and a westward shot lands short).
