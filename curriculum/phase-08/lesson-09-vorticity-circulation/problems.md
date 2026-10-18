# P08.09  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact (for the Lamb–Oseen vortex, Kelvin's theorem, Biot–Savart) and which come from a model whose assumptions you should name (a horseshoe vortex, laminar diffusion).

Several problems use your personal seed **s**: the day of the month you were born (1–31).

Conventions as in the lesson: Γ counterclockwise positive in plane problems; in 3-D, a filament's Γ follows the right-hand rule about its direction. Air at sea level: ρ = 1.225 kg/m³, ν = 1.46 × 10⁻⁵ m²/s; at 3000 m: ρ = 0.90913 kg/m³, ν = 1.863 × 10⁻⁵ m²/s.

## 1. (Warm-up: a Lamb–Oseen vortex)

A Lamb–Oseen vortex has Γ = 10 + s m²/s and core radius r_c = 0.5 m (so r_c² = 0.25 m²).

(a) Find the radius of peak swirl, 1.1209 r_c, and the peak swirl speed.
(b) Find the circulation within r = 0.5, 1 and 2 m from the closed form, as a percentage of Γ, and check each with `circulation_polygon` on a 2000-point circle.
(c) Find the vorticity on the axis. Show that near the axis the fluid rotates as a solid body at half that rate.
(d) Check Stokes's theorem on the off-centre disk of radius 0.4 m centred at (0.5, 0.3) m: line integral against area integral.

Answer form: a radius and a speed; three circulations, three percentages and three checks; a vorticity and a rotation rate; two numbers.

## 2. (Derivation: Kelvin, Lamb–Oseen and Biot–Savart)

(a) Prove Kelvin's theorem: for a material loop C(t) in an inviscid fluid with p a function of ρ only (barotropic) and a conservative body force, dΓ/dt = 0. Differentiate Γ = ∮ u·dl following the particles: dΓ/dt = ∮ (Du/Dt)·dl + ∮ u·du, the second term because each element dl of a material loop is stretched at the rate du. Use the Euler equation for Du/Dt and show that each term is the integral of an exact differential round a closed loop. Say which term a viscous fluid adds, and why it is not zero in general.
(b) In plane flow the vorticity equation is Dω/Dt = ν∇²ω. Show that for an axisymmetric ω(r, t) with no radial velocity the convective term vanishes, and verify by substitution that ω = Γ/(π r_c²) exp(−r²/r_c²) with r_c² = r_c0² + 4νt is a solution. Derive u_θ from ω and the circulation within r.
(c) From the Biot–Savart law for a straight segment, derive the speed Γ/(2πd) at distance d from an infinite line vortex, and Γ/(4πd) at distance d abeam the end of a semi-infinite one. Then derive the horseshoe vortex's downwash at the centre of its bound vortex, −Γ/(πb′), and far downstream, −2Γ/(πb′).

Answer form: three derivations, each with the assumptions listed.

## 3. (Unfamiliar scenario: crossing an airliner's wake)

An airliner of 70 t on approach at 75 m/s at sea level, wing span 34 m, sheds two trailing vortices with an elliptic-loading spacing b′ = πb/4.

(a) Find Γ₀ = W/(ρVb′) and the descent rate Γ₀/(2πb′) of the vortex pair (each vortex carries the other down).
(b) The light aeroplane (span 10 m, 60 m/s) flies along the axis of one vortex. Treating the vortex as a line vortex, find the vertical velocity at its wing tips, 5 m either side of the axis, and the change in angle of attack it induces at each tip. Why is this dangerous?
(c) Treat each vortex as a laminar Lamb–Oseen vortex starting with zero core. Find the core radius of peak swirl and the peak swirl after 60 and 120 s, and the time for laminar diffusion alone to halve the circulation within 5 m of the axis. Do the numbers make physical sense?
(d) Air-traffic separation rules space aircraft on approach by minutes, not days. What decays the real vortices, and which assumption of the Lamb–Oseen model fails?

Answer form: a circulation and a speed; a velocity, an angle and a sentence; four numbers, a time and a sentence; two or three sentences.

## 4. (Numerical, uses your code: Kelvin's theorem for your own loop)

Take a material loop that is initially a circle of radius 0.6 centred at (0.5 + s/50, 0) in the lesson's vortex (Γ = 1, r_c0 = 1).

(a) Compute the turnover time 2πr/u_θ at the loop's centre.
(b) Advect the loop with `material_loop_exact` (4000 points) to 0, 1, 2 and 4 turnovers with ν = 0 and with ν = 0.005. For each, report the circulation from `circulation_polygon` and from `material_circulation_by_area`, and the loop's length.
(c) In the inviscid case the line-integral circulation drifts slightly. Show, with 1000, 2000 and 4000 points, that the drift is a discretisation error that falls as the square of the segment length, not a failure of Kelvin's theorem.
(d) Why does the viscous loop stretch less than the inviscid one?

Answer form: a time; a table of 8 rows × 3 numbers; three drifts and an order; a sentence.

## 5. (Numerical, uses your code: downwash at your tailplane)

The light aeroplane (750 kg, 10 m span, elliptic loading, 60 m/s at 3000 m) has its tailplane l_t = 3.5 + s/10 m behind the wing's bound vortex.

(a) Compute Γ₀ for the horseshoe model and the downwash and downwash angle at the tail on the centreline, in the plane of the wing and 0.5 m above it.
(b) Split the in-plane downwash into the bound vortex's and the trailing legs' contributions. Compare the total with the far-wake value −2Γ₀/(πb′).
(c) How does the downwash vary across the tailplane's span (y = 0, 1, 2 m at z = 0.5 m)? Why does it grow outboard?
(d) The tail is set at an incidence chosen for this downwash. What happens to the tail's angle of attack if the flaps are lowered and the wing's circulation rises by 30 % at the same speed?

Answer form: a circulation and four numbers; two contributions and a comparison; three velocities and a sentence; a sentence with a number.

## 6. (Judgment: how fast does lift build up?)

When the wing of Problem 5 starts impulsively at angle α = 2.733° (its cruise c_l), Kelvin's theorem requires a starting vortex of strength equal and opposite to the bound circulation Γ = πcV sin α.

(a) Using a 2-D point-vortex model (the starting vortex left at the point where the motion began, the wing a distance d downstream), find the vertical velocity the starting vortex induces at the wing as a fraction of the free-stream angle, for d = 1, 5, 10 and 50 chords.
(b) After how many chords of travel has the starting vortex's influence fallen to 10 % of α?
(c) The real build-up of lift is Wagner's problem, which follows the whole sheet of vorticity shed while the circulation grows, not a single point vortex shed at once. What does the point-vortex model get right, and what would you check before using it for a gust-response calculation?

Answer form: four fractions; a distance; two or three sentences.
