# P10.20  Worked solutions

Each solution states assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15. Lattice numbers come from `solution.py` and a probe script that calls it (flat rectangular lattice, 4 chordwise × 40 spanwise panels with cosine spanwise spacing, α = 5° from the zero-lift line, unless stated). "Wieselsberger" is the classical empirical fit quoted in the lesson; it is not a model of this wing.

## 1. One horseshoe and its image

(a) Looking downstream: the real trailing vortices at (±b_v/2, h) rotate so that the air descends between them and rises outboard of them; the images at (±b_v/2, −h) rotate the opposite way, so at the wing's height between the tips they induce an *upwash*.

(b) Real pair at the centre: each vortex at distance b_v/2 gives Γ/(2π b_v/2) downward; total w₀ = 2Γ/(π b_v). Each image at distance r = √((b_v/2)² + (2h)²) gives Γ/(2πr), of which the vertical part is the fraction (b_v/2)/r, upward. So w/w₀ = 1 − (b_v/2)²/r² = (2h)²/((b_v/2)² + (2h)²) = **(4h/b_v)²/(1 + (4h/b_v)²)**.

| h/b_v | 0.05 | 0.1 | 0.25 | 0.5 | 1 |
|---|---|---|---|---|---|
| ratio | 0.0385 | 0.1379 | 0.5000 | 0.8000 | 0.9412 |

(c) With b_v = πb/4, h/b_v = 0.1/(π/4) = 0.127: ratio **0.206**. Wieselsberger's fit gives **0.515**; the lattice in the lesson gives 0.531. Wieselsberger is far closer. One horseshoe concentrates the whole load in two line vortices whose images sit directly below them; the induced velocity at the centre is dominated by them. A distributed load sheds vorticity across the span, the images are spread too, and the centre-point downwash is a poor measure of the drag, which integrates Γ·w over the span.

**Check.** h → ∞: ratio → 1; h → 0: ratio → 0 (the images cancel the trailing pair), as the tests check for the code.

**Common mistakes.** Putting the images at depth h instead of 2h below the wake; evaluating the downwash at one point and calling it the drag ratio.

## 2. Why images work

(a) A 2-D vortex along +x of strength Γ at (y₀, z₀) induces w = Γ(y − y₀)/(2π r²), v = −Γ(z − z₀)/(2π r²). On z = 0 the real vortex at z₀ = h and the image at z₀ = −h have the same r and the same y − y₀; with strengths Γ and −Γ their w cancel exactly, and their v add (the ground sees a tangential velocity, which is allowed). For a 3-D filament use Biot–Savart, dV = (Γ/4π) dl × r/|r|³, and the reflection matrix M = diag(1, 1, −1). The image element is M dl at the mirrored position, so for a point P on the plane (where MP = P) its r′ = M r. For any vectors, (Ma) × (Mb) = det(M) M(a × b) = −M(a × b); with the image circulation −Γ the image element contributes dV′ = M dV, which has the same horizontal components and the opposite vertical component. The vertical parts cancel element by element, for any orientation of the filament, so reflecting every point of a horseshoe and reversing Γ satisfies w = 0 on z = 0 exactly; the test finds a residual of 10⁻¹⁵ relative.

(b) In the Trefftz plane the image trailing vortices sit 2h below and rotate the other way: between them they induce upwash on the real wake. At a given Γ(y), w(y) is smaller, so Di = (ρ/2)∫Γ w dy is smaller. The lift (ρV∫Γ dy) is unchanged at given Γ, so the induced drag at a given lift falls: K < 1.

(c) Tangency at control point i: Σ_j w_ij Γ_j + V sin α = 0. The images add upwash, which reduces the magnitude of the induced downwash per unit circulation; to cancel the same V sin α the circulations must grow, so CL at fixed α rises. At fixed circulation the lift (Kutta–Joukowski) would not change; what the ground changes is the circulation a given angle produces.

**Common mistakes.** Reflecting the geometry but keeping Γ (doubles the downwash instead of cancelling it).

## 3. An airliner's flare (60 × 7 m flat rectangle)

(a)

| h [m] | h/b | CL/CL_free | K | Wieselsberger |
|---|---|---|---|---|
| 5 | 0.083 | 1.200 | 0.484 | 0.466 |
| 10 | 0.167 | 1.086 | 0.673 | 0.658 |
| 20 | 0.333 | 1.034 | 0.840 | 0.841 (outside its stated range) |

(b) At constant speed the lift coefficient is constant. Descending from 20 m to 5 m, the induced drag at that CL falls to 0.484/0.840 = 58 % of its value at 20 m, and the angle of attack needed falls because CL_α rises by about 16 %. With the thrust unchanged the aircraft decelerates less and needs less nose-up attitude: it tends to keep flying a few metres above the runway. That is the "float"; pilots reduce thrust to idle and hold the attitude to let the speed decay.

(c) (1) Deployed flaps load the inboard wing and put the trailing edge closer to the ground: chord-dominated ground effect grows, so the lift increase is larger than the flat-plate value at the same h (and the flap's own pressure field can reduce it, see Problem 6). (2) Sweep and dihedral put the tips higher than the root (a 5 m root height is perhaps 7–8 m at the tips): effective h is larger and the effect smaller. The fuselage and engines also block the flow under the wing.

**Common mistakes.** Using Wieselsberger outside its range without saying so; quoting CL ratio at fixed α as the lift change at fixed speed.

## 4. Your own flare (s = 15: h = 1.25 m, h/b = 0.125)

(a) CL/CL_free at fixed α = **1.126**; K = **0.593**. CL = 750 × 9.81/(½ × 1.225 × 30² × 12) = **1.112**. Induced drag in free air **314 N**, at h = 1.25 m **186 N**.

(b) α (from zero lift) for that CL: **13.60° → 12.08°**, 1.5° less.

(c) ΔD = 128 N, so the deceleration is **0.171 m/s² smaller**. Over a float of a few seconds that is a few tenths of a m/s of speed kept, and tens of metres of extra float distance if the pilot waits for the speed to decay: approach speed control matters.

(d) K = **0.586, 0.593, 0.596** for ns = 20, 40, 80. Converging at roughly first order; the spread (0.010) is smaller than the difference from Wieselsberger's fit and much smaller than the effect.

**Independent check.** Wieselsberger at h/b = 0.125 gives 0.577, within 0.015.

**Common mistakes.** Comparing induced drag at fixed α (where CL also changes) instead of at fixed lift.

## 5. The tail in ground effect

(a) Free air, tail at 4.8 m from the leading edge: dε/dα = **0.318** (tail in the wing plane) and **0.298** (0.5 m above). The far-field formula gives **0.358**. Even 100 m downstream the lattice gives 0.273: the formula assumes an elliptic load with uniform downwash, while a rectangular wing sheds strong vorticity near its tips and its wake downwash at the centre (where the tail is) is below the span average. Roll-up of the wake, not modelled, moves the tip vortices inboard and changes it again.

(b) With the wing at h and the tail at the same height: dε/dα = **0.124, 0.208, 0.293** for h = 1, 2, 5 m (0.126, 0.201, 0.275 with the tail 0.5 m higher).

(c) Tail lift slope (2.83 × 0.707 m rectangle): a_t = **3.73 /rad**; V_H = 2.0 × 4.5/(12 × 1.2) = **0.625**. At α = 12°: Δε = (0.318 − 0.124) × 12° = **2.34°** less downwash at the tail. The tail's angle of attack rises by 2.34°, its lift rises, and ΔCm = −V_H a_t Δε = −0.625 × 3.73 × 0.0408 = **−0.095**: **nose-down**. The pilot must pull the elevator further back in the last metres of the flare; this is why aircraft with marginal elevator authority at forward centre of gravity can fail to flare.

**Common mistakes.** Sign: less downwash means *more* tail lift and a nose-down moment.

## 6. Where the image model fails

(a) A thick or flapped wing near the ground forms a converging-diverging channel between its lower surface and the ground. The flow accelerates in the throat and the pressure on the lower surface can fall (a venturi effect), reducing lift. A flat lattice has no thickness, so it only sees the "ram" (chord) and span effects that raise lift. Recant's tests of flapped wings near a ground board are the kind of data that show it.

(b) A stationary ground board grows a boundary layer that thickens the effective ground and can separate under the wing's pressure field; the measured effect is distorted. Tunnels use moving belts running at the free-stream speed, boundary-layer suction ahead of the model, or a mirror-image model (the image method made physical).

(c) Descending, the wake left behind was shed at greater heights, so its images are farther away than the steady model assumes at the current height. The induced-drag reduction and lift increase lag behind the steady values: the steady result is an upper bound on the cushioning during a descent.

(d) The flat, stationary ground. Waves move the image system and change h by ±0.5 m (h/c changing by a large fraction); the lift and pitching moment of a WIG craft change strongly with h, and a wave crest can suddenly increase lift and pitch the craft up (or a trough drop it). Pitch stability with height is the classic WIG design problem.

## Check-yourself answers

1. The ground is a plane the flow cannot cross. A mirror-image wing at the same distance below it, with circulation of opposite sign, cancels the vertical velocity on the plane everywhere, so solving wing + image enforces the ground condition exactly.
2. The images induce upwash on the real wing and wake: the downwash per unit circulation falls, so (i) more circulation is needed to meet tangency at a fixed α (lift rises), and (ii) the Trefftz-plane drag at a given lift falls.
3. About half: K = 0.531 for the trainer's rectangular wing (lattice), 0.515 from Wieselsberger's fit.
4. Recover the free-air result at large h; check the span load is symmetric; compute w on the ground plane (residual 10⁻¹⁵).
5. 2CL_α/(πAR), twice the downwash angle at an elliptic wing: in the far wake of an elliptic wing the downwash is uniform and double its value at the wing. It overestimates for the trainer's rectangular wing at its tail (0.318 against 0.358) because the centre downwash of a rectangular wake is below average.
6. The ground reduces the downwash at the tail; the tail's angle of attack rises; the tail lifts more and the aircraft pitches nose-down; more up-elevator is needed in the flare.
7. Upper bound: in a descent the wake was shed higher up, so the images are farther away than the steady model assumes.
