# P06.12  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for **s = 15**; numerical results are from `solution.py` and the probe script that calls it, given to the figures shown. All results are outputs of models, not measurements.

## 1. A pendulum in phase space (s = 15: l = 1.25 m, a = 50°)

**Assumptions.** Point bob, rigid massless rod, no friction; q = θ, p = θ̇ (H per unit ml²).

(a) ω₀ = √(9.80665/1.25) = **2.80095 rad/s**; T₀ = **2.24323 s**; E = −ω₀² cos 50° = −5.0429 s⁻²; **k² = 0.178606** (= sin²25°).

(b) Elliptic integral: **T = 2.354909 s**. dA/dE: **2.354909 s** (agreement to ten digits, 2.3549086036 against 2.3549086036). Series: **2.350004 s**, 0.21 % low: the a²/16 term is the first of an infinite series and 50° is already large enough for the next term to matter.

(c) Phase area **6.43728 rad²/s**; largest p on the orbit **2.36747 rad/s** (at q = 0); largest p on the separatrix 2ω₀ = **5.60190 rad/s**, a bob speed of 2ω₀l = **7.002 m/s** at the bottom (the speed that just reaches the top, √(4gl)).

**Check.** k² = sin²(a/2) as it should; T/T₀ = 1.0498, between 1 and the separatrix's infinity.

**Common mistakes.** Passing k instead of m = k² to `ellipk`; using degrees inside cos.

## 2. The Hamiltonian toolkit

(a) dH = Σ(q̇ᵢdpᵢ + pᵢdq̇ᵢ) − Σ(∂L/∂qᵢ dqᵢ + ∂L/∂q̇ᵢ dq̇ᵢ). The definition **pᵢ = ∂L/∂q̇ᵢ** makes the dq̇ terms cancel, which is why H depends on (q, p) only. Lagrange's equation gives ∂L/∂qᵢ = d/dt(∂L/∂q̇ᵢ) = ṗᵢ, so dH = Σ(q̇ᵢdpᵢ − ṗᵢdqᵢ). Comparing with dH = Σ(∂H/∂qᵢ dqᵢ + ∂H/∂pᵢ dpᵢ): **q̇ᵢ = ∂H/∂pᵢ, ṗᵢ = −∂H/∂qᵢ**.

(b) On a closed orbit of energy E between turning points q₋ and q₊, p(q, E) = √(2(E − V(q))) on the upper branch, and A(E) = 2∫ p dq. Differentiating, the limits contribute nothing (p = 0 there), and ∂p/∂E = 1/p, so dA/dE = 2∫ dq/p. Since p = q̇, dq/p = dt, and 2∫_{q₋}^{q₊} dt is the time out and back: **dA/dE = T**.

(c) The phase velocity is v = (∂H/∂p, −∂H/∂q). div v = ∂²H/∂q∂p − ∂²H/∂p∂q = **0** (equal mixed partials). The rate of change of the area of a region R carried by the flow is dA/dt = ∫_R div v dA (the transport theorem), so **A is constant**. With damping, v = (p, −V′(q) − cp) and div v = ∂p/∂q + ∂(−V′ − cp)/∂p = −c, so dA/dt = −cA and **A(t) = A₀e^{−ct}**, for any region.

(d) Kick: J = [[1, 0], [h a′(q), 1]], **det 1**. Drift: J = [[1, h], [0, 1]], **det 1**. Explicit Euler (both from old values): J = [[1, h], [h a′(q), 1]], **det 1 − h²a′(q)**. Symplectic Euler is a kick followed by a drift (det 1·1 = **1**); Verlet is half kick, drift, half kick (**1**). So symplectic Euler and Verlet preserve area at any h; explicit Euler does not (for the pendulum, a′ = −ω₀² cos q and the factor is 1 + h²ω₀² cos q).

**Check.** The tests confirm (b) numerically at six energies and (d) by finite differences.

**Common mistakes.** In (b), forgetting that the boundary terms vanish because p = 0 at the turning points; in (d), treating the Euler step as a composition of a kick and a drift (it is not: both use old values).

## 3. A bead on a spinning hoop (s = 15: Ω = 7 rad/s)

**Assumptions.** Frictionless bead, hoop spun at exactly constant Ω by a motor, point mass.

(a) p = mR²θ̇ and **H = p²/(2mR²) − ½mR²Ω² sin²θ − mgR cos θ**. T + V = ½mR²(θ̇² + Ω² sin²θ) − mgR cos θ: the Ω² term has the **opposite sign**, so **H ≠ T + V**. H is **conserved** (L has no explicit t, so the energy function is a constant of motion; here T is not homogeneous quadratic in θ̇ because of the Ω² sin²θ part). T + V is **not** conserved: the motor that keeps Ω constant does work on the bead as it moves out and in. The simulation in (d) shows both.

(b) Hamilton: θ̇ = p/(mR²), ṗ = mR²Ω² sin θ cos θ − mgR sin θ = mR sin θ(RΩ² cos θ − g). Equilibria at p = 0 with sin θ = 0 (θ = 0, π) or **cos θ₀ = g/(RΩ²)**, which exists only when RΩ² > g. With V_eff = −½mR²Ω² sin²θ − mgR cos θ, V_eff″(0) = mR(g − RΩ²): the bottom is stable only for **Ω < Ω_c = √(g/R) = 4.4287 rad/s**.

(c) Ω = 7 rad/s > Ω_c. Bottom: **unstable**, growth rate √(Ω² − g/R) = **5.4210 s⁻¹**. Top (θ = π): unstable, √(Ω² + g/R) = **8.283 s⁻¹**. Off-axis: cos θ₀ = 9.80665/(0.5 × 49) = 0.4003, **θ₀ = ±66.40°**, stable, with V_eff″(θ₀)/(mR²) = Ω² sin²θ₀ and frequency **Ω sin θ₀ = 6.4148 rad/s** (by the formula and by the numerical second derivative, equal).

(d) Over 10 s from θ₀ + 0.2 rad: H per unit mass varies by **1.5 × 10⁻¹⁰ J/kg** (solver tolerance), T + V by **3.73 J/kg**; θ stays between **53.8° and 77.9°**, oscillating about θ₀.

**Check.** For Ω → 0, H → T + V and the bead is a pendulum with ω = √(g/R).

**Common mistakes.** Calling H "the energy" here; forgetting the second unstable branch at the top.

## 4. Your own blob (s = 15: centre (0.3, 0.8ω₀))

(a) Centre energy **−0.6353 ω₀²**; area π × 0.2 × 0.2ω₀ = **0.394784 rad²/s**. The highest boundary energy is about −0.455ω₀², so the blob is **entirely inside the separatrix** (libration only).

(b) Exact flow over 5 s: area change **−1.26 × 10⁻⁷**, which is the polygon's chord error (the boundary has lengthened), not physics; Liouville says zero.

(c) Relative area change against the exact flow (Euler against the initial area):

| h (s) | Euler (1 s) | RK4 (5 s) | Verlet (5 s) |
|---|---|---|---|
| 0.02 | +16.7 % | −1.28 × 10⁻⁷ | −2.3 × 10⁻¹⁰ |
| 0.04 | +35.0 % | −4.09 × 10⁻⁶ | −1.2 × 10⁻⁹ |
| 0.08 | +69.4 % | −1.28 × 10⁻⁴ | +4.2 × 10⁻¹⁰ |

(d) RK4's loss grows by 32 and 31 per doubling of h: **≈ h⁵ over a fixed time**, i.e. a per-step defect ≈ h⁶ (the (ωh)⁶/72 of the oscillator) times T/h steps. Verlet's column is round-off and polygon noise at 10⁻¹⁰, with no trend in h: the map preserves area exactly at every step size.

**Common mistakes.** Comparing a map's blob with the initial area instead of with the exact flow carried with the same points, which mixes the polygon error into the integrator error.

## 5. An orbit-determination error (s = 15: a = 7100 km, e = 0.015, δv = 5 mm/s)

(a) v_p = 7605.97 m/s. First order: δE = v_pδv, δa = 2a²δE/μ = **9.6191 m**; exactly from the perturbed energy: **9.6191 m** (9.619096 against 9.619080, a difference of 2 × 10⁻⁶ relative).

(b) (13): **92.03 m** after 1 orbit, **920.28 m** after 10.

(c) Simulated along-track separations: **92.028 m** and **920.282 m**, equal to (b) to 10⁻⁶ relative. A 5 mm/s error, smaller than many thruster calibration errors, moves the satellite a kilometre along its orbit in under 17 hours (period 5953.86 s).

(d) J_r = **3.76058 × 10⁷ m²/s** (small, because the orbit is nearly circular: J_r → 0 as e → 0); dJ_r/dE gives the period **5953.86 s** = 2π√(a³/μ).

**Common mistakes.** Using δv instead of δa in (13); forgetting that a speed error at periapsis keeps the satellite at periapsis, so the eccentricity factor stays √((1 + e)/(1 − e)).

## 6. Phase space in practice

(a) **Right:** the cloud's phase-space volume (for a Gaussian cloud, det Σ) is constant. **Wrong:** volume is not size; the cloud shears, and its extent along the orbit grows linearly without bound (335 m after one orbit, 6.7 km after 20 in the lesson), balanced by no growth (and some narrowing in a combination of coordinates) elsewhere.

(b) **Numerics.** The Kepler flow conserves volume, so a 30 % fall is the integrator's dissipation: RK45 is not symplectic and, at rtol 10⁻⁶, contracts phase space each orbit (as RK4 did in experiment 2, by more). Use a symplectic fixed-step method for statistics over many orbits, or tighten the tolerance until det Σ is constant to the needed level, and test by propagating forward and back.

(c) **Hardly.** Over 90 minutes the structure-preserving benefit is negligible; what matters is the trajectory's accuracy, which a high-order adaptive method at tight tolerance delivers best. Symplecticity pays off for long propagations and for statistics, not for one orbit.

(d) **Atmospheric drag.** It is not a Hamiltonian force: it removes energy, so a shrinks and the period shortens continually; the satellite runs *ahead* of a drag-free prediction, and the along-track error grows quadratically in time rather than linearly, at a rate set by the uncertain air density.

## Check-yourself answers

1. p = ∂L/∂q̇. For a pendulum angle, p = ml²θ̇ (an angular momentum about the pivot); for the polar angle of an orbit, p_θ = mr²θ̇, the angular momentum (per unit mass, h = r²θ̇).
2. q̇ = ∂H/∂p = p/m, ṗ = −∂H/∂q = −kq, so mq̈ = ṗ = −kq.
3. f is constant along every motion. Kepler: the angular momentum h (θ is cyclic) and the Laplace–Runge–Lenz vector (so the orbit closes); also H itself.
4. Closed curves around (0, 0): libration; open curves above and below: rotation; the curves p = ±2ω₀ cos(q/2) through (±π, 0): the separatrix, whose motion takes infinite time to reach the upright position.
5. T = dA/dE. For a harmonic oscillator the orbit is an ellipse of area 2π(E − E_min)/ω, so T = 2π/ω whatever the energy.
6. The phase flow of a Hamiltonian system has zero divergence, so it preserves area (volume). The pendulum's period grows with energy, so points of the blob at higher energy lag behind those at lower energy: the blob is sheared into a filament, keeping its area.
7. Explicit Euler: 1 + h²ω₀² cos q per step (> 1 near the bottom), so energy and area grow. Verlet: exactly 1, so area is preserved and the energy error stays bounded (it oscillates about a modified energy).
8. A small error in energy (semi-major axis) changes the period, so the phase along the orbit drifts at a constant rate and the along-track error grows linearly; the energy error itself is conserved, so the orbit's size and shape (radial extent) stay within a fixed band.
