# P05.12  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe script that calls it, given to the figures shown. Material values are typical, not design allowables.

## 1. A pitot boom (s = 15: L = 0.75 m)

**Assumptions.** Rigid root, Euler–Bernoulli, small deflection, loads downward (negative).

(a) I = π(25⁴ − 21⁴)/64 = **9628 mm⁴**; EI = 73.1 × 10⁹ × 9.628 × 10⁻⁹ = **703.8 N·m²**. Probe: P = −6 × 0.3 × 9.80665 = **−17.65 N**. Tube: A = 138.2 mm², q = −6 × 2780 × 138.2 × 10⁻⁶ × 9.80665 = **−23.64 N/m**.

(b) Tip deflection from the probe PL³/(3EI) = **−3.527 mm**; from the tube's weight qL⁴/(8EI) = **−1.328 mm**; total **−4.855 mm**, by superposition (the sum of two separate solutions). Tip slope PL²/(2EI) + qL³/(6EI) = **−0.540°**.

(c) Root moment PL + qL²/2 = **−19.89 N·m**; stress |M|(D/2)/I = **25.8 MPa**, a margin of more than 10 on 2024-T3 yield. The boom is sized by stiffness and vibration (a probe that moves 5 mm at 6 g also vibrates in flight), and by handling loads on the ground, not by bending strength.

**Independent check.** `fd_beam(0.75, 703.8, lambda x: q + 0*x, 400, point_loads=[(0.75, P)])` gives a tip deflection of −4.855 mm.

**Common mistakes.** Forgetting the load factor on the tube's own weight; using D instead of D/2 for c.

## 2. Derivations

(a) M(x) = q(L − x)²/2 from the outboard free body. EI w′ = −q(L − x)³/6 + C₁, EI w = q(L − x)⁴/24 + C₁x + C₂. w′(0) = 0 ⇒ C₁ = qL³/6; w(0) = 0 ⇒ C₂ = −qL⁴/24. Expanding gives **w = qx²(6L² − 4Lx + x²)/(24EI)**. From EI w⁗ = q: EI w = qx⁴/24 + c₃x³/6 + c₂x²/2 + c₁x + c₀; w(0) = w′(0) = 0 ⇒ c₀ = c₁ = 0; EI w‴(L) = qL + c₃ = 0 and EI w″(L) = qL²/2 + c₃L + c₂ = 0 ⇒ c₃ = −qL, c₂ = qL²/2; the same polynomial. The second route needed no reaction, only the free-end conditions.

(b) Superposition: the tip of the free cantilever rises qL⁴/(8EI); a downward R at the tip lowers it RL³/(3EI); the prop requires zero: **R = 3qL/8**. Integration with R unknown: M = q(L − x)²/2 − R(L − x); EI w′ = −q(L − x)³/6 + R(L − x)²/2 + C₁; EI w = q(L − x)⁴/24 − R(L − x)³/6 + C₁x + C₂. w′(0) = 0 ⇒ C₁ = qL³/6 − RL²/2; w(0) = 0 ⇒ C₂ = −qL⁴/24 + RL³/6; w(L) = 0 ⇒ C₁L + C₂ = 0 ⇒ qL⁴/8 = RL³/3 ⇒ R = 3qL/8 ✓. Superposition used only the compatibility condition at the prop, w(L) = 0; the clamp conditions were built into the two standard cases.

(c) A point load P at η lifts the tip by Pη²(3s − η)/(6EI) and rotates it by Pη²/(2EI). With dP = w(η)dη, w = (4L/πs)√(1 − η²/s²), t = η/s:
θ_tip = (4L/πs)(s³/2EI)∫₀¹ t²√(1 − t²) dt = (2Ls²/πEI)(π/16) = **L s²/(8EI)**;
w_tip = (4L/πs)(s⁴/6EI)∫₀¹ t²(3 − t)√(1 − t²) dt = (2Ls³/3πEI)(3π/16 − 2/15) = **(L s³/EI)(1/8 − 4/(45π))**.
Proof of the second integral: u = 1 − t², du = −2t dt, t² = 1 − u: ∫₀¹ t³√(1 − t²) dt = ½∫₀¹(1 − u)√u du = ½(2/3 − 2/5) = **2/15**. (The first follows from t = sin θ: ∫₀^{π/2} sin²θ cos²θ dθ = π/16.)

**Common mistakes.** Measuring the moment from the root instead of the free end; forgetting that the prop's force acts against the load.

## 3. A strut-braced wing (s = 15: a = 2.90 m)

(a) Unknowns: the root force, the root moment and the strut force (three); equations: vertical force and moment equilibrium (two). **One redundant**: one compatibility condition is needed (w = 0 at the strut).

(b) q = 13 974.5/5 = 2794.9 N/m. The free cantilever deflects at a by `cantilever_udl(a, 5, q, EI)`; a force R at a deflects it there by R a³/(3EI). Equating: **R = 12 096 N downward on the wing** (the strut is in tension), 86.6 % of the lift.

(c) Root moment q s²/2 − R a = 34 936.2 − 35 079.4 = **−143 N·m** (almost zero, slightly hogging) against **34 936 N·m** unbraced; moment at the strut q(s − a)²/2 = **6163 N·m**; tip deflection **37.9 mm** against **627.7 mm** unbraced (EI = 3.479 × 10⁵ N·m²). Between the root and the strut the wing hardly moves (largest |w| 2.74 mm). The strut turns a 35 kN·m root moment into a 6 kN·m moment at the strut and lets the spar be much lighter, at the price of the strut's drag.

(d) `fd_beam(5, 3.479e5, lambda x: q + 0*x, 1000, point_loads=[(2.90, -12096.4)])` gives w = 0.0000 mm at the strut node, a tip deflection of 37.93 mm and a root moment of −143 N·m.

(e) Pinned root: moments about the root, R a = q s²/2 ⇒ **R = 12 047 N**, and the root carries 1928 N, also downward. Within 0.4 % of the clamped value: the root fitting's moment capacity hardly matters once the strut is there, which is why real strut-braced wings use pinned root fittings.

**Common mistakes.** Treating the strut as carrying all the lift; forgetting that the overhang outboard of the strut still bends the spar.

## 4. The tip deflection of your own wing (s = 15)

m = 730 kg, semi-span 5.25 m, h = 0.170 m, b_f = 46 mm: I = **6.739 × 10⁶ mm⁴**, EI = **4.832 × 10⁵ N·m²**, L = 13 601.8 N.

(a) `elliptic_tip`: **393.90 mm (7.50 % of the semi-span)**, slope **5.56°**; integration 393.90 mm; FD (n = 400) 393.84 mm.

(b) Taper 0.6: **448.57 mm** by integration, 448.49 mm by FD (**8.54 %**).

(c) Factor **1.709** on the root EI (same taper). Deflection is proportional to 1/EI when the spanwise shape of EI is unchanged, so scaling EI by k scales the deflection by 1/k: k = 448.57/(0.05 × 5250) = 1.709. A cap sized for strength (MS = 0 at ultimate) is not stiff enough here: this wing would be sized by its stiffness requirement.

**Common mistakes.** Using the semi-span of the lesson; forgetting the taper in (c).

## 5. Testing the finite-difference solver

(a) Clamped–clamped: EI w = qx⁴/24 + c₃x³/6 + c₂x²/2 with w(L) = w′(L) = 0 ⇒ **w = qx²(L − x)²/(24EI)**. Maximum FD errors 5.21 × 10⁻⁷, 1.30 × 10⁻⁷, 3.26 × 10⁻⁸, 8.14 × 10⁻⁹ m for n = 20, 40, 80, 160: observed orders **2.00, 2.00, 2.00**.

(b) Load on a node (a = 0.3 m): errors 4.37 × 10⁻⁶ … 1.71 × 10⁻⁸ m, orders **2.00**. Load between nodes (a = 0.3333 m, snapped): errors 3.15 × 10⁻⁵, 1.40 × 10⁻⁵, 7.47 × 10⁻⁶, 3.55 × 10⁻⁶, 1.89 × 10⁻⁶ m, orders **1.17, 0.90, 1.08, 0.91**: first order, erratic. Snapping moves the load by up to h/2, an O(h) error in the load position and hence in w, which swamps the O(h²) truncation error; the erratic orders come from the snapping distance varying with n. Fix: split the load between the two neighbouring nodes in proportion to distance, or put a node under it.

(c) Under a unit load at 0.3 m, w(0.7) = **4.180156 × 10⁻⁷ m**; under a unit load at 0.7 m, w(0.3) = **4.180156 × 10⁻⁷ m**; relative difference 2.9 × 10⁻¹³. The discrete system is symmetric (it is a symmetric discretisation of a self-adjoint operator), so reciprocity holds to round-off, not just to O(h²).

**Common mistakes.** Comparing at points that are not nodes; reading erratic orders as a bug in the solver rather than in the load placement.

## 6. Why is the measured tip deflection different?

Prediction 242.8 mm (spar alone, rigid root); measured 0.85 × 242.8 = **206.4 mm**.

(a) Shear deformation: M₀/(G A_web) = 29 655/(26.9 × 10⁹ × 4.8 × 10⁻⁴) = **2.30 mm, 0.9 %**, and it *adds* deflection: it makes the discrepancy slightly larger, not smaller.

(b) The measured tip slope is about 3.06°, w′² = 2.9 × 10⁻³: the small-slope error is of that order, **under 0.5 %**. Not the explanation.

(c) EI factor = 242.8/206.4 = **1.18** ignoring shear, **1.19** including it. The wing's skins and stringers, which form a wing box with the spar caps, can easily supply 20 % more EI.

(d) With a root of k_θ = 10⁷ N·m/rad adding 14.8 mm, the bending part must be 206.4 − 2.3 − 14.8 = 189.3 mm: factor **1.28**.

(e) (1) An inclinometer or dial gauges at the root to measure the root rotation directly; (2) deflections at several stations, since a root rotation adds a straight line in y while extra EI scales the whole cubic-like shape; (3) strain gauges on the caps to measure the curvature (and hence EI) at a station directly.

**Common mistakes.** Assuming every correction reduces the deflection; attributing a 15 % discrepancy to "experimental error" without testing each explanation.

## Check-yourself answers

1. EI w″ = M: E is Young's modulus, I the second moment of area about the neutral axis, w the deflection (positive up), M the bending moment (positive sagging, compressing the top).
2. Clamped: w = 0, w′ = 0. Pinned: w = 0, M = 0. Free: M = 0, V = 0.
3. Doubling L multiplies it by 2⁴ = 16. Doubling the depth of a rectangle multiplies I by 8, so the deflection falls to 1/8.
4. Because the beam equation is linear in w and q (linear elasticity, small slopes, fixed supports). It fails with large deflections, plasticity, contact or gaps at supports, and loads that change direction with the deflection.
5. Remove the prop, replace it by an unknown force, and require the total deflection at the prop to be zero: compatibility supplies the extra equation.
6. The deflection at point B due to a unit load at A equals the deflection at A due to a unit load at B.
7. The truncation error of central differences is proportional to the fourth derivative of the solution, which is smooth for a uniform load; the elliptic load behaves like √(s − y) at the tip, its derivatives are infinite there, and the error falls only as h^1.5.
8. The tank's weight opposes the lift, so moving it outboard gives more relief: the root bending moment falls further (the relief is the tank weight times its distance from the root), and the tip deflection falls further (the relief grows as y²(3s − y) with the tank station y). The catch is a heavier load case on the ground, where the tank weight is not balanced by lift.
