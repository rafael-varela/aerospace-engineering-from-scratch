# P06.05  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe scripts that call it, given to the figures shown. Dimensions are illustrative.

## 1. Classify, assemble, move (s = 15)

**Assumptions.** Rigid links, ideal pins, angles counterclockwise from O₂→O₄.

(a) (100, 40, 120, 80): 40 + 120 < 100 + 80, shortest is the crank: **crank-rocker**. (60, 100, 90, 110): 60 + 110 < 100 + 90, shortest is the ground: **double-crank**. (100, 70, 30, 90): 30 + 100 < 70 + 90, shortest is the coupler: **double-rocker**. (100, 80, 90, 95): 80 + 100 < 90 + 95, shortest is the crank: **crank-rocker**. (65, 20, 60, 40): 20 + 65 < 60 + 40, shortest is the crank: **crank-rocker**.

(b) θ₂ = 150°. Branch +1: **θ₃ = 27.25°, θ₄ = 110.46°**; branch −1: **θ₃ = −44.15°, θ₄ = −127.36°**. The transmission angle is **83.20° on both** (the branches are mirror images about the line A→O₄, so the angle between coupler and rocker is the same).

(c) Branch +1, ω₂ = 10 rad/s: **ω₃ = 2.137 rad/s, ω₄ = 4.235 rad/s**; |v_A| = ω₂r₂ = **0.400 m/s**, |v_B| = ω₄r₄ = **0.339 m/s**. Check: A = (−34.64, 20.00) mm, so v_A = 10 × (−0.0200, −0.03464) = (−0.200, −0.346) m/s; r_BA = (106.68, 54.95) mm, ω₃ × r_BA = 2.137 × (−0.05495, 0.10668) = (−0.1174, 0.2280) m/s; v_A + ω₃ × r_BA = **(−0.3174, −0.1184) m/s**, equal to ω₄ × r_B/O₄ = (−0.3174, −0.1184) m/s ✓.

**Common mistakes.** Measuring θ₄ from O₄→O₂ instead of +x; forgetting that ω × r turns r by +90°.

## 2. Derivations

(a) Real part: r₂ cos θ₂ + r₃ cos θ₃ − r₄ cos θ₄ − r₁ = 0; imaginary part: r₂ sin θ₂ + r₃ sin θ₃ − r₄ sin θ₄ = 0. ∂/∂θ₃ = (−r₃ sin θ₃, r₃ cos θ₃), ∂/∂θ₄ = (r₄ sin θ₄, −r₄ cos θ₄), giving (4). det J = r₃r₄ sin θ₃ cos θ₄ − r₃r₄ cos θ₃ sin θ₄ = **r₃r₄ sin(θ₃ − θ₄)**. It is zero when the coupler and the rocker are parallel or anti-parallel: B, A and O₄ in line (toggle). A small change of θ₂ then cannot be absorbed by small changes of θ₃ and θ₄.

(b) d/dt: r₂iω₂e^{iθ₂} + r₃iω₃e^{iθ₃} − r₄iω₄e^{iθ₄} = 0. Real and imaginary parts: −r₃ sin θ₃ ω₃ + r₄ sin θ₄ ω₄ = r₂ω₂ sin θ₂ and r₃ cos θ₃ ω₃ − r₄ cos θ₄ ω₄ = −r₂ω₂ cos θ₂: (5). d²/dt²: Σ ± r(iα − ω²)e^{iθ} = 0, so J[α₃, α₄]ᵀ = −K with
K_x = −r₂α₂ sin θ₂ − r₂ω₂² cos θ₂ − r₃ω₃² cos θ₃ + r₄ω₄² cos θ₄,
K_y = r₂α₂ cos θ₂ − r₂ω₂² sin θ₂ − r₃ω₃² sin θ₃ + r₄ω₄² sin θ₄.

(c) Bodies j, k, l in plane motion relative to j. A point of body k moves (relative to j) perpendicular to its line to I_jk; as a point of l, perpendicular to its line to I_jl. At I_kl the two velocities must be equal, so their directions must agree: the lines I_jk–I_kl and I_jl–I_kl must coincide, i.e. the three centres are collinear. For the four-bar, I₂₄ lies on line I₁₂–I₁₄ (the ground line) and on line I₂₃–I₃₄ (AB). As a point of 2 its velocity is ω₂ × (I₂₄ − O₂); as a point of 4, ω₄ × (I₂₄ − O₄). On the x-axis these are ω₂x₂₄ and ω₄(x₂₄ − r₁) (times e_y), so **ω₄/ω₂ = x₂₄/(x₂₄ − r₁)**.

(d) With r₁ = r₃ and r₂ = r₄, θ₃ = 0 and θ₄ = θ₂ satisfy (3): r₂e^{iθ₂} + r₁ − r₂e^{iθ₂} − r₁ = 0. Substituting in (5): J = [[0, r₂ sin θ₂], [r₁, −r₂ cos θ₂]] and the right side (r₂ω₂ sin θ₂, −r₂ω₂ cos θ₂); the first row gives ω₄ = ω₂, then the second r₁ω₃ = 0, **ω₃ = 0**. The coupler translates.

## 3. A flapping-wing drive

**Assumptions.** Constant crank speed 2π × 5 rad/s; ideal pins; branch with the lever above the ground line.

(a) s + l = 15 + 60 = 75 < 25 + 55 = 80, shortest is the crank: a **crank-rocker**; the motor can turn fully.

(b) Extended (O₂B = 75 mm) and folded (O₂B = 45 mm) triangles with sides 55 and 25 mm: rocker angles **θ₄ = 44.10°** (extended) and **126.22°** (folded); amplitude **82.13°**. A 0.1° sweep over a full crank turn gives 44.0954° to 126.2215°, the same to 10⁻⁵ degrees.

(c) Crank at the extended extreme **13.41°**, at the folded **206.63°** (−153.37°). From extended to folded the crank turns **193.2°**; back again **166.8°**: a **time ratio of 1.158**, so one stroke takes 16 % longer than the other (a "quick return").

(d) At 5 Hz: largest |ω₄| = **27.6 rad/s**, largest |α₄| = **1630 rad/s²**, smallest transmission angle **28.96°**. The transmission angle falls below the 40° rule of thumb: the drive works, but it is close to toggle at one end of the stroke and will be sensitive to pin clearance there; a longer crank-to-rocker ratio or a different ground length would raise it.

**Common mistakes.** Taking the extremes at the crank's horizontal positions; forgetting the folded case is r₃ − r₂, not r₂ − r₃.

## 4. Your own servo arm (s = 15: r₂ = 15 mm)

(a) 15 + 120 < 25 + 120: **crank-rocker**. Neutral: pushrod **4.779°**, horn **90.956°**.

(b) Deflection **−24.47° to +25.62°**; gearing **0.5992** at neutral (r₂/r₄ = 0.600), range **0.4344 to 0.6001**; best straight line **0.5738 °/°**, worst departure **1.18°**; minimum transmission angle **60.60°** at −45°. Compared with the 20 mm arm: less travel, a lower and less constant gearing, more transmission margin.

(c) At servo +20°, ω₂ = 2 rad/s, α₂ = 3 rad/s²: ω₄ = 1.16992 rad/s, α₄ = 1.32031 rad/s². Errors of central differences:

| h (ms) | \|ω₄ FD − ω₄\| (rad/s) | \|α₄ FD − α₄\| (rad/s²) |
|---|---|---|
| 10 | 8.59 × 10⁻⁵ | 2.85 × 10⁻⁴ |
| 5 | 2.15 × 10⁻⁵ | 7.13 × 10⁻⁵ |
| 2.5 | 5.37 × 10⁻⁶ | 1.78 × 10⁻⁵ |

Observed orders **2.00, 2.00** for both.

## 5. Choosing the horn

(a) Over ±45° (servo arm 20 mm):

| horn (mm) | travel ≥ ±25° | gearing max/min ≤ 1.25 | μ_min > 40° | torque 0.3 × max gearing ≤ 0.35 N·m |
|---|---|---|---|---|
| 15 | yes (−67.8/+76.9) | no (4.90) | no (7.5°) | no (1.96) |
| 20 | yes (±45.0) | yes (1.00) | yes (45.0°) | yes (0.300) |
| 25 | yes (−34.2/+34.7) | yes (1.20) | yes (52.9°) | yes (0.240) |
| 30 | yes (−27.5/+28.6) | no (1.34) | yes (57.3°) | yes (0.200) |
| 35 | no (−22.8/+24.6) | no (1.47) | yes (59.9°) | yes (0.172) |
| 40 | no (−19.4/+21.6) | no (1.60) | yes (61.4°) | yes (0.150) |

Only the **20 mm and 25 mm horns** meet all four.

(b) Free play at the surface = 1° × the neutral gearing: **1.00°** (20 mm) and **0.80°** (25 mm); 0.67° for 30 mm. A longer horn reduces free play and servo torque: backlash pushes toward the longer acceptable horn.

(c) The 20 mm horn makes a parallelogram, a change-point linkage: it can pass through a flat, folded configuration into its crossed branch, where the surface would move the wrong way. It also gives no torque advantage: the servo carries the full 0.3 N·m hinge moment, within 15 % of its 0.35 N·m rating.

(d) Choose the **25 mm horn**: it meets travel (±34°), linearity (gearing within 20 %), transmission (above 52°) and torque (0.24 N·m against 0.35 N·m) with margin. It cuts servo backlash at the surface to 0.8°, and it is a crank-rocker far from any folding position. Set the flight controller's gain with the neutral gearing 0.80 and accept the 0.83° nonlinearity, or correct it in software.

## Check-yourself answers

1. v_B = v_A + ω × r_B/A. The relative motion of B about A is a rotation, so it is perpendicular to r_B/A (the distance AB cannot change in a rigid body).
2. The point with zero velocity at this instant; it is where the perpendiculars to two points' velocity directions meet.
3. The three instant centres of any three bodies in plane motion lie on one straight line.
4. The Jacobian's determinant r₃r₄ sin(θ₃ − θ₄) vanishes when coupler and rocker are in line: the linkage cannot transmit motion, and an infinitesimal crank motion would need an infinite rocker rate.
5. s + l = 30 + 110 = 140 < 100 + 60 = 160 and the shortest is the crank: crank-rocker.
6. At a known position, the loop's derivatives are linear in the rates (ω and α appear to the first power; the ω² terms are known once the velocities are found).
7. ω₃ = 0 and ω₄ = 3 rad/s.
8. 0.5 × 0.6 = 0.30 N·m.
