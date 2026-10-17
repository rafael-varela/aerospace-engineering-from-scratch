# P10.17  Worked solutions

Each solution states: assumptions, units, the derivation or method, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15; numerical results come from `solution.py` and probe scripts that call it, given to the figures shown. The blade geometry and section model are the lesson's stated assumptions, not measured data, and nothing here is a measurement.

## 1. The bound for a bigger UAV (s = 15: T = 2.25 N)

**Assumptions.** Actuator disk, uniform inflow, no losses; ρ = 1.225 kg/m³; A = π(0.30)²/4 = 0.07069 m².

(a) `actuator_disk(2.25, 15, A)`: **v = 0.821 m/s**, slipstream **16.64 m/s** (V + 2v), ideal power **35.60 W**, **η_i = 0.9481**.

(b) At 6000 rpm (n = 100 rev/s): C_T = 2.25/(1.225 × 100² × 0.30⁴) = **0.02268**, J = 15/(100 × 0.30) = **0.500**, and (5) gives **0.9481**, the same as (a) to round-off: (5) is (3) in other variables.

(c) **P = 10^{3/2}/√(2 × 1.225 × 0.07069) = 75.99 W**, v_h = √(10/(2 × 1.225 × 0.07069)) = **7.60 m/s**.

**Independent check.** Momentum: ṁ = ρA(V + v) = 1.370 kg/s, ṁ × 2v = 2.25 N ✓. Energy: ½ṁ((V + 2v)² − V²) = 35.60 W ✓.

**Common mistakes.** Using the diameter instead of the radius in A; forgetting that P = T(V + v), not TV.

## 2. Derivations

(a) Mass: ṁ = ρA(V + v). Momentum: T = ṁw. Energy: P = ½ṁ((V + w)² − V²) = ½ṁw(2V + w). Work at the disk: P = T(V + v) = ṁw(V + v). Equating, ½(2V + w) = V + v ⇒ **w = 2v**. Then T = 2ρA(V + v)v. With a = v/V, A = πD²/4 and V = JnD: T = (π/2)ρn²D⁴J²a(1 + a), so C_T = (π/2)J²a(1 + a), a = [−1 + √(1 + 8C_T/(πJ²))]/2, and **η_i = V/(V + v) = 1/(1 + a) = 2/(1 + √(1 + 8C_T/(πJ²)))**.

(b) Axial balance with W = (V + u_a)/sin φ: 4πrρF(V + u_a)u_a = ½ρ(V + u_a)²Bc C_n/sin²φ ⇒ u_a/(V + u_a) = Bc C_n/(8πrF sin²φ) = σ′C_n/(4F sin²φ) = **k**, so **u_a = Vk/(1 − k)** and V + u_a = V/(1 − k). Torque balance with W² = (V + u_a)(Ωr − u_t)/(sin φ cos φ) (product of the two expressions W = (V + u_a)/sin φ = (Ωr − u_t)/cos φ): 4πr²ρF(V + u_a)u_t = ½ρ(V + u_a)(Ωr − u_t)Bc C_t r/(sin φ cos φ) ⇒ u_t/(Ωr − u_t) = σ′C_t/(4F sin φ cos φ) = **k′**, so **u_t = Ωr k′/(1 + k′)** and Ωr − u_t = Ωr/(1 + k′). The triangle (V + u_a) cos φ = (Ωr − u_t) sin φ becomes V cos φ/(1 − k) = Ωr sin φ/(1 + k′); multiplying by (1 − k)(1 + k′)/(Ωr): **R(φ) = sin φ(1 − k) − (V/(Ωr)) cos φ(1 + k′) = 0**.

(c) With cd = 0, dT = Bc ½ρW² cl cos φ dr and dQ = Bc ½ρW² cl sin φ r dr, so dP = Ω dQ = dT Ωr tan φ = dT (V + u_a) (no swirl: tan φ = (V + u_a)/(Ωr)). Total: P = ∫(V + u_a) dT with dT = 4πρ(V + u_a)u_a r dr. Minimise P − λT over the function u_a(r): ∂/∂u_a[(V + u_a)²u_a − λ(V + u_a)u_a] = 0 ⇒ (V + u_a)(V + 3u_a) = λ(V + 2u_a), whose solution u_a is the same on every annulus (the equation does not contain r, and its left side divided by V + 2u_a increases monotonically with u_a, so the root is unique). Uniform inflow is therefore the only stationary point; it is the minimum, because a non-uniform slipstream carrying the same momentum flux always carries more kinetic energy. That proves the bound: no blade beats η_i at its own thrust.

**Common mistakes.** Taking W² = (V + u_a)² + (Ωr)² (forgetting the swirl) in the torque balance; dividing by V in the residual and then failing at V = 0.

## 3. A quadcopter in hover

**Assumptions.** Four isolated rotors (no interference), static BEM with the lesson's blade and section, ISA.

(a) T = 1.2 × 9.81/4 = **2.943 N** per rotor. `rpm_for_thrust(blade, 2.943, 0.0)` gives **5695 rpm**, and `bem` gives **21.03 W** of shaft power per rotor.

(b) Ideal power 2.943^{1.5}/√(2 × 1.225 × 0.05067) = **14.33 W**; **FM = 0.681**. Profile drag is the largest loss (as in the lesson's loss budget); larger chord at lower rpm, or better low-Re sections, would raise FM.

(c) Electrical power 4 × 21.03/0.8 = **105.2 W**; usable energy 11.1 × 2.2 × 3600 × 0.8 = 70 330 J; **11.1 min** of hover.

(d) Hover needs 4 × 21.03 = 84.1 W of shaft power against 21.4 W for winged cruise: **3.9 times** more. Momentum theory explains it: to hold 11.8 N the rotors accelerate the air through only 4 × 0.0507 m² at hover, and the ideal power T^{3/2}/√(2ρA) is large for a small A. In cruise the wing supports the weight by deflecting the air that flows past its whole 1.5 m span at 12 m/s, a very large mass flow, so its induced power is small, and the propeller only has to overcome drag, a tenth of the weight.

**Independent check.** The static FM from the lesson's Experiment 3 at 6000 and 8000 rpm is 0.681, the same: the static coefficients do not depend on rpm in this model.

**Common mistakes.** Using the 12 m/s cruise solve for hover; forgetting the motor efficiency.

## 4. Choosing the pitch

**Assumptions.** Same chord, constant geometric pitch, same section model.

(a) p = 0.152 m: **5158 rpm, J = 0.550, η = 0.753, P = 19.18 W**.

(b) Peak efficiency: 10 × 6, **0.753 at J ≈ 0.55** (it works at its peak in cruise); 10 × 4.5, **0.684 at J ≈ 0.42** (the cruise J = 0.457 is just past its peak).

(c) Static thrust at each propeller's own cruise rpm: 10 × 4.5 **3.49 N** at 6203 rpm, 10 × 6 **2.76 N** at 5158 rpm. Fit the 10 × 6 if cruise endurance matters most: it saves 10 % of shaft power. The price is take-off and climb thrust at the same motor rpm, and the higher torque per rpm the motor must deliver.

(d) For s = 15, p = 0.130 m: **5702 rpm, J = 0.497, η = 0.715, P = 20.20 W**, between the two.

**Independent check.** Each η is below `ideal_efficiency(C_T, J)` (the test `test_real_blade_never_beats_momentum_theory` is the general version).

**Common mistakes.** Changing the pitch in the docstring but not in `beta`; comparing propellers at the same rpm instead of the same thrust.

## 5. When does the classical iteration fail?

(a) Converged (✓) or not (✗) after at most 4000 iterations at the cruise point:

| n_el | 0.5 | 0.3 | 0.2 | 0.1 | 0.05 |
|---|---|---|---|---|---|
| 40 | ✓ | ✓ | ✓ | ✓ | ✓ |
| 80 | ✗ | ✓ | ✓ | ✓ | ✓ |
| 160 | ✗ | ✗ | ✓ | ✓ | ✓ |
| 320 | ✗ | ✗ | ✗ | ✓ | ✓ |
| 640 | ✗ | ✗ | ✗ | ✗ | ✓ |

The largest safe relaxation halves each time the resolution doubles: 0.5, 0.3, 0.2, 0.1, 0.05. Converged cases need 43–420 iterations.

(b) F in the last three of 160 elements: **0.0333, 0.0204, 0.0090**. In (8) the new u_a is roughly (dT/dr)/(4πrρF(V + u_a)): an error δ in the element's loads, caused by an error in φ, is amplified by 1/F, about 110 in the last element. Halving the elements near the tip pushes the last midpoint closer to R, halves F there, and doubles the amplification: the iteration's contraction factor exceeds one unless the relaxation shrinks in proportion.

(c) `bem` converges for every n_el (bisection cannot fail once R(φ) is bracketed). Thrust errors against n_el = 1280: 2.09 × 10⁻³, 6.51 × 10⁻⁴, 1.87 × 10⁻⁴, 5.04 × 10⁻⁵, 1.30 × 10⁻⁵ N for n_el = 10–160, observed orders **1.68, 1.80, 1.89, 1.95**.

**Common mistakes.** Blaming the physics for a solver failure; reporting a "converged" fixed point that hit max_iter.

## 6. What does the model not know?

(a)

| cd0 | rpm | η | P (W) |
|---|---|---|---|
| 0.010 | 6189 | 0.779 | 18.56 |
| 0.020 | 6203 | 0.676 | 21.37 |
| 0.030 | 6217 | 0.597 | 24.22 |

The rpm is robust (±0.2 %): it is set by the thrust, which depends on lift. The power and η are not (±13 %).

(b) (1) Laminar separation bubbles: they raise cd, sometimes greatly, below Re ≈ 10⁵, and lower η. (2) The Reynolds-number dependence of lift: a lower lift slope and an earlier, softer stall at low Re; at the cruise point the sections are far from stall, so the slope matters more, moving the rpm up for the same thrust. Static operation, with the inner sections at high α, is where stall matters.

(c) At 0.75R: **α = −1.22°, cl = 0.271, cl/cd = 13.5**. The section model's best cl/cd is **31.1 at cl = 0.94**: the blade works at less than half its best lift-to-drag ratio. A designer would reduce the blade area (narrower chord) or the rpm and increase pitch, so that each section works at a higher cl, as the 10 × 6 does.

(d) (1) Thrust and torque against J on a propeller test stand in a wind tunnel: tests the whole model, especially the section drag (through C_P) and α₀ (through the zero-thrust J). (2) A traverse of the slipstream behind the disk (axial and swirl velocity against r): tests the annulus-independence assumption and Prandtl's tip loss directly.

**Common mistakes.** Treating BEM's agreement with momentum theory as validation of the section data; it validates the code.

## Check-yourself answers

1. J = V/(nD) (flight distance per revolution in diameters); C_T = T/(ρn²D⁴), C_P = P/(ρn³D⁵), η = J C_T/C_P. C_P = 2πC_Q, because P = 2πnQ.
2. Half: v = w/2. Mass, momentum and energy conservation together (with the work done by the disk, P = T(V + v)).
3. At fixed thrust the kinetic energy left in the slipstream is least when its velocity increase is uniform (the calculus of variations in Problem 2(c)), and the Froude efficiency is exactly that uniform case; every real loss (non-uniform inflow, swirl, tip loss, profile drag) adds to the power.
4. The finite number of blades: air leaks round each tip, the load falls to zero there and an annulus near the tip receives less momentum than the disk would give it. F = 0 at the tip, F → 1 for infinitely many blades.
5. Doubling the thrust multiplies it by 2^{3/2} = 2.83. Doubling the diameter quadruples A and divides the power by 2.
6. The sections are cambered (α₀ = −4°): they still lift when the geometric angle of attack is zero, so the blade must advance further per turn before the lift vanishes. The aerodynamic pitch at 0.75R is 0.622 D.
7. The fixed-point update of u_a near the tip divides by a small F and amplifies errors (Problem 5(b)), so its contraction factor exceeds one unless the relaxation is small enough; bisection only needs the sign change of R(φ) between the bracket ends and halves the interval every step, whatever F is.
