# P05.04  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` (and the probe scripts that call it) and given to the figures shown. Friction coefficients and the hinge-moment model are illustrative, not measured.

## 1. An aircraft parked on a ramp (s = 15: μ = 0.55)

**Assumptions.** Rigid aircraft, Coulomb friction at the braked tyres, one coefficient, no wind.

(a) Weight W; normal force W cos α, component down the slope W sin α. Braked wheels can supply up to μ W cos α along the slope, so the aircraft holds while W sin α ≤ μ W cos α, that is **tan α ≤ μ**: α_max = arctan 0.55 = **28.81°**.

(b) Only the mains brake: the friction available is μ × 0.9 W cos α, so tan α ≤ 0.9μ = 0.495: α_max = **26.34°**.

(c) Both the driving component and the friction are proportional to W, so W cancels. The 90 % share ignores the shift of normal load between nose and mains that the slope itself causes (the CG is above the ground, so the weight's component along the slope moves load toward the downhill wheels), and it ignores the rolling resistance of the free nose wheel.

**Check.** At μ = 0 the aircraft rolls on any slope; at 0.9μ the angle is smaller, as it must be.

**Common mistakes.** Using W instead of the normal force W cos α; applying friction at a free-rolling wheel.

## 2. Derivations

(a) **Capstan.** An element of cable subtending dθ on a fixed drum has tension T at one end and T + dT at the other. Their components toward the centre are each T sin(dθ/2), so the drum pushes the element with N = T dθ (first order). If the cable is about to slide toward the T + dT end, friction μN acts the other way: dT = μT dθ. Separate and integrate from T₁ at θ = 0 to T₂ at the full wrap: ln(T₂/T₁) = μθ, so **T₂/T₁ = e^{μθ}**. Check: μ = 0 gives T₂ = T₁ (a frictionless pulley); the ratio for two wraps is the square of the ratio for one, as it must be for elements in series.

(b) **Friction circle.** When a joint turns, the pin rolls up the bearing until friction is about to slip: at the contact the reaction leans by the friction angle φ_f = arctan μ from the normal, which passes through the centre. The reaction's line therefore passes at r sin φ_f from the centre: **r_f = r μ/√(1 + μ²)**, and the moment about the centre is r_f |R|, opposing the rotation. **Pulley.** Tight side T_t, slack side T_s, wrap θ. The pin reaction balances the vector sum of the tensions, of magnitude |R| = √(T_t² + T_s² − 2T_tT_s cos θ) (the angle between the two pulls is π − θ). Moments about the centre: (T_t − T_s)R = r_f |R|. Divide by T_s, write k = T_t/T_s and ρ = r_f/R, square: (k − 1)² = ρ²(k² + 1 − 2k cos θ), that is **(1 − ρ²)k² − 2(1 − ρ² cos θ)k + (1 − ρ²) = 0**; take the root above 1. At θ = π the tensions are parallel and |R| = T_t + T_s, so k − 1 = ρ(k + 1): **k = (1 + ρ)/(1 − ρ)** (and the quadratic gives the same root). If ρ ≥ 1 the pin's friction circle reaches the rim and no tension ratio can turn the pulley.

(c) **Virtual work.** The linkage has one degree of freedom. In a virtual rotation δφ of the stick the grip moves L_g δφ along the grip force's line (F is perpendicular to the stick), and the elevator turns by δδ = (dδ/dφ) δφ, against the hinge moment H. Pin reactions do no work (their points do not move) and the pushrods are internal. For a frictionless machine in equilibrium, δW = F L_g δφ + H δδ = 0 for every δφ, so **F = −H (dδ/dφ)/L_g**.

**Common mistakes.** In (a), putting the friction in the wrong direction (it gives e^{−μθ}, which is the ratio for sliding the other way); in (b), using the sum of the tension magnitudes for every wrap angle.

## 3. An over-centre toggle lock

**Assumptions.** Frictionless pins and slider, rigid links, quasi-static.

(a) The centre pin carries P and the two link forces. By symmetry the links carry the same force T (compression, pushing outward along the links). Perpendicular to the line: 2T sin θ = P, so **T = P/(2 sin θ)**. At the slider, the link pushes along itself; its component along the line is the output **Q = T cos θ = P/(2 tan θ)**; the component perpendicular is taken by the slider guide.

(b) With P = 20 N:

| θ | Q/P | Q (N) | link force (N) |
|---|---|---|---|
| 10° | 2.836 | 56.71 | 57.59 |
| 5° | 5.715 | 114.30 | 114.74 |
| 1° | 28.645 | 572.90 | 572.99 |

(c) The slider is at x = 2l cos θ, the centre pin at y = l sin θ. Virtual work: P δy = Q |δx|, so Q = P (dy/dθ)/|dx/dθ|. Central differences at θ = 5° with a step of 10⁻⁷ rad give **Q = 114.30 N**, equal to (b).

(d) As θ → 0 the force ratio grows without bound: a tiny actuator force holds a large load, and the reverse is also true: a load Q along the line produces only a perpendicular force 2Q tan θ on the centre pin, which vanishes in line. Past centre (θ slightly negative, against a stop) any load on the slider pushes the centre pin *further* against the stop, so the lock holds without the actuator. In the elevator linkage, a pushrod in line with an arm is the same toggle: the gearing goes to zero, the pilot's force needed for a given hinge moment goes to infinity (or the control locks), and small rigging errors give large force errors. Designs keep the arms near perpendicular to the pushrods over the travel.

**Common mistakes.** Writing Q = P/(2 sin θ) (that is the link force, not the output); forgetting that the guide takes the perpendicular component at the slider.

## 4. Change the gearing (s = 15: b₂ = 0.095 m)

(a) Small-angle gearing −(a_s/b₁)(b₂/h) = −(0.08/0.08)(0.095/0.08) = **−1.1875**; `gearing` gives −1.18750 at neutral.

(b) δ(+20°) = **−23.78°**, δ(−20°) = **+24.15°**; F(+20°) = **74.33 N**, against −25.19°, +25.44° and 84.19 N for the lesson's b₂ = 0.10 m. Less gearing: less elevator travel and a lighter stick, because F ∝ H × gearing and H itself is proportional to δ.

(c) `brentq` gives **b₂ = 0.09934 m**; then δ(−20°) = **+25.27°**, F(+20°) = **82.83 N**, F(−20°) = **−87.28 N**, neutral gearing **−1.2417**. The travel is asymmetric because the horn and the bellcrank's second arm swing through arcs: at neutral the second pushrod is perpendicular to both, and moving either way tilts it, but the geometry of the two circles is not symmetric about neutral (the second pushrod's end on the horn moves toward or away from the bellcrank pivot's level), so the elevator moves a little more on the pull side.

**Check.** At neutral, F ≈ k × gearing² × φ/L_g for small φ (since H = −kδ and δ ≈ gearing × φ): the force scales with the square of the gearing, (1.1875/1.25)² = 0.903, and 74.33/84.19 = 0.883, close (the rest is the nonlinearity at 20°).

**Common mistakes.** Changing b₂ without recomputing the second pushrod's length, which leaves the linkage not at neutral at φ = 0.

## 5. Pushrods against cables

F_ideal = **84.192 N** (frictionless linkage at φ = 20°).

(a) Plain-bushing pulley ratios: **1.0279** (90°), **1.0150** (45°), **1.0101** (30°); fairlead, 3° at μ = 0.2: **1.01053**.

(b)

| system | ratio k | driving (N) | back-driven (N) | band (N) |
|---|---|---|---|---|
| pushrods, pins μ = 0.2 | — | 87.31 | 81.19 | 6.12 |
| cable, plain-bushing pulleys | 1.1179 | 94.12 | 75.31 | 18.81 |
| cable, ball-bearing pulleys | 1.0347 | 87.11 | 81.37 | 5.75 |

The plain-bushing run is three times worse than the pushrods; with ball bearings the cable run is comparable to the pushrods. Of the ball-bearing run's 3.5 %, the three fairleads alone contribute 3.2 % (1.01053³ = 1.0319).

(c) The seized 90° pulley becomes a capstan with e^{0.3π/2} = 1.602, and the run's ratio becomes **1.656**: **139.4 N** to push and **50.8 N** to hold, a band of 88.6 N. The stick feels very heavy one way and the elevator will not return the stick toward neutral against a light hand: a control that feels broken, which is why pulleys are inspected for free rotation.

**Check.** A run's ratio is the product of its elements (they are in series); `cable_run_ratio` on a single element returns that element's ratio.

**Common mistakes.** Adding the percentages instead of multiplying the ratios (close for small values, wrong for the seized case); treating a pulley as a capstan when it turns.

## 6. Is friction a problem here?

(a) **Pulleys.** From Problem 5, nylon plain bushings (μ_b ≈ 0.15, illustrative) instead of ball bearings raise the run ratio from 1.035 to 1.118 and the band at full travel from 5.7 N to 18.8 N, more than three times. **Stick bolt.** The friction moment is proportional to the pin radius. Doubling the stick-pivot bolt alone (6.35 to 12.7 mm, μ = 0.2) changes the pushrod band from 6.12 N to **7.69 N** (efficiency 0.964 to 0.955); doubling all three pivots would give 12.25 N (efficiency 0.930). The premise is also wrong: the stick pivot's reaction (715.5 N) is not the largest; the elevator hinge carries 968.6 N, and the friction moments at full travel are 0.45, 0.50 and 0.60 N·m at stick, bellcrank and hinge.

(b) The band is proportional to the load, so it grows with q (four times at twice the airspeed; at the lesson's 50 m/s × 2 the pushrod band would be 24.5 N) while the efficiency is unchanged. Static friction exceeds running friction, so the breakout force exceeds the computed driving force. With friction the stick may not return to neutral when released at low speed, where the aerodynamic centring is weak; the pilot then has a dead band near neutral exactly where precise control matters. Plain bushings wear, change with temperature (nylon swells with moisture and stiffens in the cold), and collect dirt.

(c) A stick-force test on the installed system: apply known hinge-moment loads (or fly), move the stick slowly through the travel in both directions and record force against position; the loop's width is the band and its start is the breakout. Repeat cold and after a wear test. Compare against the applicable certification limits for control-system friction before accepting either change.

## Check-yourself answers

1. A member loaded at three or more points, or between its pins: it carries bending as well as axial force, so its pin forces are not along its axis and each pin has two unknown components; two-force members have one.
2. Newton's third law at the pin: whatever the pin exerts on one part, the part exerts back. A pin between two parts introduces two unknowns in a plane problem (one force vector), shared with opposite signs.
3. dδ/dφ ≈ −(a_s/b₁)(b₂/h): stick arm over first bellcrank arm, times second bellcrank arm over horn; the sign depends on the arrangement.
4. For a frictionless machine in equilibrium, the work of all external forces in any admissible virtual displacement is zero: F L_g δφ + H δδ = 0, so F = −H (dδ/dφ)/L_g.
5. The circle about the pin centre, radius r_f = rμ/√(1 + μ²), to which the bearing reaction is tangent when the joint turns; its moment r_f |R| opposes rotation. For r = 5 mm, μ = 0.2: **0.981 mm**.
6. Friction always opposes the motion: when the pilot drives the elevator against the air, friction adds to the force needed; when the air pushes the elevator back, friction helps the pilot hold. The difference between the driving and back-driven forces is the friction band (6.12 N in the lesson at full travel).
7. Element dθ: N = T dθ, dT = μN, so dT/T = μ dθ and T₂/T₁ = e^{μθ}. Over a full turn at μ = 0.2: e^{0.4π} = **3.514**.
8. A fairlead is a capstan: the cable slides, and the ratio grows exponentially with the wrap (a 90° fairlead at μ = 0.2 costs 37 %, a 90° plain pulley 2.8 %). A pulley turns, so only the pin friction, scaled down by r_f/R, acts. A seized pulley becomes a capstan, and the control becomes heavy one way and fails to centre.
