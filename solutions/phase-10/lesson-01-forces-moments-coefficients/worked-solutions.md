# P10.01  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15. Numbers come from `solution.py` and the probe script that calls it, given to the figures shown. Every balance reading here is **SYNTHETIC**: invented to teach the reduction, not measured.

## 1. The trainer, heavier and faster (s = 15: m = 760 kg, V = 57.5 m/s)

**Assumptions.** ISA sea level; level flight, so lift equals weight; the wing's C_L is referred to S = 12.0 m².

(a) W = 760 × 9.81 = 7455.6 N. q = ½ × 1.225 × 57.5² = **2025.1 Pa**. C_L = 7455.6/(2025.1 × 12.0) = **0.3068**. Re = 57.5 × 1.2/1.461 × 10⁻⁵ = **4.72 × 10⁶**. M = 57.5/340.3 = **0.169**.

(b) Room air has ν = 1.81 × 10⁻⁵/1.204 = 1.503 × 10⁻⁵ m²/s. Matching Re with a 0.24 m chord needs V = 4.72 × 10⁶ × 1.503 × 10⁻⁵/0.24 = **295.8 m/s**, which is **M = 0.862**. A low-speed tunnel cannot do it, and if it could, the flow would be transonic, with shock waves the aeroplane never sees.

(c) Match Re where the boundary layer matters (drag, maximum lift, stall), because the flow at M = 0.17 is nearly incompressible and its small compressibility effect can be corrected (Phase 11). For lift slope and pitching moment at low angles, neither matters much, and the model's Re mismatch is usually accepted, with a trip strip to make the boundary layer turbulent where it would be on the aeroplane.

**Check.** Re and the model speed scale exactly as the scale factor: 295.8 × 0.24 = 71.0 = 57.5 × 1.2 × (1.503/1.461). Correct.

**Common mistakes.** Using ISA ν for the room-air tunnel (a 3 % error); forgetting that matching Re at small scale drives the Mach number up.

## 2. Axes, moments and the aerodynamic centre

(a) In chord axes the stream direction is ê_D = (cos α, sin α) (x aft, z up) and the lift direction, 90° anticlockwise from it in this view, is ê_L = (−sin α, cos α). The force is F = (A, N). Projecting, D = F · ê_D = A cos α + N sin α, and L = F · ê_L = −A sin α + N cos α. The two unit vectors are orthonormal, so L² + D² = A² + N²: a rotation. The inverse rotation by −α gives N = L cos α + D sin α and A = −L sin α + D cos α.

(b) A force F = (F_x, F_z) at r has the moment (r − r_P) × F about P, whose y-component is (z − z_P)F_x − (x − x_P)F_z. With x aft and z up, a positive y-component tilts the nose (at small x) upward, so this is the nose-up moment. With the forces N, A acting through O plus a couple M_O, M_P = M_O + (z_O − z_P) A − (x_O − x_P) N = **M_O + (x_P − x_O) N − (z_P − z_O) A**. Check: N > 0 acting behind P (x_O > x_P) gives −(x_O − x_P)N < 0, nose-down about P, as a picture shows.

(c) If M_ac does not change with α, then M_ref = M_ac + (x_ref − x_ac) N with both points on the chord line. Divide by qSc: **c_m,ref = c_m,ac + (x_ref − x_ac)/c · c_n**, so dc_m,ref/dc_n = (x_ref − x_ac)/c and x_ac = x_ref − c dc_m,ref/dc_n. The lever arm multiplies the normal force, so c_n is the exact abscissa; c_l = c_n cos α − c_a sin α differs from c_n by a fraction of order α² (and α c_d/c_l), which is why the textbook "slope against c_l" is an approximation that is good at low angles.

(d) At fixed N and A, ∂D/∂α = N cos α − A sin α = L. So δc_d = c_l δα, and a standard uncertainty u_α in α gives **u(c_d) ≈ c_l u_α**. For c_l = 0.64 and u_α = 0.05° = 8.73 × 10⁻⁴ rad: **u(c_d) = 5.6 × 10⁻⁴**, about 8 % of a typical section drag coefficient of 0.007. The same reasoning gives ∂L/∂α = −D: angle errors barely touch the lift.

**Common mistakes.** Rotating the wrong way (the drag of a pure normal force at positive α is +N sin α, not −); measuring moment arms from the wrong point; forgetting that α must be in radians in u(c_d) ≈ c_l u_α.

## 3. A balance mounted at an angle

**Assumptions.** Rigid sting; the readings are already tare-corrected; the offsets are measured along the balance axes. S = 0.25 × 0.50 = 0.125 m²; qS = 75.0 N; qSc = 18.75 N·m.

(a) The balance axes are at α_b = 6.0 − 2.0 = 4.0° to the stream, so L = 74.9 cos 4° + 4.33 sin 4° = **75.02 N** and D = 74.9 sin 4° − 4.33 cos 4° = **0.905 N**: **c_l = 1.0003, c_d = 0.01207**. The quarter chord is 50 mm *ahead of* and 30 mm *above* the balance centre in balance axes, so M_c/4 = 2.68 + (−0.050)(74.9) − (0.030)(−4.33) = **−0.9351 N·m**, **c_m,c/4 = −0.0499**.

(b) Rotating by 6.0° instead: c_l = 0.9992 (−0.1 %) but **c_d = 0.0470, nearly four times too large**: D picks up N sin 2° ≈ 2.6 N, three times the real drag. Forgetting the 30 mm offset: c_m,c/4 = **−0.0568** instead of −0.0499, a 14 % error, from the axial force acting on a 30 mm lever. The mount angle spoils drag; the vertical offset spoils moment; neither noticeably spoils lift.

(c) The axial force is the component along the chord. At α_b = 4° the lift vector, tilted forward by 4° relative to the balance normal, has a forward component L sin 4° = 5.2 N, larger than the backward component of drag (0.90 N). On the real airfoil this forward pull is leading-edge suction: low pressure on the nose pulls the section forward.

**Check.** The inverse rotation of (L, D) = (75.02, 0.905) by 4° returns (74.90, −4.33) to rounding.

**Common mistakes.** Using α (chord) for a balance that does not share the chord's axes; measuring the offset in chord axes but applying it in balance axes; reading the negative A as an error.

## 4. Your own synthetic balance run (s = 15)

**Assumptions.** The data are **SYNTHETIC** (`synthetic_balance_data(seed=15)`); the fit assumes independent Gaussian scatter in c_l and c_m and no error in α or c_n, which is not quite the noise model, but nearly.

(a) Lift line: **a₀ = 6.007 ± 0.018 /rad**, **α_L0 = −2.092 ± 0.020°** (the intercept's uncertainty propagated through −intercept/slope with the full covariance). Moment line against c_n: **x_ac/c = 0.24805 ± 0.00015**, **c_m,ac = −0.05000 ± 0.00011**.

(b) Over seeds 1–20 the standard deviation of x_ac/c is **1.70 × 10⁻⁴** and the mean fit uncertainty is **1.70 × 10⁻⁴**: the fit's covariance describes the real scatter. That works here because the scatter in c_m (balance moment noise, q noise and the N noise transferred through the lever arm) is close to independent from point to point, as the fit assumes. If the noise were correlated (a drifting q, say), the fit would claim more precision than it has; repeating the experiment, as here, is the honest test.

(c) The truth is 0.248; 0.24805 − 0.248 = 4.5 × 10⁻⁵ = **0.30σ**. Yes, well within 2σ.

**Check.** a₀ and α_L0 agree with the generator's 6.0 /rad and −2.1° within 1σ.

**Common mistakes.** Fitting c_m against α instead of c_n; reporting 0.25 − slope with the wrong sign; quoting the fit uncertainty of a single run as the reproducibility without the seed study.

## 5. Where does the drag uncertainty come from?

Seed 42, first run at each angle, the instrument uncertainties of `NOISE` (N 0.20 N, A 0.05 N, M 0.010 N·m, α 0.05°, q 0.5 %).

(a)

| nominal α | c_l | c_d | u(c_d) | N | A | M | α | q |
|---|---|---|---|---|---|---|---|---|
| 0° | 0.2221 | 0.00638 | 0.00043 | 0.0 % | 78.8 % | 0.0 % | 20.7 % | 0.6 % |
| 4° | 0.6354 | 0.00688 | 0.00068 | 2.4 % | 30.4 % | 0.0 % | 66.9 % | 0.3 % |
| 8° | 1.0584 | 0.00952 | 0.00102 | 4.2 % | 13.5 % | 0.0 % | 82.1 % | 0.2 % |

At 4°, Monte Carlo (20 000 samples, seed 1) gives 6.73 × 10⁻⁴ against the linear 6.78 × 10⁻⁴, 0.8 % apart: the reduction is close enough to linear for first-order propagation.

(b) ∂c_d/∂α = c_l (Problem 2d), and c_l grows with α, so the angle's contribution grows with lift. q scales c_d by its relative error, 0.5 % of a small number; it matters for c_l (82 % of u(c_l) at 4°) but hardly for c_d.

(c) At 8° the sensitivity is 0.01847 per degree and the A contribution is 3.75 × 10⁻⁴, so **u_α = 0.020°** equalises them, giving **u(c_d) = 5.7 × 10⁻⁴**. Even with a perfect α, the A channel (and N) leave **3.9 × 10⁻⁴ at 4° and 4.3 × 10⁻⁴ at 8°**: no improvement in α alone can reach 0.0002. That needs a better axial channel too, or a wake survey, which measures drag directly and is how Ladson's NACA 0012 drag (used in `P10.06`) was obtained.

**Common mistakes.** Using the noise of q in pascals as a relative value; omitting the cross-coupling of α into drag; reading a Monte Carlo standard deviation from too few samples.

## 6. What transfers from the tunnel to the aeroplane?

| quantity | expectation | reason |
|---|---|---|
| a₀ | small correction | set by the potential flow; the boundary layer's displacement thickness, thicker at low Re, lowers it a few per cent |
| α_L0 | small correction | set by camber; viscous decambering at the trailing edge shifts it slightly |
| x_ac | transfers | near the quarter chord for any thin section; viscous effects move it by a few thousandths of chord |
| c_m,ac | small correction | camber-dominated; the boundary layer reduces the effective camber a little |
| c_d,min | does not transfer | skin friction falls with Re and transition moves forward: a factor of six in Re changes c_d by tens of per cent |
| c_l,max | does not transfer | stall is a boundary-layer separation; it depends strongly on Re and on surface condition |

The largest risk is c_d and c_l,max. The single most useful addition is a **Reynolds-number sweep with and without a transition trip**, so that the trend can be extrapolated and the boundary-layer state is known (`P10.06` shows how strongly the trip changes c_d).

**Common mistakes.** Assuming a coefficient is a property of the shape alone; extrapolating c_l,max by 6× in Re without data.

## Check-yourself answers

1. q = ½ρV², in pascals (N/m²). It is the pressure rise when the stream is brought to rest isentropically at low speed, and the scale on which every aerodynamic force is measured.
2. D = N sin α + A cos α and L = N cos α − A sin α; at α = 0, L = N and D = A.
3. Because the moment about a point depends on where the point is: M_P = M_O + (x_P − x_O)N. Moving the reference aft by Δx adds Δx·N to the nose-up moment, so the slope of c_m against c_n changes by Δx/c.
4. The point about which the pitching moment does not change with incidence. For a thin section it is near the quarter chord; it is found as x_ac = x_ref − c dc_m,ref/dc_n.
5. ∂c_d/∂α = c_l: a 0.1° error at c_l = 0.64 is 0.0011 in c_d, 16 % of the drag (here exactly 0.00111 for the seed-42 point at 4°).
6. Re = 6.98 × 10⁵ for the tunnel model and 4.52 × 10⁶ for the trainer: the model's boundary layer is more laminar and relatively thicker, so its drag and maximum lift do not transfer directly.
7. Because the data are generated by a stated program from a stated seed, with invented "true" values and invented noise; they test the reduction, not the airfoil. Calling them measurements would be presenting a simulation as a measurement.
