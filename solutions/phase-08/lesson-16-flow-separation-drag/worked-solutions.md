# P08.16  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe scripts that call it, given to the figures shown. Every drag coefficient is a correlation unless it is labelled exact.

## 1. A pitot mast (s = 15: d = 17.5 mm)

**Assumptions.** ISA 3000 m (ρ = 0.9091 kg/m³, μ = 1.694 × 10⁻⁵ Pa·s); a 2-D cylinder (end effects ignored, which overestimates the drag of a short mast a little); smooth surface.

(a) Re_d = 0.9091 × 60 × 0.0175/1.694 × 10⁻⁵ = **56 360**; c_d = 1 + 10 Re^(−2/3) = **1.007** (White's correlation; published curves near 1.2 would add 20 %); drag ½ × 0.9091 × 60² × 1.007 × 0.0175 × 0.3 = **8.65 N**; power DV = **519 W**.

(b) St = 0.212(1 − 12.7/56 360) = 0.212 (Roshko's correlation, used beyond its fitted range of 2000, where St stays at 0.20–0.21); f = 0.212 × 60/0.0175 = **727 Hz**.

(c) f = 90 Hz at **V = 7.4 m/s** (Re = 6990, St = 0.212). Every take-off and landing roll passes through it, but at 7.4 m/s the dynamic pressure is 25 Pa, 65 times smaller than at cruise, so the alternating force is small. Check that the structural damping is adequate, or raise the natural frequency above the shedding frequency at the highest speed at which lock-in could matter; a short fairing on the mast removes both the drag and the regular shedding.

**Check.** The drag scales as d: 8.65 N × (25/17.5) × (1/0.3) = 41.2 N per metre for the lesson's 25 mm strut, close to the lesson's 41.1 N (c_d is almost the same).

**Common mistakes.** Using the frontal area of a sphere (πd²/4) for a cylinder; forgetting the altitude density.

## 2. Derivations

(a) A point on the cylinder at angle θ from the front stagnation point is (−R cos θ, R sin θ), with outward normal n = (−cos θ, sin θ). The force of the fluid on the body is −∮(p − p∞) n R dθ per unit span; its x component is ∮(p − p∞) cos θ R dθ. Dividing by q d = q 2R: **c_d,p = ½∫₀^{2π} c_p cos θ dθ**. If c_p(π − θ) = c_p(θ), substitute θ' = π − θ: cos θ' = −cos θ, so the integral equals minus itself and is zero. The side force uses sin θ, which is odd under θ → −θ, so a c_p even in θ gives zero side force.

(b) Front, |θ| < θ_s, c_p = 1 − 4 sin²θ: ½∫_{−θs}^{θs}(cos θ − 4 sin²θ cos θ) dθ = ½[2 sin θ_s − (8/3) sin³θ_s] = sin θ_s − (4/3) sin³θ_s. Wake: ½ c_pb ∫_{θs}^{2π−θs} cos θ dθ = ½ c_pb[sin(2π − θ_s) − sin θ_s] = −c_pb sin θ_s. Sum: **c_d = sin θ_s − (4/3) sin³θ_s − c_pb sin θ_s**. Check: θ_s = π gives 0 (no separation); θ_s = 82°, c_pb = −1.2 gives 0.884, as the code does.

(c) c_d = 24/Re gives D = ½ρ_f V²(24μ/(ρ_f V d))(πd²/4) = 3πμVd. Equating with (ρ_p − ρ_f) g πd³/6: **V = (ρ_p − ρ_f) g d²/(18μ)**. In general, (ρ_p − ρ_f) g πd³/6 = ½ρ_f V² c_d πd²/4 ⇒ c_d V² = (4/3)(ρ_p − ρ_f) g d/ρ_f. Multiply by ρ_f²d²/μ²: **c_d Re² = (4/3)(ρ_p − ρ_f)ρ_f g d³/μ² = (4/3)Ar**. Since c_d is a function of Re only, Re at terminal velocity is a function of Ar only.

**Common mistakes.** A sign error in the normal (giving negative drag); forgetting buoyancy.

## 3. The ping-pong ball

(a) Morrison: **V_t = 9.34 m/s**, Re = 24 860, c_d = 0.395. Schiller–Naumann: **8.85 m/s**, Re = 23 550, c_d = 0.44. They differ by **5.3 %** in speed (11 % in c_d): the honest uncertainty of a correlation-based terminal velocity.

(b) 95 % of V_t after **1.77 s** and **10.6 m** of fall.

(c) At z = 2 m: **5.59 m/s**, **0.60** of V_t, at t = 0.667 s against 0.639 s in a vacuum: **29 ms later**. The drag at 2 m is (0.60)² = 36 % of the weight, so drag is clearly present, but the ball is far from terminal velocity. The experiment measures c_d through a 29 ms delay in landing time; with 1/240 s video timing (4 ms) that is about a 15 % measurement of the drag effect, and the result depends on the assumed c_d over a range of speeds. A drop of 10 m or more, or a lighter object, would measure the terminal velocity directly.

**Check.** The time integration reaches 9.344 m/s, the root of the force balance, after 14 s (the lesson).

## 4. A sounding-rocket payload coming down

(a) mg = ½ρV²c_d(πD²/4) ⇒ D = √(8mg/(πρc_dV²)) = **0.90 m** at sea level. At 3000 m (ρ = 0.9091) the same canopy descends at 6 × √(1.225/0.9091) = **6.96 m/s**.

(b) Smooth 0.12 m sphere, 2 kg (ρ_p = 2210 kg/m³): Morrison **142.6 m/s** (Re = 1.17 × 10⁶, c_d = 0.139); Schiller–Naumann **80.2 m/s** (Re = 6.6 × 10⁵, c_d = 0.44). The payload falls right through the drag crisis: Morrison puts the terminal speed on the low-drag side, Schiller–Naumann's constant 0.44 knows nothing of the crisis. Neither is inside its range: Schiller–Naumann's formula is fitted to Re ≤ 1000 and its 0.44 is a subcritical value; Morrison's fit stops at Re = 10⁶. The force balance has a single root with Morrison's curve, at 142.6 m/s. A real payload is not a smooth sphere: roughness, fins and tumbling all move the crisis.

(c) The higher figure, about **140 m/s** (and a note that it is uncertain by tens of per cent), because the impact energy (½ × 2 × 143² = 20 kJ) and the hazard scale with it, and the low-drag branch is physically possible. The 80 m/s estimate would understate the hazard by a factor of three in energy.

## 5. How thick should a fairing be?

(a) 25 mm strut, 60 m/s, 3000 m:

| t/c | chord (mm) | Re_c | c_d (on thickness) | drag (N/m) | saving |
|---|---|---|---|---|---|
| 0.10 | 250 | 805 000 | 0.118 | 4.81 | 8.5× |
| 0.15 | 167 | 537 000 | 0.094 | 3.84 | 10.7× |
| 0.25 | 100 | 322 000 | 0.081 | 3.32 | 12.4× |
| 0.33 | 76 | 244 000 | 0.089 | 3.64 | 11.3× |
| 0.50 | 50 | 161 000 | 0.155 | 6.33 | 6.5× |

(b) On the thickness, c_d,t = 2C_f(c/t)(1 + 2t/c + 60(t/c)⁴). The friction part grows like c/t as the fairing gets thinner (more wetted area per unit thickness); the form part 60(t/c)⁴ grows quickly as it gets thicker. The minimum is near t/c = 0.25.

(c) At t/c = 0.167: c_d = 0.0893, 3.65 N/m, **0.33 N per metre more** than the optimum (10 % more drag), plus a longer, heavier fairing.

## 6. What the wake model can and cannot tell you

(a) (1) The real c_p ahead of separation never reaches the potential-flow −3 (the subcritical minimum is nearer −1.2), because the displacement of the boundary layer and the wake change the outer flow; the model's front contribution is therefore wrong. (2) The wake pressure is neither constant nor steady: the shed vortices make it oscillate and the time-averaged base pressure varies round the rear. (3) Friction drag is left out (a few per cent at Re = 10⁵); and (4) White's fit itself is about 20 % below other published values near 10⁵.

(b) The base pressure, c_pb: it carries most of the drag (−c_pb sin θ_s is about 1.19 of the 0.88). One pressure tap at θ = 180° connected to a manometer referenced to the free-stream static pressure measures it; a ring of taps then gives the whole c_p(θ) and θ_s.

(c) At 60 m/s a golf ball (42.7 mm) has Re = 1.7 × 10⁵, below a smooth sphere's crisis (Morrison c_d = 0.44). The dimples trip the boundary layer to turbulence, which separates later and leaves a narrower wake: the golf ball has the lower drag, roughly half. For the strut, tripping would help only if Re were near the crisis (2–5 × 10⁵); at 8 × 10⁴ the gain is uncertain and small compared with the twelvefold saving of a fairing. Fair it.

## Check-yourself answers

1. Pressure drag is the stream-wise resultant of the pressure on the surface; friction drag that of the wall shear. A flat plate edge-on to the flow is friction-dominated; a plate across the flow (or a parachute) is pressure-dominated.
2. Its c_p is symmetric front to back, so the pressure on the rear pushes forward exactly as hard as that on the front pushes back. A real boundary layer cannot climb the rear adverse pressure gradient, separates, and leaves a low-pressure wake.
3. The wake term −c_pb sin θ_s: the low base pressure.
4. D = 3πμVd (c_d = 24/Re); V = (ρ_p − ρ_f) g d²/(18μ). Valid for Re ≪ 1 (within a few per cent below Re ≈ 0.1–0.5).
5. A sudden fall of c_d (from about 0.4–0.5 to 0.1 for a smooth sphere) near Re = 3–4 × 10⁵, when the boundary layer becomes turbulent before separating; a turbulent layer has more momentum near the wall, separates further back, and leaves a narrower wake with higher base pressure.
6. f = St V/d = 0.2 × 50/0.006 ≈ 1.7 kHz.
7. A correlation is a fit to someone's measurements of other objects in other conditions; it has scatter and a range. State its name (source) and the range of Re (and conditions) it was fitted to, and whether you are inside it.
8. Thinner fairings have more chord and so more friction per unit thickness; thicker ones have adverse pressure gradients steep enough to raise the form drag (and eventually separate). The two effects balance near t/c = 0.25.
