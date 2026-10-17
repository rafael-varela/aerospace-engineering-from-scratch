# P08.04  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical. Before you use Bernoulli's equation, name the streamline and check each of its four conditions (steady, inviscid, incompressible, no shaft work between the two points).

Several problems use your personal seed **s**: the day of the month you were born (1–31).

Shared numbers: water at 20 °C ρ = 998.2 kg/m³; air at 20 °C, 101 325 Pa, ρ = 1.204 kg/m³. LOX density 1141 kg/m³ (typical). Knots appear because pilots use them: 1 kt = 1852/3600 m/s.

## 1. (Warm-up: a pitot tube on a water manometer)

A pitot-static probe in a small wind tunnel (air at 20 °C) is connected to a vertical U-tube water manometer that reads Δh = 30 + s mm.

(a) Find the pressure difference and the airspeed.
(b) If you can read the manometer to ±1 mm, what is the uncertainty in the airspeed? Show that the relative uncertainty in V is half the relative uncertainty in Δh.
(c) Which of Bernoulli's four conditions is least well satisfied along the stagnation streamline in front of the probe, and does it matter here?

Answer form: a pressure and a speed; an uncertainty with a derivation; a sentence.

## 2. (Derivation: Bernoulli, the pitot tube and the energy equation)

(a) From the momentum balance on a small fluid element moving along a streamline, derive ∂p/∂s + ρV ∂V/∂s + ρg ∂z/∂s = 0 for steady inviscid flow, and integrate it for constant ρ. Name each assumption where you use it.
(b) Show that a pitot-static tube measures p₀ − p = ½ρV² in incompressible flow, and that the airspeed indicator's calibration law (the lesson's q_c relation) expands as q_c = ½ρV²(1 + M²/4 + …) at low Mach number.
(c) Define EAS and show that an indicator calibrated with sea-level density reads EAS when the flow is incompressible and the ports are perfect. Why do pilots fly by indicated airspeed for take-off and landing rather than by true airspeed?
(d) Starting from the steady-flow energy equation for a pipe, show that it reduces to Bernoulli when there is no shaft work, no heat transfer and no loss, and identify the head loss h_L as the drop in total head.

Answer form: a derivation; a derivation and an expansion; a definition and a sentence; a derivation.

## 3. (Unfamiliar scenario: a rocket turbopump)

A turbopump raises liquid oxygen from the tank outlet (3 bar, 3 m/s) to the injector manifold (50 bar, 20 m/s, 0.5 m higher), at a mass flow of 20 + s/2 kg/s.

(a) Find the volume flow and the pump head with `shaft_head`. Which term dominates?
(b) Find the shaft power for a pump efficiency of 70 %. Compare with the simple estimate ṁΔp/(ρη).
(c) Where does the 30 % that is not head go? What happens to the liquid's temperature, and why does that matter for a cryogenic propellant?
(d) The turbine that drives the pump is fed by the same propellants. Why do engine designers care about every kilowatt of pump power?

Answer form: a volume flow, a head and a sentence; two powers; two sentences; a sentence.

## 4. (Numerical, uses your code: your own pitot-static system)

Your aeroplane flies at pressure altitude h = 2000 + 100s m (ISA) at a true airspeed of 50 + s m/s.

(a) Compute the Mach number, q, q_c, EAS, CAS and IAS (perfect ports) with `airspeeds_from_tas`, in m/s and knots.
(b) The static port sits where c_p = −0.03. Compute the IAS, the IAS − CAS error in knots, and the altimeter error.
(c) The pilot reads 100 kt IAS with that port. Find the true airspeed (with a root finder on `airspeeds_from_tas(..., cp_static=-0.03)`), and compare it with what you get if you ignore the port error and treat 100 kt as CAS (`tas_from_cas`). Compare also with the rule of thumb TAS = IAS/√σ.
(d) Which of these errors would matter for navigation (time over a fix)? For the stall speed?

Answer form: six numbers; three numbers; three speeds; two sentences.

## 5. (Numerical, uses your code: where Bernoulli's airspeed fails, and a contraction)

(a) With `airspeeds_from_tas` and `bernoulli_airspeed`, find the true airspeed at which applying incompressible Bernoulli (with the true density) to the impact pressure q_c gives a 1 % error at sea level and at 3000 m, and the speed for a 5 % error. Show that the error depends on Mach number alone and give the two critical Mach numbers.
(b) A wind-tunnel contraction runs from A_in = 0.5(6 + s/5) m² to a 0.5 m² test section in 2 m, with 30 m/s in the test section. Compute the pressure along it with `euler_pressure_along_duct` for n = 25, 50, 100, 200 and the error against `bernoulli_pressure_along_duct`, with observed orders.
(c) The tunnel's speed is set from the pressure drop across the contraction (a "piezometer ring"). Show that the drop is q_test(1 − 1/CR²) with CR the contraction ratio, and invert it with `venturi_flow_rate`. Why do tunnels use this rather than a pitot in the test section?

Answer form: four speeds and two Mach numbers with a sentence; a four-row table with orders; a derivation, a flow rate and a sentence.

## 6. (Judgment: where does Bernoulli fail?)

Say whether Bernoulli's equation can be applied between the two points, and if not which condition fails and what to use instead: (a) far upstream of the aeroplane and the stagnation point on the pitot tube; (b) across the propeller disc; (c) from the free stream to a point inside the boundary layer on the wing; (d) from the free stream to the separated wake behind a fixed landing-gear leg; (e) along the axis of a pitot tube on an aircraft at M = 0.8; (f) along a streamline in a gust that changes the flow in 0.1 s around a 1.5 m wing chord at 60 m/s. For the static port, explain why it is placed on the fuselage side where c_p ≈ 0 and not on the nose.

Answer form: six labelled lines and a paragraph.
