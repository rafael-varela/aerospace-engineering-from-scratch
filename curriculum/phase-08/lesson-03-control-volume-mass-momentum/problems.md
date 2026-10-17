# P08.03  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical. For every control-volume problem, draw the control volume, mark every face the fluid crosses with its outward normal, and list the forces before you write a balance.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

Engine and vehicle numbers in this set are illustrative values chosen to be realistic, not the data of a particular engine.

## 1. (Warm-up: a rocket engine at sea level and in vacuum)

A small liquid rocket engine has a propellant mass flow of 12 kg/s, a nozzle exit velocity of 2800 m/s, an exit pressure of 40 kPa and an exit area of A_e = 0.05 + s/1000 m².

(a) Draw the control volume and derive the thrust from the momentum balance.
(b) Compute the thrust at sea level (101 325 Pa), at 11 km (22 632 Pa) and in vacuum, and the effective exhaust velocity c = F/ṁ and specific impulse I_sp = c/g₀ for sea level and vacuum.
(c) At what ISA altitude is the pressure thrust zero? What is the nozzle called below that altitude, and what can happen to the flow inside it?

Answer form: a derivation; three thrusts, two velocities and two impulses; an altitude, a name and a sentence.

## 2. (Derivation: the box around a jet engine)

(a) Starting from the Reynolds transport theorem for a system property B with intensive value b, dB_sys/dt = d/dt ∫_CV ρb dV + ∮_CS ρb (u_rel·n) dA, write the mass and x-momentum balances for a steady flow through a fixed control volume.
(b) Draw a cylindrical box of cross-section A_box around a turbojet in the engine's frame. Write the mass flow through each face (inlet, jet, outer exit, side, fuel line). Show that the side face must carry ṁ_side = ρ₀V₀A_e − ṁ_a, and interpret its sign.
(c) Write the x-momentum balance and show that the force through the pylon is F = ṁ_eV_e − ṁ_aV₀ + (p_e − p_a)A_e, independent of A_box. Where did the side face's momentum go?
(d) Repeat (c) in a frame moving with velocity −W along x, so every x-velocity increases by W while the normal velocities relative to the faces do not change. Show that the force is unchanged, and that this needs Σρu_nA = 0.

Answer form: two equations; five mass flows and a sentence; a derivation and a sentence; a derivation.

## 3. (Unfamiliar scenario: a quadcopter in hover)

A quadcopter of mass m = 1.5 + s/10 kg hovers with four rotors of diameter 0.25 m. Model each rotor as an actuator disc: a thin disc that adds momentum to the air passing through it, with the air far above at rest and the wake far below moving at w.

(a) Apply mass and momentum balances to a control volume around one rotor's stream tube to show T = ṁ w. Apply the energy balance (kinetic energy flux leaving = power supplied) to show P = ½ṁw². With P = T v_i at the disc, show v_i = w/2.
(b) Compute v_i, w and the ideal power per rotor and in total at sea level (ρ = 1.225 kg/m³).
(c) Repeat at 3000 m (ρ = 0.9091 kg/m³). By what fraction does the hover power rise? Why does thin air cost power?
(d) The measured electrical power in hover is about twice the ideal total. Name two losses the actuator disc ignores.

Answer form: three derivations; three numbers per rotor and a total; two numbers and a sentence; two losses.

## 4. (Numerical, uses your code: your own turbojet in cruise)

Use the lesson's cruise engine at 11 km and 236 m/s (`CRUISE` in the code) with the air mass flow changed to ṁ_a = 15 + s/5 kg/s.

(a) Compute A_e and the gross, ram-drag, pressure and net thrusts with `engine_case`.
(b) Build `turbojet_control_volume` boxes of 0.5, 5 and 50 m² and report the mass residual and the difference between the box force and the formula.
(c) The nozzle is underexpanded (p_e > p_a). Suppose a longer, divergent nozzle expanded the exhaust isentropically to p_a, with γ = 1.33 and c_p = 1150 J/(kg·K) for the hot gas (typical values): V_e² gains 2c_pT_e[1 − (p_a/p_e)^((γ−1)/γ)]. Compute the new V_e, A_e and net thrust. By how much does thrust rise? Why is "same V_e, p_e = p_a" not the right comparison?
(d) In one sentence: why is the best nozzle the one with p_e = p_a, even though pressure thrust is then zero?

Answer form: an area and four forces; a three-row table; two numbers, a thrust, a percentage and a sentence; a sentence.

## 5. (Numerical, uses your code: non-uniform jets and gravity losses)

(a) The engine of Problem 4 has a non-uniform exhaust profile. Using `profile_fluxes`, compute the momentum-flux correction factor β of a 1/7th-power-law profile and of a laminar (parabolic) profile, and the error in gross thrust if you write ṁV̄ instead of βṁV̄ for each. Which is realistic for a turbojet exhaust, and why?
(b) The lesson's sounding rocket (m₀ = 100 kg, 60 kg of propellant, c = 2000 m/s, vertical, no drag) can burn its propellant in 4, 8, 16, 32 or 60 s. For each, compute the initial thrust-to-weight ratio, the burnout speed by `rocket_burn` and by Tsiolkovsky with gravity loss, and the burnout altitude.
(c) Why does a short burn give a higher burnout speed? What limits how short it can be in practice?

Answer form: two factors and two thrust errors with a sentence; a five-row table; two sentences.

## 6. (Judgment: what does the test stand measure?)

A turbojet is run on a ground test stand: the engine hangs from a frame on flexures and a load cell measures the horizontal force. The engineer reports the load-cell reading as "the thrust".

(a) Draw the control volume that turns the load-cell reading into a thrust. Which terms of the turbojet thrust formula are zero on the stand, and which are not?
(b) Name three effects that make the installed thrust in flight differ from the stand reading at the same engine setting.
(c) The intake draws air from the test cell, which then moves past the engine at 5 m/s. Does that change the reading? Which term does it bring in?
(d) The exhaust profile is non-uniform and the engineer computes thrust as ṁV̄ from a mass-flow meter and a mean exhaust speed instead of reading the load cell. Which error would you check first?

Answer form: a sketch and a list; three effects; a sentence and a term; a sentence.
