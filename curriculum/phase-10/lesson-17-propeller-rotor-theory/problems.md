# P10.17  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact (closed form) and which are numerical. Draw the velocity triangle before writing any blade-element equation.

Several problems use your personal seed **s**: the day of the month you were born (1–31). ISA sea level: ρ = 1.225 kg/m³.

**Safety.** Every problem here is a calculation. Do not run propellers or rotors to check your answers.

## 1. (Warm-up: the bound for a bigger UAV)

A 2 kg UAV cruises at 15 m/s on a 0.30 m (12 in) propeller that must deliver T = 1.5 + s/20 N.

(a) Use `actuator_disk` to find the induced velocity at the disk, the slipstream speed far behind, the ideal power and the ideal efficiency.
(b) Express the same operating point at 6000 rpm as C_T and J and check (5) of the lesson against (a).
(c) What is the least power that could hold 10 N of static thrust with this disk? What is the induced velocity?

Answer form: four numbers; two coefficients and an efficiency; a power and a velocity.

## 2. (Derivation: from conservation laws to a single residual)

(a) Derive w = 2v for the actuator disk from mass, momentum and energy conservation, and then η_i = 2/(1 + √(1 + 8C_T/(πJ²))).
(b) Starting from the blade-element loads (7) and the annulus balances (8) of the lesson, derive k = σ′C_n/(4F sin²φ), u_a = Vk/(1 − k), and the corresponding expressions for k′ and u_t. Then derive the residual R(φ) = sin φ(1 − k) − (V/(Ωr)) cos φ (1 + k′).
(c) Show that with cd = 0, u_t = 0 and F = 1 each annulus absorbs dP = dT(V + u_a). Then prove that, at fixed total thrust, uniform u_a gives the least total power. (Hint: minimise ∫ dT (V + u_a) subject to fixed ∫ dT with dT ∝ r(V + u_a)u_a.)

Answer form: three derivations.

## 3. (Unfamiliar scenario: a quadcopter in hover)

Use four of the reference 10 × 4.5 propellers to hover the 1.2 kg airframe as a quadcopter (no wing, all weight on the rotors).

(a) Thrust per rotor? Find the rpm with `rpm_for_thrust(blade, T, 0.0)` and the shaft power per rotor with `bem`.
(b) Compute the figure of merit FM = P_ideal/P for one rotor. Which loss would you attack first to improve it?
(c) With 80 % motor-and-controller efficiency, what electrical power does the hover need? How long would an 11.1 V, 2.2 A·h battery last if 80 % of its energy is usable?
(d) Compare with the fixed-wing cruise power of the lesson (21.4 W shaft). Explain the ratio in one paragraph using momentum theory.

Answer form: a thrust, an rpm and a power; a number and a sentence; a power and a time; a ratio and a paragraph.

## 4. (Numerical, uses your code: choosing the pitch)

Replace the blade's pitch by p = 0.152 m (a 10 × 6 propeller), keeping the chord.

(a) Find the rpm, J, η and shaft power at the UAV's cruise (12 m/s, the lesson's thrust).
(b) Sweep J at that rpm and find the peak efficiency and where it is. Compare with the 10 × 4.5.
(c) Compute the static thrust of each propeller at its own cruise rpm. Which propeller would you fit, and what does the choice cost?
(d) Repeat (a) for your own pitch p = 0.100 + s/500 m.

Answer form: four numbers; a peak η and its J for each propeller; two thrusts and a recommendation; four numbers.

## 5. (Numerical, uses your code: when does the classical iteration fail?)

(a) Run `bem_fixed_point` at the cruise point for n_el = 40, 80, 160, 320, 640 and relaxation 0.5, 0.3, 0.2, 0.1, 0.05 (max_iter = 4000). Tabulate which cases converge. What is the largest safe relaxation at each n_el?
(b) Print F for the last three elements at n_el = 160. Explain, using (8), why the update of u_a becomes very sensitive in those elements.
(c) Check that `bem` (bisection on R(φ)) converges in every case, and compute the observed order of the thrust error for n_el = 10–160 against n_el = 1280.

Answer form: a table and a rule; three numbers and an argument; a statement and four orders.

## 6. (Judgment: what does the model not know?)

(a) The section model's drag level is uncertain by a factor of 1.5. Run the cruise point with cd0 = 0.010, 0.020 and 0.030 and report rpm, η and power. Which output is robust and which is not?
(b) The sections work at Re ≈ 6 × 10⁴. Name two phenomena from `P10.14` that the section model does not represent, and say in which direction each would move η.
(c) At the cruise point, at what angle of attack and cl does the 0.75R section work? Is the blade operating near its best lift-to-drag ratio? What would a designer change?
(d) List two measurements (made in a facility, not at home) that would test the BEM prediction, and say which model assumption each one would test.

Answer form: a table and a sentence; two phenomena with directions; an angle, a cl and a design change; two measurements.
