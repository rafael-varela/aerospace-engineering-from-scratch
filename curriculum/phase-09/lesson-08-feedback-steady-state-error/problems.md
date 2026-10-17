# P09.08  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical. Simulations are of a linear model, not of an aircraft: say so where you report them.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: your auto-throttle gains)

For the lesson's cruise plant, take K_p = 200 + 20s N per m/s and, for part (c), K_i = 20 + 2s N/m.

(a) Proportional control: give L(0), the steady-state error after a 2 m/s speed command, and the speed lost after a 200 N drag increase.
(b) Is the loop stable? Compare K_p with the critical gain.
(c) PI control: give the type, K_v, the error to a 0.5 m/s² speed ramp and the closed-loop poles.

Answer form: three numbers with units; yes/no with the ratio; a type, a constant, an error and four poles.

## 2. (Derivation: errors from the loop structure)

(a) Derive E = SR − GSD for unity feedback with the disturbance at the plant input, and show S + T = 1.
(b) For L = K₀N₀(s)/(sⁿD₀(s)) with N₀(0) = D₀(0) = 1, derive the table of steady-state errors for steps, ramps and parabolas for n = 0, 1, 2 from the final-value theorem.
(c) Derive the steady-state speed loss F/(c + C(0)) after a drag increase F for a controller without an integrator, and show it is zero with one.
(d) Give a loop for which lim sE(s) exists and is finite but e(t) does not converge, and say which condition of the theorem fails.

Answer form: four derivations, one counterexample.

## 3. (Unfamiliar scenario: a satellite's attitude loop is type 2)

A satellite axis has J = 10 kg m² and is controlled by a reaction wheel with a PD law τ = −(K_p θ + K_d θ̇), K_p = 0.1 N m/rad, K_d = 1 N m s/rad.

(a) Write L(s), give its type and the error constants. Where are the closed-loop poles?
(b) The attitude command follows a constant angular acceleration of 10⁻⁴ rad/s² (a slew profile). What steady-state error remains?
(c) Solar radiation pressure applies a constant torque of 10⁻⁵ N m. What steady attitude error does it cause? Why does the type-2 loop not reject it?
(d) Add an integral term K_i/s with K_i = 0.002 N m/(rad s). Is the loop stable, what is its type, and what is the disturbance error? Confirm with a 3000 s simulation.

Answer form: a transfer function, a type, constants and poles; an angle; an angle and a sentence; a verdict and two numbers.

## 4. (Numerical, uses your code: gain sweeps and the theorem's limits)

(a) For K_p = 100, 400, 1000, 2000, 4000 N per m/s, tabulate the step error and the least damping ratio of the closed loop.
(b) With K_p = 400, tabulate for K_i = 10, 40, 100, 200 N/m the ramp error for 0.5 m/s², the least damping ratio and the slowest real pole.
(c) Run `fvt_vs_simulation` for the PI ramp with dt = 0.1, 0.05, 0.025 s and explain the difference from the theorem.
(d) Find K_p for which the theorem would give a step error of 0.65 % and show that the loop is then unstable.

Answer form: two tables, three pairs and an explanation, one gain and one sentence.

## 5. (A climb is a disturbance)

The aircraft enters a 2° climb with the thrust unchanged. Treat the weight component mg sin γ as a drag increase.

(a) What force is that, and how much speed would the aircraft lose with no auto-throttle? How long would it take?
(b) Simulate the K_p = 400 and the PI loops. Report the final speed loss, the peak loss and its time, and the time to come within 0.1 m/s of the final value.
(c) What thrust must the engine supply in steady state, and why could a real auto-throttle fail to achieve the PI result?

Answer form: a force, a speed and a time; a table; a force and two sentences.

## 6. (Judgment: what the steady-state answer leaves out)

(a) The PI loop's sensitivity peaks at 1.405 at 0.627 rad/s. Give \|S\| at 0.01, 0.1, 0.3, 0.627, 1 and 3 rad/s and explain what it means for a gusty atmosphere.
(b) The 0.5 m/s² ramp to +10 m/s needs how much extra thrust? What happens to the PI loop if the engine cannot supply it?
(c) The model ignores the phugoid. Which of the lesson's numbers would you trust in a flight test, and which not?

Answer form: six numbers and two sentences; a force and two sentences; a short argument.
