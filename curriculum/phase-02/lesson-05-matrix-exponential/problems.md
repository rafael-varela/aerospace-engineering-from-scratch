# P02.05  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: exponentials by hand)

(a) For N = [[0, 1], [0, 0]] compute e^{Nt} from the series. (b) For A = [[0, ω], [−ω, 0]] with ω = s/10 rad/s, sum the series in closed form (group even and odd powers) and show that e^{At} is a rotation matrix; give e^{A·1 s} to four figures. (c) For D = diag(−1, −3) write e^{Dt}. (d) Check all three with `expm_taylor` and `expm_scaling_squaring`, and say for which of them `expm_eig` works and why.

Answer form: three matrices as functions of t, one numerical 2×2 matrix, one sentence on `expm_eig`.

## 2. (Derivation: the properties of Φ)

(a) Show from the series that dΦ/dt = AΦ = ΦA and Φ(0) = I. (b) Prove Φ(t)Φ(s) = Φ(t + s) by the uniqueness theorem (consider both sides as functions of t), and deduce Φ⁻¹(t) = Φ(−t). (c) Prove det Φ(t) = e^{tr(A)t} (hint: differentiate det Φ using the fact that d(det X)/dt = det X · tr(X⁻¹ dX/dt), or use the eigen form when A is diagonalisable and argue by continuity). (d) Show that if AB = BA then e^{A+B} = e^A e^B by multiplying the two series and using the binomial theorem, and compute both sides for A = [[0, 1], [0, 0]], B = [[0, 0], [1, 0]] to show the converse fails. (e) Derive the augmented-matrix identity (4) by computing the powers of M = [[A, B], [0, 0]].

Answer form: five derivations.

## 3. (Unfamiliar scenario: nutation of a spinning spacecraft)

A spin-stabilised spacecraft's transverse angular rates (ω₁, ω₂) obey ω̇₁ = −dω₁ + λω₂, ω̇₂ = −λω₁ − dω₂, where λ = (I₃ − I₁)Ω/I₁ is the nutation frequency set by the spin rate Ω and the inertia ratio, and d is the damping rate of a nutation damper. Take λ = 0.5 + s/50 rad/s and d = 0.02 s⁻¹.

(a) Write A and show, using Problem 2(d), that e^{At} = e^{−dt} R(λt) where R is a rotation matrix; state which two commuting matrices you split A into. (b) Give the nutation period, the time for the nutation amplitude to halve, and the number of nutation cycles in that time. (c) Evaluate e^{A·10 s} by hand and check with `state_transition`; verify det = e^{−2d·10}. (d) A transverse rate of 0.1 rad/s is imparted by a thruster misfire. How long until it is below 0.001 rad/s, and does the answer depend on λ? Is the system asymptotically stable according to `is_stable`?

Answer form: one matrix identity, three numbers, one 2×2 matrix with a determinant, one time with two yes/no answers.

## 4. (Numerical, uses your code: your rendezvous)

A chaser is at the same altitude as the station, (x₀, y₀) = (0, −(500 + 20s)) m, at rest relative to it, in the 6778 km orbit of `iss_orbit`.

(a) Where is the chaser after one orbit if nothing is done, and why? (b) Plan a two-impulse rendezvous with T = 0.4 orbit: give Δ**v**₁, Δ**v**₂ and the total, and the transfer duration in minutes. (c) Tabulate x and y every T/8 along the transfer with `propagate` and report the maximum radial excursion; explain its sign using the Clohessy–Wiltshire physics (which way does a chaser have to go to catch up?). (d) Compare the total Δv for T = 0.25, 0.5, 0.75 orbit. (e) Show that `rendezvous_delta_v` raises for T = 1 orbit and explain why from det Φ_rv.

Answer form: one position with a reason, two vectors and two scalars, a nine-row table with one explanation, three numbers, one explanation.

## 5. (Judgment: choosing the method and trusting the model)

(a) You need Φ(t) for a 50-state structural model at 10 000 time steps. Which of the three constructions do you use and why? What if the model is the 4-state Clohessy–Wiltshire matrix at a single T? (b) A colleague computes e^{At} for the short-period matrix at t = 60 s by `expm_taylor` with 30 terms and reports entries of order 10²⁰. Diagnose the problem in two sentences and say what the true entries are of order. (c) The two-impulse plan arrives at the origin to 10⁻¹³ m in the model. List, with magnitudes, three reasons the real chaser will miss by metres, in order of importance for a 1 km, half-orbit transfer at 400 km. (d) A controller updates its thrust command every 10 s. Explain why the zero-order-hold discretisation is exactly right for simulating it and why RK4 with the same 10 s step is not, and give one situation in which RK4 would still be the right choice.

Answer form: four short paragraphs with specific numbers or mechanisms.

## 6. (Hard: exact discretisation versus Euler, and the Taylor bound)

A damped oscillator ẍ + cẋ + kx = u has k = s² s⁻² and c = 0.2s s⁻¹ (so ζ = 0.1 for every s).

(a) Give the eigenvalues of A and show that the eigenvalues of the exact A_d = e^{Ah} are e^{λh}, with |e^{λh}| = e^{−ζω_n h} < 1 for every h > 0. (b) Forward Euler uses A_E = I + hA. Show its eigenvalues are 1 + λh and derive the largest h for which |1 + λh| < 1 (the stability limit), as a formula and a number for your s. (c) Tabulate |eigenvalue|_max of A_d and of A_E for h = 0.01, 0.05, 0.1 and 0.2 s, using `discretize` for A_d. Interpret. (d) The Taylor remainder after N terms is bounded by ‖A‖^N/N! · e^{‖A‖}. For ‖A·1 s‖₁ and N = 10, 20, 30, compute the bound and the actual error of `expm_taylor` against `scipy.linalg.expm`, and explain why the actual error is much smaller than the bound and yet still unacceptable. How many squarings does `expm_scaling_squaring` use for this matrix?

Answer form: two eigenvalues and a proof, a formula and a number, an eight-entry table with interpretation, a six-entry table with two explanations and a count.
