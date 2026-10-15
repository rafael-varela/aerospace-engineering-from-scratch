# P02.08  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
(a) F = 3/s + 2/(s + a) + 4/(s² + 16); substitute s = 2; `laplace_numeric(lambda t: 3 + 2*np.exp(-a*t) + np.sin(4*t), 2.0)`. (b) Complete the square; match (s + 2)/((s + 2)² + 3²) with the damped-cosine entry (a = 2, ω = 3); `residues([1, 2], [1, 4, 13])` gives −2 ± 3i with residues ½ each; `inverse_laplace([1, 2], [1, 4, 13], t)`. (c) The integrand of L{3} does not decay at s = 0; F has a pole at s = 0 from the constant term.

## 2.
(a) u = e^{−st}, dv = f′dt; boundary term f(t)e^{−st} → 0 needs |f| ≤ Me^{at} with Re s > a; apply once more to f′. (b) Use (a) term by term; collect Y. (c) Cover-up for the pole at 0; residue at p from N(p)/D′(p) with D′ = 3s² + 4ζω_n s + ω_n² or the cover-up form ω_n²/[p(p − p̄)]; combine the conjugate pair into 2|r|e^{−ζω_n t}cos(ω_d t + arg r); convert to the sin/cos form. (d) ẏ = e^{−ζω_n t}(ω_n/√(1 − ζ²)) sin ω_d t; zero at ω_d t = π; substitute into y.

## 3.
(a) P_true = R/s²; P_ind = R/(s²(τs + 1)); poles 0 (double), −1/τ; the simple-pole formula divides by D′(0) = 0. (b) R/(s²(τs + 1)) = R[1/s² − τ/s + τ²/(τs + 1)]; invert: R[t − τ + τe^{−t/τ}]; for t ≫ τ, p_ind ≈ R(t − τ). (c) Pressure lag Rτ divided by −12 Pa/m is the altitude lag; during a climb the true pressure is lower than the indicated one, so the indicated altitude is lower than the true one. (d) `solve_ivp(lambda t, p: [(R*t - p[0])/tau], (0, 10), [0.0], t_eval=[1, 5, 10])`. (e) E(s) = R/s² − R/(s²(τs + 1)) = Rτ/(s(τs + 1)); sE → Rτ as s → 0.

## 4.
(a) `poles(transfer_function(25, z)[1])`, `overshoot(z)`, `peak_time(25, z)`, `settling_time(t, step_response(25, z, t), 1.0)` on `np.linspace(0, 1, 100001)`. (b) Compare with 0.10 and 0.250 s. (c) `y - delayed_step_response(25, z, t, 0.05)`; `max` and `argmax`. (d) `frequency_response(num, den, np.array([3, 10, 25]))`; |H| and −angle in degrees; a ratio within a few per cent of 1 and a lag of a few degrees is "undistorted".

## 5.
(a) Valve flow saturation gives a rate limit; the normalised response changes with amplitude; superposition fails. (b) Overshoot and peak time become rate-limit numbers; measure the maximum rate (deg/s) and the time to reach 90 % at full amplitude. (c) Overshoot < 5 % needs ζ > 0.69, i.e. pole angle from the negative real axis below arccos 0.69; settling 150 ms needs ζω_n > 4/0.15; current poles have ζω_n = 15: increasing ζ at fixed ω_n moves along the circle of radius 25 toward the real axis; check whether the point with ζ ≈ 0.69 has real part beyond the line. (d) Stiffness k enters ω_n² ∝ k; ω_n falls by √0.7; poles move toward the origin along the same angle; settling time and bandwidth suffer, overshoot does not.

## 6.
(a) sin ω_p t + sin ω_p(t − t_p)u_{t_p} is zero for t > t_p since sin(ω_p t − π) = −sin ω_p t; transform the two terms with the shift rule. (b) U_num = [ω_p], U_den = [1, 0, ω_p²]; `response_with_initial_conditions(num, den, U_num, U_den, 0, 0)`; y_sin = `inverse_laplace` on t and on t − t_p (causal, so zero before t_p); add. (c) `solve_ivp(actuator_rhs, ..., args=(wn, zeta, pulse), max_step=5e-4)`; `convolve_causal(impulse_response(...), u_samples, dt)`. (d) Loop over t_p; compare short-pulse peaks with (2t_p/π) × max g and long-pulse peaks with 1 + a fraction of the overshoot.
