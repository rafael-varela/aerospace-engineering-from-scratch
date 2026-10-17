# Phase 02 Exam: Mathematics II: Differential Equations, Transforms & Probability

Integrated assessment of the whole phase. Suggested self-assessment mode: closed reference, a time limit chosen and recorded in advance (the exam has not been piloted, so no duration is suggested), no code except in Part C, no computer algebra in Part B; a scientific calculator (or Python used only as a calculator) is allowed for the arithmetic in Part B, since the probability questions need values of the normal and t distributions, which you may take from a table. There is no enforced lock; the point is to find out what you can do unaided. Solutions and pass thresholds are in `solutions/phase-02/exam/`.

Report your result per competency: (1) mathematical and physical reasoning, (2) analytical problem solving, (3) computational competence, (4) engineering interpretation. Questions are tagged with the competencies they assess. Every question is written over an aerospace situation; the mathematics is what is being tested, the framing is there so that you also practise deciding which mathematics applies. Angles are in radians unless a degree sign is written. SI units throughout. Every data set in this exam is invented for the exam and is labelled synthetic; nothing here is a measurement.

## Part A: Concepts

**1.** (1, 4) A colleague linearises the undamped pendulum θ̈ + (g/ℓ) sin θ = 0 about its two equilibria. At θ = π the Jacobian has eigenvalues ±√(g/ℓ), and they conclude "unstable, so the inverted position cannot be held". At θ = 0 the eigenvalues are ±i√(g/ℓ), and they conclude "a centre, so the hanging pendulum oscillates forever, even with a little bearing friction". In eight sentences or fewer: write the system as two first-order equations, give the Jacobian at each equilibrium, say which of the two conclusions linearisation justifies and which it does not, and why (what theorem is being used and what its hypothesis is). Then say what adding a damping term c θ̇ does to the classification at θ = 0 and at θ = π, using the trace and determinant.

**2.** (1) A wing accelerometer is sampled at f_s = 500 Hz for N = 1000 samples and the signal is written as a phasor, a(t) = Re(Â e^{iωt}) with Â complex. (a) State what |Â| and arg Â mean physically and why the derivative ȧ(t) has phasor iωÂ. (b) A 380 Hz vibration source is present. At what frequency does it appear in the discrete spectrum, and why? Give the general rule for where a frequency f > f_s/2 lands. (c) What is the frequency resolution Δf of the DFT of this record, and what two changes to the acquisition would each halve it? (d) Say in one sentence what anti-aliasing filter you would place before the sampler and why no digital processing after sampling can replace it.

**3.** (4) A supplier's batch of aluminium-lithium plate must have a mean yield strength of at least 450 MPa. A test engineer draws eight coupons, finds a sample mean of 441 MPa, runs a one-sided t-test and reports "p = 0.019, so there is a 1.9% probability that the batch is in spec". (a) State the null and alternative hypotheses this test uses, and say precisely what p = 0.019 is the probability *of*, and why the engineer's sentence is wrong. (b) Define the Type I and Type II errors for this acceptance test and say which one the supplier worries about and which one the airframe designer worries about. (c) The engineer proposes to "fix" a borderline batch by testing eight more coupons and keeping whichever eight look better. Say what this does to the Type I error rate and name the correct procedure.

**4.** (1) The Joukowski map z = ζ + a²/ζ turns the potential flow around a circle into the flow around an airfoil. (a) State the Cauchy–Riemann equations for f(ζ) = u + iv and show in two lines that u satisfies Laplace's equation ∇²u = 0. (b) Explain why that fact is what makes conformal mapping useful for inviscid incompressible flow: if W(z) is a complex potential satisfying Laplace's equation in the airfoil plane, why does W(z(ζ)) satisfy it in the circle plane? (c) State where the map fails to be conformal (where dz/dζ = 0), what happens to angles there, and what feature of the airfoil that produces.

## Part B: Analysis

Work every derivation by hand and show the steps. Numbers to three significant figures unless more are needed to show a check.

**5.** (2) First-order ODEs: the vacuum burn. A sounding rocket has lift-off mass m₀ = 250 kg, dry mass 100 kg, mass flow ṁ = 7.5 kg/s and exhaust velocity v_e = 2200 m/s; take g = 9.81 m/s² constant and neglect drag.
(a) Write and solve the separable equation for the mass m(t) and give the burn time t_b.
(b) From Newton's second law for the rocket, ṁ v_e − m g = m dv/dt, derive v(t) during the burn by direct integration (the integral ∫ dt/(m₀ − ṁt) does it) and show it is the rocket equation with a gravity-loss term. Evaluate the ideal Δv = v_e ln(m₀/m_dry), the gravity loss, and the burnout speed.
(c) After burnout the rocket coasts with linear drag, m dv/dt = −m g − c v, with m = 100 kg and c = 25 N·s/m. Solve this linear equation by integrating factor for a general v(0) = v_b, identify the time constant and the terminal velocity, and show that the solution approaches the terminal velocity for every initial condition (the equilibrium of the autonomous equation is stable). With v_b from (b), how long after burnout is the apex (v = 0)? Say in one sentence why the linear-drag model is the wrong one for this rocket at that speed.

**6.** (2) Second-order linear ODEs: landing gear. A landing-gear leg is modelled as m ẍ + c ẋ + k x = 0 with m = 600 kg, k = 2.4 × 10⁵ N/m and c = 6000 N·s/m.
(a) Write the characteristic equation and compute ω_n, the damping ratio ζ, the critical damping c_c and the damped frequency ω_d. Classify the response.
(b) Derive the solution for x(0) = 0.05 m, ẋ(0) = 0, in the form e^{−ζω_n t}(A cos ω_d t + B sin ω_d t), and give A and B.
(c) Give the decay time constant 1/(ζω_n), the damped period, the logarithmic decrement δ = ln(x_n/x_{n+1}), and show that the first minimum occurs at t = π/ω_d with x = −x₀ e^{−πζ/√(1−ζ²)}; evaluate it. Evaluate x(0.2 s) as a check.
(d) The designer wants the first rebound to be no more than 10% of the initial compression. Find the ζ required and the c that gives it.

**7.** (2) Forced response and resonance: a panel. A fuselage panel behaves as m ẍ + c ẋ + k x = F₀ cos ωt with m = 2 kg, k = 8 × 10⁴ N/m, ζ = 0.02 and F₀ = 10 N.
(a) Using the complex-gain (phasor) method, derive the steady-state amplitude X(ω) = F₀/|k − mω² + icω| and the phase lag. Give c and ω_n.
(b) Evaluate the static deflection F₀/k, the amplitude and phase at ω = ω_n, and show that at resonance X = F₀/(cω_n) exactly (the spring and inertia forces cancel and only the damper resists). Give the magnification factor 1/(2ζ).
(c) Show that the amplitude peak is at ω = ω_n √(1 − 2ζ²), not at ω_n, and evaluate how different the two frequencies and the two amplitudes are here. Give the half-power bandwidth 2ζω_n.
(d) Evaluate the amplitude and phase at ω = 0.9 ω_n, and say in one sentence what this tells you about how precisely a flutter-test frequency sweep must be spaced to find a peak with ζ = 0.02.

**8.** (1, 2) Systems and phase portraits. Two linear systems ẋ = A x model, respectively, a short-period pitch mode and a divergent coupled mode:

```
A₁ = [ −1   1 ]        A₂ = [ 0   1 ]
     [ −4  −2 ]             [ 2   1 ]
```

(a) For each matrix compute the trace, the determinant and the discriminant tr² − 4 det, classify the equilibrium at the origin (node, saddle, spiral, centre; stable or unstable), and sketch the phase portrait.
(b) For A₁ give the eigenvalues, the equivalent damping ratio and undamped natural frequency of the mode, and the period of the oscillation.
(c) For A₂ give both eigenvalues and eigenvectors, say along which direction trajectories leave and along which they approach, and explain why a saddle is unstable for every initial condition except a set of measure zero.

**9.** (2) Matrix exponential. The undamped oscillator ẋ = A x with A = [[0, 1], [−ω², 0]] (state x = (position, velocity)).
(a) Compute A², A³ and A⁴, and sum the series e^{At} = Σ (At)ⁿ/n! to show that Φ(t) = e^{At} = [[cos ωt, sin(ωt)/ω], [−ω sin ωt, cos ωt]].
(b) Verify Φ(0) = I, Φ̇(t) = A Φ(t) and det Φ(t) = 1, and say why det Φ = e^{t·tr A} in general and what det Φ = 1 means for areas in the phase plane.
(c) For ω = 2 rad/s, t = 0.5 s and x(0) = (0.1 m, 0), evaluate Φ(0.5) and x(0.5). State the composition property Φ(t₁ + t₂) = Φ(t₁)Φ(t₂) and say why it is the reason the state-transition matrix is used for propagating a state at fixed steps.

**10.** (2) Laplace transform: a servo actuator. An actuator obeys ÿ + 5ẏ + 6y = 6u(t), with u the unit step and y(0) = ẏ(0) = 0.
(a) Transform the equation, give the transfer function G(s) = Y(s)/U(s), its poles and its DC gain.
(b) Expand Y(s) = 6/(s(s + 2)(s + 3)) in partial fractions (show the cover-up or the simultaneous equations) and invert to y(t).
(c) Check y(0) and ẏ(0) from your y(t), and the final value both from y(t) and from the final-value theorem. Evaluate y(1 s) and y(2 s), and give the time to reach 95% of the final value to two figures.
(d) Say which pole dominates the settling and why, and what the response would look like if the 6 in the denominator were replaced by 25 (poles become complex): name the two quantities that then describe it.

**11.** (2) Fourier series: pulsed thrust. A pulsed thruster delivers f(t) = A for 0 < t < T/2 and 0 for T/2 < t < T, repeating with period T, with A = 400 N.
(a) Derive a₀, a_n and b_n for n ≥ 1 by integrating against cos(2πnt/T) and sin(2πnt/T), and state which harmonics vanish and why (write f as a constant plus a square wave).
(b) Evaluate the mean thrust and the amplitudes of the first, second, third and fifth harmonics.
(c) Compute the mean-square value of f directly from the time waveform, then by Parseval's identity for the partial sum through n = 5, and give the fraction of the mean-square captured. Repeat for the sum through n = 1.
(d) State what happens at t = T/2 to the partial sums as more terms are added (name the phenomenon) and what value the series converges to there.

**12.** (2) PDEs: a re-entry tile. The temperature excess u(x, t) in a slab of thickness L = 0.05 m with both faces held at the reference temperature obeys u_t = α u_xx with α = 5 × 10⁻⁷ m²/s (a ceramic insulator), u(0, t) = u(L, t) = 0 and u(x, 0) = u₀ [sin(πx/L) + 0.2 sin(3πx/L)].
(a) Substitute u = X(x)T(t), separate, and derive the two ODEs with the separation constant. Explain why the constant must be negative for non-trivial solutions that satisfy both boundary conditions, and find the eigenvalues and eigenfunctions.
(b) Write the general solution and use the orthogonality of the sines to pick out the coefficients for this initial condition.
(c) Give the time constant τ₁ of the fundamental mode and that of the third mode. How large is the third-mode contribution relative to the first at t = τ₁? How long until the centre-line temperature excess falls to 10% of u₀ (fundamental only, and say why that is justified)?
(d) State how τ₁ scales with L and with α, and what that implies for doubling the tile thickness.

**13.** (1, 2) Conformal maps: the Joukowski airfoil. Take the map z = ζ + 1/ζ (a = 1).
(a) Show that the unit circle ζ = e^{iθ} is mapped onto the segment −2 ≤ x ≤ 2 of the real axis (a flat plate) and that each interior point of the segment is covered twice.
(b) A circle of radius R passes through ζ = 1 and is centred at ζ_c = −0.1 + 0.1i. Give R, the angle β = arcsin(0.1/R), the trailing-edge point z(1), the leading-edge point (image of the other crossing of the real axis), and the chord.
(c) Show that dz/dζ = 0 at ζ = ±1 and deduce the interior angle of the airfoil at the trailing edge.
(d) The Kutta condition makes the circulation Γ = 4πUR sin(α + β). For U = 30 m/s and α = 4°, evaluate Γ and the lift per unit span L′ = ρUΓ with ρ = 1.225 kg/m³. Divide by ½ρU²·chord to get a lift coefficient and compare it with the thin-airfoil value 2π(α + β).

**14.** (2) Calculus of variations: the minimum-time descent. A bead slides without friction from the origin under gravity along a curve y(x) (y measured downward) to a point (x₁, y₁).
(a) Starting from J[y] = ∫ F(y, y′) dx and a variation y + εη with η = 0 at the ends, derive the Euler–Lagrange equation ∂F/∂y − d/dx (∂F/∂y′) = 0.
(b) Show that when F has no explicit x-dependence, F − y′ ∂F/∂y′ is constant along a solution (the Beltrami identity).
(c) Apply (a) to F = √(1 + y′²) (shortest path) and show the solution is a straight line. Apply (b) to the descent-time integrand F = √((1 + y′²)/(2gy)) and show that y(1 + y′²) = C.
(d) Verify that the cycloid x = a(φ − sin φ), y = a(1 − cos φ) satisfies y(1 + y′²) = 2a. For a = 0.5 m the cycloid reaches (πa, 2a) at φ = π in time T = π√(a/g). Evaluate it, evaluate the time along the straight line to the same point (uniform acceleration along an incline), and give the ratio.

**15.** (2) Index notation and the stress tensor.
(a) Evaluate δ_ii, δ_ij δ_jk, and ε_ijk ε_ijk in three dimensions, showing the reasoning.
(b) Using ε_ijk ε_imn = δ_jm δ_kn − δ_jn δ_km, derive the component form of a × (b × c) and hence the identity a × (b × c) = b(a·c) − c(a·b).
(c) Plane stress at a point in a wing skin is σ = [[80, 20], [20, 20]] MPa in the (x, y) axes. Write the transformation law σ′_ij = Q_ik Q_jl σ_kl for a rotation of the axes by θ and evaluate all components of σ′ for θ = 30°. Verify that the trace and the determinant are unchanged and say why (what kind of quantities they are).
(d) Give the principal stresses, the principal angle and the maximum in-plane shear stress, and say which of them is a tensor invariant and which is not.

**16.** (2) Probability and distributions: landing dispersion. A capsule's landing-point errors are independent Gaussians with σ_x = 150 m (cross-range) and σ_y = 400 m (down-range), both zero-mean.
(a) Give the probability that the cross-range error is within ±300 m.
(b) Show that the normalised radius r² = (x/σ_x)² + (y/σ_y)² has P(r² ≤ k²) = 1 − e^{−k²/2} (integrate the joint density in polar coordinates). Find k for the 95% ellipse and its semi-axes.
(c) The landing sequence relies on eight independent pyrotechnic initiators each with failure probability 0.02. Give the probability that all eight fire, and that at most one fails. Name the distribution.
(d) Micrometeoroid impacts on the heat shield arrive as a Poisson process at 0.3 per year. Over a 5-year mission, give P(no impacts) and P(two or more). Then show in two lines why the exponential waiting time between impacts is memoryless.

**17.** (2) Hypothesis testing: material acceptance. The eight coupons of Question 3 (synthetic data) have sample mean x̄ = 441 MPa and sample standard deviation s = 10 MPa; the specification is μ ≥ 450 MPa.
(a) State H₀ and H₁, compute the standard error and the t statistic, state the degrees of freedom, and give the one-sided p-value and the critical value at α = 0.05 (from a t-table). Decide.
(b) Give the 95% two-sided confidence interval for the batch mean and say what it does and does not tell you about the individual coupons.
(c) The engineer's spreadsheet used a z-test with σ = 10 MPa treated as known. Give the p-value it would report and explain why it is too small.
(d) Does the decision change at α = 0.01? Say what a Type II error would mean here and what two things reduce its probability.

**18.** (2, 3) Regression and confidence bands: a static-fire motor. Five static firings (synthetic data) of a small motor at chamber pressures p = 1.0, 1.5, 2.0, 2.5, 3.0 MPa gave mean thrusts F = 210, 330, 425, 540, 640 N. Fit F = a + b p by least squares.
(a) Compute p̄, F̄, S_xx = Σ(p − p̄)² and S_xy = Σ(p − p̄)(F − F̄) by hand and give b and a.
(b) Give the fitted values, the residuals, SSE, the residual standard deviation s = √(SSE/(n − 2)) and R².
(c) Give the standard error of b and its 95% confidence interval (t with 3 degrees of freedom).
(d) Give the 95% confidence band for the mean thrust at p = 2.5 MPa and the 95% prediction interval for a single new firing at the same pressure, and say why they differ. Then evaluate the confidence half-width at p = 4.0 MPa and say what two things make a prediction there untrustworthy.

**19.** (2, 4) Random processes: the Dryden gust spectrum. The longitudinal gust velocity has one-sided power spectral density Φ(Ω) = σ² (2L/π) / (1 + (LΩ)²) in spatial frequency Ω (rad/m), with σ = 1.5 m/s and L = 533 m.
(a) Show by integrating Φ from 0 to ∞ that the area under the spectrum is σ², and say what that means for a time history of the gust.
(b) For an aircraft at V = 100 m/s, convert to temporal frequency ω = VΩ, write Φ(ω) (keep the area equal to σ²), and give the break frequency V/L in rad/s and Hz and the value of Φ at zero frequency.
(c) Show that the fraction of the gust variance below a frequency ω_c is (2/π) arctan(Lω_c/V), and evaluate it for ω_c = 1 rad/s and ω_c = 0.1 rad/s.
(d) A 60 s flight record is used to estimate this spectrum by a single periodogram. Give the frequency resolution, say how many resolution bins lie below the break frequency, say what the variance of a single-periodogram estimate is relative to the true spectrum and why averaging segments (Welch) trades resolution for variance, and state whether this record can validate the Dryden low-frequency shape.

**20.** (2) Bayesian inference: a failure rate. A component's failures are a Poisson process with rate λ (per flight-hour). The prior is Gamma(α₀ = 2, β₀ = 1000 h) (density ∝ λ^{α₀−1} e^{−β₀λ}), with mean α₀/β₀ = 2 × 10⁻³ per hour. In service, 3 failures are recorded in 2500 flight-hours.
(a) Write the Poisson likelihood for k failures in T hours and show that the posterior is Gamma(α₀ + k, β₀ + T).
(b) Give the posterior mean, mode and standard deviation, the maximum-likelihood estimate k/T, and the prior mean, and explain the ordering of the three estimates.
(c) How many further failure-free hours would bring the posterior mean below 10⁻³ per hour?
(d) Show that the posterior predictive probability of zero failures in the next τ hours is (β₁/(β₁ + τ))^{α₁}, evaluate it for τ = 1000 h, and compare it with the plug-in value e^{−λ̂τ} using the posterior mean; say why they differ.

## Part C: Computation

**21.** (3) Write a Python function `rk4_step(f, t, y, h)` that advances the ODE ẏ = f(t, y) one step of the classical fourth-order Runge–Kutta method, with `y` a NumPy array, and a function `integrate(f, t0, y0, t_end, n_steps)` that uses it. Apply it to the vacuum-burn equation of Question 5(b), dv/dt = ṁv_e/(m₀ − ṁt) − g, from t = 0 to burnout. Write three pytest tests: one that compares v(t_b) against the closed form v_e ln(m₀/m_dry) − g t_b, one that shows fourth-order convergence (the error falls by a factor of about 16 when the step halves), and one limiting case (ṁ = 0 gives free fall v = v₀ − gt exactly). Then explain in two sentences why the tests assert a convergence order and a closed-form comparison rather than a hard-coded number, and in one sentence why `scipy.integrate.solve_ivp` is still what you would call in production.

## Part D: Interpretation

**22.** (4) The posterior of Question 20 is Gamma(5, 3500 h): mean 1.43 × 10⁻³ per hour. The requirement is λ ≤ 10⁻³ per hour. A programme manager reads the mean and says "the requirement is not met; ground the fleet". An engineer replies "the maximum-likelihood rate is 1.2 × 10⁻³, and anyway the prior is dragging the number up; with a flat prior we would be fine".
(a) Compute the posterior probability that λ exceeds 10⁻³ per hour (you may use the Gamma survival function from a table or `scipy.stats.gamma.sf`), and the 95% credible interval. Say what each statement above gets right and wrong, and what the honest one-sentence summary is.
(b) Repeat the probability with a vague prior Gamma(0.5, 0) (posterior Gamma(3.5, 2500)). Does the engineer's "we would be fine" survive? Say what the comparison tells you about how much the data versus the prior drive the conclusion at this sample size.
(c) The decision to ground a fleet is not the same as the estimate. Say what additional quantity (not in the posterior) is needed to make the decision, and give two engineering reasons the Poisson constant-rate model itself might be wrong for this component, with the observation that would reveal each.
(d) How many failure-free hours would be needed before the posterior probability that λ > 10⁻³ falls below 5%? Compare with the 1500 h from Question 20(c) and explain why the two answers differ so much.
