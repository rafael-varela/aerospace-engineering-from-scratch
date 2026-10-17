# P02.08  Laplace Transform

> Differentiation becomes multiplication by s, initial conditions become part of the input, and an ODE becomes a fraction. Where the fraction blows up is where the system lives.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the Laplace transforms of the standard signals and the derivative, shift and s-shift rules, and prove each numerically by quadrature. 2. Implement and test a transfer-function solver: poles, residues by the cover-up rule, the inverse transform as a sum of exponentials, and step, impulse, delayed-step and initial-condition responses, against `scipy.signal.lti`, `solve_ivp` and SymPy as oracles. 3. Measure the overshoot, peak time and settling time of an elevator actuator from its response and compare with the pole-location formulas. 4. Decide whether an actuator meets a step-response specification and locate the pole positions that would. |
| **Prerequisites** | `P02.02` (second-order ODEs, damping ratio), `P02.03` (forced response), `P02.07` (complex numbers, the complex gain), `P01.04` (integration by parts, improper integrals and partial fractions), `P00.09` (SymPy). Skill check: you can integrate ∫₀^∞ e^{−st} dt, do a partial-fraction split with two linear factors, and find the roots of a quadratic. |
| **Engineering question** | An electro-hydraulic elevator actuator behaves as y″ + 30y′ + 625y = 625u (ω_n = 25 rad/s, ζ = 0.6); commanded a 1° step, how far does it overshoot, when does it peak, when is it within 2 % of the command, and what does it do to a 100 ms command pulse? |
| **Assumptions and limits** | Linear, time-invariant actuator with constant coefficients; inputs that are zero before t = 0 and grow no faster than an exponential; simple (non-repeated) poles in the residue code, repeated poles left to SymPy; the transform is one-sided, so the method knows nothing about t < 0 and nothing about a system whose coefficients change with time or amplitude (a rate-limited valve, hydraulic stiffness that varies with oil temperature). The 2 % settling time is read from a sampled response and depends on the grid. |
| **Analytical work** | The transform of a step, ramp, exponential, sinusoid and damped sinusoid; the derivative rule by integration by parts; the shift and s-shift rules; the transfer function of the actuator; poles and the characteristic polynomial; partial fractions by the cover-up rule and the residue formula N(p)/D′(p); the inverse transform of the step response; overshoot and peak time from the pole positions; the impulse response as the derivative of the step response and the convolution integral; the frequency response as G(iω). |
| **Computational deliverable** | `starter/main.py`: `laplace_numeric`, `transform_table`, `transfer_function`, `evaluate`, `poles`, `is_stable`, `residues`, `inverse_laplace`, `response_with_initial_conditions`, `step_response`, `step_response_closed_form`, `impulse_response`, `delayed_step_response`, `actuator_rhs`, `convolve_causal`, `overshoot`, `peak_time`, `settling_time`, `frequency_response`, `laplace_symbolic_step_response`. Tests prove: quadrature reproduces every table entry to 10⁻⁷; the linearity, derivative, shift and s-shift rules hold numerically; the residue expansion reproduces N/D at random complex points and raises on improper or repeated-pole inputs; the inverse transform is real, causal, reproduces the table and round-trips through quadrature; the step response from residues matches the closed form to 10⁻¹², `scipy.signal.lti` to 10⁻⁹ and `solve_ivp` to 10⁻⁸; initial conditions enter correctly; the impulse response is the derivative of the step response and integrates to G(0); a pulse is the difference of two steps; overshoot and peak time match the pole formulas; G(iω) equals the complex gain of P02.07; SymPy agrees. |
| **Practical activity** | **computational experiment**: a sweep of overshoot and settling time against damping ratio and natural frequency, and of the peak pulse response against pulse length (a shock-response curve). |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `laplace_transform.png`, the hand derivation of the step response by partial fractions, the sweep tables, and a short note stating whether the actuator meets a 10 % overshoot / 250 ms settling specification and what pole positions would. |
| **Assessment** | `assessments/phase-02/grading-rubric.md`; competencies 1–4. Pass: the derivative rule and the step response are derived by hand; all public tests pass; the hand overshoot and peak time agree with the code to three figures; the specification note names the pole region that satisfies it. |
| **References** | `@mit-ocw-18-03-differential-equations` L19 "Introduction to the Laplace Transform", L20 "Derivative Formulas", L22 "Using Laplace Transform to Solve ODEs with Discontinuous Inputs", L23 "Use with Impulse Inputs", notes 26 "Laplace transform: basic properties", 27 "Application to ODEs", 28 "Second order equations; completing the squares", 29 "The pole diagram", 30 "The transfer function and frequency response" (required-free); L21 "Convolution Formula", notes 23 "Step functions and delta functions", 24 "Step response, impulse response", 25 "Convolution" (supplementary-free); `@openstax-calculus-volume-2` §3.4 "Partial Fractions" (required-free, the inverse-transform algebra). |
| **Extensions** | ★ Extend `residues` to repeated poles with the derivative formula and verify on 1/(s + 1)². ★ Close a proportional loop around the actuator, y″ + 30y′ + 625y = 625K(r − y), and find the gain K at which the overshoot reaches 20 %. ★ Add a pure transport delay e^{−τs} to the loop and show with `delayed_step_response` why a delay destabilises feedback. |
| **Time** | unmeasured |

## 1. The engineering problem

The elevator of a fly-by-wire aircraft is moved by an electro-hydraulic actuator. The flight-control computer commands an angle u(t); the surface reaches y(t) a little later, with some overshoot. Tested on the bench, the actuator behaves to good approximation as a second-order system,

y″ + 2ζω_n y′ + ω_n² y = ω_n² u,   ω_n = 25 rad/s, ζ = 0.6,   (1)

with unit static gain (a constant command is eventually reached exactly). The control engineer must answer four questions before the actuator is accepted: for a 1° step command, how far does the surface overshoot, when does it peak, when is it within 2 % of the command, and how much of a 100 ms pulse (the shortest command the computer will send) actually arrives at the surface? The specification says less than 10 % overshoot and settled within 250 ms.

P02.02 and P02.03 can answer these one at a time: solve the homogeneous equation, add a particular solution, fit the initial conditions, repeat for each input. The **Laplace transform** does all of them with one algebraic object, the **transfer function** G(s), and it handles the inputs that P02.03 could not: steps that switch on at t = 0.2 s, pulses, impulses, and any input whose transform is a rational function. It also makes the stability question a geometric one: the system is stable exactly when every **pole** of G lies in the left half of the complex plane, and the overshoot and settling time are read directly from where the poles sit. This is the mathematics of every control-system course that follows.

Success criteria, written down first: every entry of our transform table is reproduced by numerical quadrature to 10⁻⁷; the step response built from residues agrees with the textbook closed form to 10⁻¹², with `scipy.signal.lti` to 10⁻⁹ and with a direct integration of (1) to 10⁻⁸; the overshoot and peak time read from the response agree with the pole-location formulas to 10⁻⁴; and the frequency response G(iω) reproduces the complex gain of P02.07 exactly.

## 2. Understand and derive

### 2.1 The transform and the table

For a function f(t) defined for t ≥ 0, the **Laplace transform** is

**F(s) = L{f}(s) = ∫₀^∞ f(t) e^{−st} dt**,   (2)

a function of the complex variable s wherever the integral converges (`@mit-ocw-18-03-differential-equations` notes 26). For every f that grows no faster than e^{at} it converges for Re s > a. Units: if t is in seconds, s is in s⁻¹, and F has the units of f times seconds. Five integrals give the table that does most of the work:

| f(t), t ≥ 0 | F(s) | region |
|---|---|---|
| 1 (unit step) | 1/s | Re s > 0 |
| t | 1/s² | Re s > 0 |
| e^{−at} | 1/(s + a) | Re s > −a |
| cos ωt, sin ωt | s/(s² + ω²), ω/(s² + ω²) | Re s > 0 |
| e^{−at} cos ωt, e^{−at} sin ωt | (s + a)/((s + a)² + ω²), ω/((s + a)² + ω²) | Re s > −a |

The first: ∫₀^∞ e^{−st} dt = [−e^{−st}/s]₀^∞ = 1/s. The third is the same integral with s + a. The sinusoids come from cos ωt = Re e^{iωt} and the third line with a = −iω: 1/(s − iω) = (s + iω)/(s² + ω²), whose real and imaginary parts are the two entries. The last line is the **s-shift rule**: multiplying f by e^{−at} replaces s by s + a in F, because e^{−at}e^{−st} = e^{−(s+a)t}. The code's `laplace_numeric` evaluates (2) with `scipy.integrate.quad` for real s and the test checks every row to 10⁻⁷; the table is not something to memorise but something to verify once and then use.

### 2.2 The derivative rule: the reason the transform exists

Integrate by parts (P01.04) with u = e^{−st}, dv = f′ dt:

L{f′} = [f e^{−st}]₀^∞ + s ∫₀^∞ f e^{−st} dt = **s F(s) − f(0)**,   (3)

the boundary term at ∞ vanishing in the region of convergence. Apply it twice: **L{f″} = s² F(s) − s f(0) − f′(0)**. Differentiation has become multiplication by s, and the initial conditions appear as polynomial terms. With linearity (L is an integral, so L{af + bg} = aF + bG), a constant-coefficient ODE becomes an algebraic equation in F.

Two more rules, both one-liners from (2). **Shift in t** (notes 23): if u_c(t) is the unit step at t = c, then L{f(t − c) u_c(t)} = e^{−cs} F(s), by the substitution t′ = t − c. A delay multiplies the transform by e^{−cs}. **Convolution** (notes 25): L{∫₀^t g(t − τ) u(τ) dτ} = G(s) U(s). Products in s are convolutions in t, which is what makes the transfer function mean something in the time domain. The test checks (3), the shift rule and the s-shift rule numerically.

### 2.3 The transfer function and its poles

Transform (1) with zero initial conditions: (s² + 2ζω_n s + ω_n²) Y(s) = ω_n² U(s), so

**Y(s) = G(s) U(s)**,   **G(s) = ω_n²/(s² + 2ζω_n s + ω_n²)** = 625/(s² + 30s + 625).   (4)

G is the **transfer function**: output transform over input transform (notes 30). It is the ODE with the time derivatives replaced by powers of s, so its denominator is the characteristic polynomial of P02.02 and its roots, the **poles**, are the exponents of the free response: p = −ζω_n ± iω_n√(1 − ζ²) = **−15 ± 20i s⁻¹** for the actuator. Every output of the system contains terms e^{pt}, so the system is **stable** (every free motion decays) exactly when every pole has Re p < 0: the left half of the **pole diagram** (notes 29). The real part −ζω_n = −15 s⁻¹ is the decay rate; the imaginary part 20 rad/s is the damped frequency. Move the poles left and the response is faster; move them toward the real axis (larger ζ) and it oscillates less. The figure in `main()` shows the poles of (4) for ζ = 0.2, 0.6, 1, 2: a semicircle of radius ω_n for ζ < 1, then two real poles splitting along the axis.

With non-zero initial conditions the derivative rule adds terms: (s² + 2ζω_n s + ω_n²) Y = ω_n² U + (s + 2ζω_n) y(0) + y′(0). So

**Y(s) = [ω_n² U(s) + (s + 2ζω_n) y₀ + ẏ₀] / (s² + 2ζω_n s + ω_n²)**.   (5)

The initial state is just another input, with the same denominator: that is `response_with_initial_conditions`, and the test checks it against `solve_ivp` started from (0.4, −3).

### 2.4 Inverting: partial fractions and residues

Every Y we meet is a rational function N(s)/D(s) with deg N < deg D. Factor D by its roots p_k and, if they are distinct, write (`@openstax-calculus-volume-2` §3.4; notes 27)

$$
\frac{N(s)}{D(s)} = \sum_k \frac{r_k}{s - p_k}, \qquad r_k = \lim_{s\to p_k}(s - p_k)\frac{N(s)}{D(s)} = \frac{N(p_k)}{D'(p_k)} .
$$

The first form of r_k is the **cover-up rule** (cover the factor s − p_k and evaluate the rest at p_k); the second follows by l'Hôpital and is what the code computes, since D′ is `np.polyder`. Each term inverts by the table: **y(t) = Σ r_k e^{p_k t}**. Complex poles come in conjugate pairs with conjugate residues, so their two terms add to a real damped sinusoid 2|r| e^{−at} cos(ωt + arg r).

The step response is the case U = 1/s:

Y(s) = ω_n²/[s(s² + 2ζω_n s + ω_n²)],   poles 0, p, p̄,   residues 1, r, r̄ with r = ω_n²/[p(p − p̄)]; for the actuator **r = −0.5 − 0.375i**.

So y(t) = 1 + 2 Re(r e^{pt}) = 1 − e^{−15t}(cos 20t + 0.75 sin 20t). The general form, which `step_response_closed_form` holds for comparison, is

**y(t) = 1 − e^{−ζω_n t}[cos ω_d t + (ζ/√(1 − ζ²)) sin ω_d t]**,   ω_d = ω_n√(1 − ζ²).   (6)

The residue arithmetic and (6) agree to 4 × 10⁻¹⁶. SymPy's `inverse_laplace_transform` returns the same expression, (e^{15t} − ¾ sin 20t − cos 20t) e^{−15t}, and the test evaluates it at fifty points.

A repeated pole (D with a double root) breaks the simple-pole formula; the residue of the second-order term needs a derivative. `residues` detects the case and raises `ValueError` rather than returning nonsense; problem 3 meets one and sends it to SymPy or to hand algebra.

### 2.5 Reading the response from the poles

Differentiate (6): ẏ = e^{−ζω_n t} ω_n sin(ω_d t)/√(1 − ζ²), which vanishes first at ω_d t = π. So the **peak time** is t_p = π/ω_d = π/20 = **157.1 ms**, and inserting it in (6) gives the **overshoot**

**M_p = e^{−πζ/√(1 − ζ²)}** = e^{−0.6π/0.8} = **9.48 %**.   (7)

It depends on ζ alone: the pole *angle*. The envelope e^{−ζω_n t} reaches 2 % at t = ln 50/(ζω_n) = 3.9/(ζω_n), which is the origin of the rule of thumb **t_s ≈ 4/(ζω_n)** = 267 ms; it depends on the pole's *real part* alone. The sampled response actually enters the 2 % band for good at 238 ms, earlier than the rule, because the rule bounds the envelope and the oscillation inside the envelope happens to be near a zero crossing when the envelope hits 2 %. The actuator meets both halves of the specification (9.5 % < 10 %, 238 ms < 250 ms) with little margin on overshoot: a ζ of 0.59 would fail it. Problem 4 and the sweep in part 4 map the margins.

### 2.6 Impulse, pulse and delay

The **impulse response** g(t) is the inverse transform of G itself (U = 1, the transform of the delta function, L23): g = L⁻¹{G} = (ω_n/√(1 − ζ²)) e^{−ζω_n t} sin ω_d t. Since G = s · (G/s), the derivative rule gives g = d/dt (step response), which the test checks by finite differences, and ∫₀^∞ g dt = G(0) = 1, the static gain. The convolution rule then says that *any* response is y(t) = ∫₀^t g(t − τ) u(τ) dτ: the output is the input's past, weighted by the impulse response. `convolve_causal` computes this sum on a grid and reproduces the step response to 6 × 10⁻⁴ with a first-order rule at 0.1 ms.

A **delayed step** u_c(t) has transform e^{−cs}/s, so its response is the step response shifted: y(t − c) u_c(t). A **pulse** of length c is a step minus a delayed step, and its response is the difference of the two shifted step responses. For the 100 ms command pulse the surface reaches 0.969 of the commanded angle, 9 ms after the pulse has already ended, and the second step then pulls it back; a 50 ms pulse reaches only 0.67 (the test integrates (1) with a pulse as input and checks the superposition to 10⁻⁶). This is the practical meaning of bandwidth: commands shorter than about 1/ω_n never arrive in full.

### 2.7 The frequency response: the previous lesson in one line

Put s = iω in (4): G(iω) = ω_n²/(ω_n² − ω² + 2iζω_nω) = 1/(1 − r² + 2iζr) with r = ω/ω_n, which is the complex gain H(r) of P02.07. The transfer function evaluated on the imaginary axis *is* the frequency response: the steady-state amplitude ratio and phase lag for a sinusoid at ω. At ω = ω_n the actuator has |G| = 1/(2ζ) = 0.833 and lags by 90°; at 50 rad/s, 0.26 and 141°. The time-domain test drives (1) with cos 20t and reads an amplitude of |G(20i)| after the transient. Poles near the imaginary axis are what make |G(iω)| peak (ζ = 0.2 in the figure); poles far from it give a flat, slow response.

### 2.8 Where the model stops being valid

- **Linearity.** A hydraulic valve has a maximum flow, so the actuator has a maximum rate (typically 40–60°/s). A 1° step at 25 rad/s needs rates up to ω_n × 1° × 0.5 ≈ 12°/s, fine; a 10° step would be rate-limited and the response would be a ramp, not (6). No transfer function describes that.
- **Time invariance.** Hydraulic stiffness depends on oil temperature and trapped air; the coefficients drift during a flight. The transform assumes they are constant.
- **Simple poles.** Critical damping (ζ = 1) gives a double pole; the code hands it to SymPy.
- **One-sidedness.** The transform starts at t = 0 and encodes the past in two numbers, y(0) and ẏ(0). A system with memory beyond its state (hysteresis, backlash) is outside its reach.

## 3. Build the computational model

Open `starter/main.py`. `elevator_actuator()` holds ω_n and ζ. Polynomials are coefficient arrays in descending powers, the NumPy and SciPy convention.

### Step 1: the transform by quadrature

```python
def laplace_numeric(f, s, t_max=None):
    if s <= 0.0: raise ValueError(...)
    t_max = 60.0 / s if t_max is None else t_max
    g = lambda t: float(f(t)) * np.exp(-s * t)
    head, _ = quad(g, 0.0, t_max, limit=400)
    tail, _ = quad(g, t_max, np.inf, limit=200)
    return head + tail
```

Splitting at 60/s matters: `quad` on [0, ∞) alone misses a sinusoid's oscillations. `transform_table` returns the closed forms of 2.1.

### Step 2: poles and residues

```python
def residues(num, den):
    p = poles(den)
    r = np.polyval(num, p) / np.polyval(np.polyder(den), p)      # N(p) / D'(p)
    return p, r

def inverse_laplace(num, den, t):
    p, r = residues(num, den)
    f = sum(rk * np.exp(pk * t) for pk, rk in zip(p, r))
    return np.where(t >= 0.0, f, 0.0).real                       # after checking the imaginary part is rounding
```

Check that deg N < deg D and that no two poles coincide (relative to the largest pole modulus) and raise `ValueError` otherwise; check that the imaginary part is at rounding level before dropping it.

### Step 3: responses

`transfer_function` returns (4); `response_with_initial_conditions` builds (5) with `np.polymul` and `np.polyadd`, taking U as its own (num, den) pair so that a step (1, [1, 0]) and a sinusoid ([ω], [1, 0, ω²]) are handled alike; `step_response`, `impulse_response` and `delayed_step_response` are two-line wrappers. `overshoot`, `peak_time` are (7); `settling_time` finds the last sample outside the band and returns the next time.

### Step 4: the oracles

`actuator_rhs` is (1) in first-order form for `solve_ivp` with an arbitrary input function; `frequency_response` is `evaluate(num, den, 1j * omega)`; `laplace_symbolic_step_response` calls `sympy.inverse_laplace_transform` on G/s.

### Verify

Run `python -m pytest curriculum/phase-02/lesson-08-laplace-transform`. The tests:

1. **The table and the rules.** Eight table entries by quadrature to 10⁻⁷; linearity, the derivative rule with f(0) ≠ 0, the shift rule with c = 0.4 and the s-shift rule.
2. **Algebra.** Unit static gain; the poles equal −ζω_n ± iω_d; negative or zero damping is reported unstable; an overdamped pair has product ω_n² and sum −4ω_n; the cover-up rule on 1/((s + 1)(s + 3)); a cubic's residue expansion reproduces N/D at five random complex points; improper and repeated-pole inputs raise.
3. **Inversion.** Three table entries recovered to 10⁻¹², real dtype, zero for t < 0, and a quadrature round trip.
4. **Agreement.** Step response against the closed form (10⁻¹²), `lti.step` (10⁻⁹), `solve_ivp` (10⁻⁸); initial conditions against `solve_ivp`; impulse response against `np.gradient` of the step, `lti.impulse`, its integral G(0) = 1, and convolution with a step; a pulse against `solve_ivp` with a discontinuous input.
5. **Formulas.** Overshoot and peak time from the response match (7) to 10⁻⁴ and 10⁻³ for ζ = 0.2, 0.6, 0.9; the settling time lies inside 1.2 × 4/(ζω_n); M_p(0) = 1.
6. **Frequency response.** G(iω) equals the P02.07 complex gain; a cos 20t drive gives the amplitude |G(20i)| after 1 s to 10⁻³; SymPy's inverse transform agrees at fifty points.

### Compare

`scipy.signal.lti(num, den).step()` and `.impulse()` are the library versions of Steps 2–3. They convert the transfer function to state-space form and use the matrix exponential of P02.05 rather than residues, which is why they also handle repeated poles and why they are what you call in practice. Our residue code agrees with them to 4 × 10⁻¹⁴ on the step and 10⁻⁸ ω_n on the impulse; the difference is rounding in `np.roots`. We wrote the residue version because it shows *where each term of the response comes from*: one exponential per pole, with the pole's real part as its decay rate and its imaginary part as its frequency, which is the picture every later control lesson uses. `solve_ivp` is the third oracle and the only one that also handles the rate-limited valve of 2.8.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Three sweeps with `main.py`; all records are synthetic.

1. **Damping sweep.** For ζ = 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 at ω_n = 25 rad/s: overshoot from (7) and from the sampled response, peak time, 2 % settling time. **Prediction to write down first:** overshoot falls monotonically; settling time is *not* monotone in ζ (the envelope rule says it falls, but the oscillation's phase at the band edge makes it jump).
2. **Frequency sweep.** For ω_n = 15, 20, 25, 30, 40 rad/s at ζ = 0.6: peak time and settling time. Predict both to scale as 1/ω_n exactly.
3. **Pulse sweep (a shock-response curve).** Half-sine command pulses of length 25, 50, 100, 200, 400 ms: the peak surface deflection as a fraction of the pulse amplitude. Predict first where the peak reaches 0.9 (a pulse of about one natural period, 250 ms) and the behaviour for short pulses (peak ∝ pulse area, the impulse limit).

Raw data format: three CSV tables, `zeta, Mp_formula, Mp_sampled, tp, ts`; `wn, tp, ts`; `pulse_ms, peak`, each with a header naming the fixed parameters.

What you should find: overshoot 37.2 % at ζ = 0.3 down to 1.5 % at 0.8; the settling time is 449, 336, 323, 238, 239, 150 ms for ζ = 0.3 … 0.8, against the rule's 533, 400, 320, 267, 229, 200 ms: it is not monotone, it can exceed the rule (ζ = 0.5, 0.7) and it can beat it by 25 % (ζ = 0.8); every time in sweep 2 scales as 1/ω_n to all digits; the pulse peaks are 0.197, 0.383, 0.686, 0.960, 1.028 for ζ = 0.6, with the 25 ms value within 1 % of (pulse area) × max g = 0.0159 × 12.5 = 0.198.

## 5. Validate and improve

| quantity (actuator) | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| poles | −15 ± 20i | −15 ± 20i | 0 | exact quadratic |
| residue of G/s at p | −0.5 − 0.375i | −0.5 − 0.375i | 10⁻¹⁶ | — |
| step response vs closed form | (6) | residue sum | 4 × 10⁻¹⁶ | rounding |
| overshoot | 9.48 % | 9.48 % at 157.1 ms (grid) | < 10⁻⁴ | 0.1 ms grid |
| settling time | 4/(ζω_n) = 267 ms (envelope) | 238 ms (sampled) | 11 % | rule bounds the envelope, not the response |
| ∫ g dt | 1 | 1.00069 on [0, 0.5 s] | 7 × 10⁻⁴ | truncation: e^{−7.5} tail |
| G(25i) | 1/(2ζ) at −90° = 0.833 | 0.8333 at −90.0° | 0 | — |

Each discrepancy is named and small; the only one of engineering size is the settling time, and it is the rule of thumb that is approximate, not the solution. What the agreement does not show: whether the bench actuator is linear. A real step test at 1°, 3° and 10° is the check: if the three normalised responses coincide, (1) is a model; if the 10° response is a straight ramp, the valve is saturating and the specification must be met with the nonlinear model (extension of P02.06's thinking: a rate limit is a nonlinearity whose effect depends on amplitude).

Improvements worth making: repeated-pole residues (extension); a `pulse_response` that takes an arbitrary piecewise-constant command and superposes shifted steps; closing a loop around G and tracking the poles as the gain rises (the root locus of a later phase). Improvement worth not making: trusting (7) for a rate-limited step.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivation of (3) and of the step response (6) by partial fractions; `laplace_transform.png`; the three CSV tables; and a one-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) stating whether the actuator meets the 10 % / 250 ms specification, the margin, and the region of the pole diagram that would meet it.

### Adversarial self-checks

- **What if ω_n doubles?** The poles move to −30 ± 40i: same angle, same overshoot, every time halved. What if ζ → 1? The poles meet at −25: no overshoot, double pole, `residues` raises; what if ζ → 0? Poles on the imaginary axis, no decay, G(iω_n) infinite.
- **What if the input is a ramp instead of a step?** U = 1/s², a double pole at the origin: the steady-state error for a ramp command of rate R is 2ζR/ω_n (problem 3 does the first-order version). Derive it from the final-value theorem lim_{s→0} s Y(s).
- **Which assumption fails first in the real world?** Linearity: the valve's rate limit, at commands of a few degrees.
- **What measurement could falsify the model?** Step responses at three amplitudes that do not scale; or a measured frequency response whose phase lag at ω_n is not 90°.
- **Does the result survive a change of solver or grid?** Residues, closed form, `lti` and `solve_ivp` agree to 10⁻⁸; the overshoot is grid-independent to 10⁻⁴; the settling time moves by one grid step.
- **Which decision changes if the uncertainty doubles?** If ζ is known as 0.6 ± 0.05 the overshoot is between 7.4 % and 12 % and the specification is not demonstrably met; the decision would be to measure ζ from the log decrement of the bench step (P02.06, problem 4) before signing off.

## Check yourself

1. Write the definition of the Laplace transform and compute L{e^{−3t}} from it.
2. State the derivative rule for f′ and f″ and say where the initial conditions go.
3. What is L{t e^{−3t}}? (Use the s-shift rule on 1/s².)
4. The denominator s² + 4s + 29: where are the poles, is the system stable, and what are its damped frequency and decay rate?
5. For ζ = 0.5, what is the step-response overshoot, and does it depend on ω_n?
6. Give the 2 % settling-time rule for ω_n = 25 rad/s, ζ = 0.6, and say why the real response may settle sooner.
7. What is the transform of a step that starts at t = 2 s, and what is the response of G to it in terms of the ordinary step response?
8. G(iω) at ω = ω_n for the actuator is −0.833i. What are the amplitude ratio and phase lag, and which lesson's formula is this?
