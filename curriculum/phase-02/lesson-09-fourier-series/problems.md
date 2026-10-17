# P02.09  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a one-sided pulse train)

A thruster fires with duty cycle d = (20 + s)/100: the torque is 1 for 0 < t < dT and 0 for dT < t < T, period T = 1 s.

(a) Derive a₀, a_n and b_n by hand from (2) of the lesson, and show that the amplitude of harmonic n is 2|sin(πnd)|/(πn).
(b) Compute a₀ and the first six harmonics with `fourier_coefficients` (8192 samples) and compare with your formulas.
(c) Check Parseval: the mean square of the pulse train is d. How close do six harmonics come, and how close do 200?
(d) Unlike the square wave, this train has a non-zero mean and *even* harmonics. Which harmonics vanish for your d, if any, and for which d would the second harmonic vanish?

Answer form: three formulas, a table of seven numbers with the comparison, two Parseval sums, two short answers.

## 2. (Derivation: orthogonality, the square wave, and Parseval)

(a) Prove that ∫₀^T sin(mω₀t) sin(nω₀t) dt = (T/2)δ_mn for positive integers m, n, and that ∫₀^T sin(mω₀t) cos(nω₀t) dt = 0.
(b) Derive the coefficient formulas (2) from (1) using (a).
(c) Derive b_n = 4A/(πn) for odd n (and zero for even n) for the odd square wave of amplitude A.
(d) Derive Parseval's identity (5) and use it with (c) to show Σ_{odd n} 1/n² = π²/8. Deduce Σ_{all n} 1/n² = π²/6.

Answer form: four derivations.

## 3. (Unfamiliar scenario: ripple on a rectified 400 Hz bus)

An aircraft alternator delivers v(t) = V₀ |sin(2πf₀t)| with f₀ = 400 Hz after full-wave rectification. The avionics need a smooth DC supply.

(a) What is the period of v? Show that v is even and derive a₀ = 4V₀/π and a_n = −4V₀/(π(4n² − 1)) for the harmonics of *that* period.
(b) Confirm with `fourier_coefficients` (give the first six a_n and the largest |b_n|).
(c) The **ripple** is the alternating part. Give the amplitude of the lowest ripple harmonic (at what frequency?) as a fraction of the DC level, and the rms ripple as a fraction of the DC level, both from the coefficients and from Parseval with the exact mean square ½V₀².
(d) A first-order RC filter attenuates a sinusoid at ω by 1/√(1 + (ωτ)²). What time constant τ brings the lowest ripple harmonic to 1 % of DC? Why is it enough to design for the lowest harmonic alone?

Answer form: a period and a derivation, seven numbers, four fractions, a time constant with a justification.

## 4. (Numerical, uses your code: choosing a pulse period)

The thruster schedule is changed to period T = 1 + s/100 seconds; the boom mode stays at 5.0 Hz, ζ = 0.01.

(a) List the frequencies of harmonics 3, 5 and 7 and say which is nearest the mode. Compute the response coefficients with 51 harmonics, give the dominant harmonic and its amplitude, and the peak of the summed response (static deflections). Compare with the 14.33 of the 1 s pattern.
(b) Sweep T from 0.90 to 1.30 s in 0.01 s steps and tabulate the peak response at 0.90, 1.00, 1.10, 1.15 (or your T), 1.20 s. Which T in the range gives the least response?
(c) Verify your T's result by integrating the modal equation with `solve_ivp` for 40 s and comparing the last five seconds with the harmonic sum.
(d) The mode's frequency is uncertain by ± 0.2 Hz. Over what range of T is your chosen period safe from the 5th harmonic, and from the 3rd and 7th?

Answer form: three frequencies, a harmonic, two amplitudes, a five-row table and a best T, two agreement numbers, three T ranges.

## 5. (Judgment: what the series does and does not tell you)

(a) The partial sum of the square-wave torque overshoots by 9 % near each firing edge. A colleague proposes adding this "Fourier overshoot" to the thruster's rated torque. Explain why that is wrong, and what the overshoot actually is.
(b) The predicted boom response at the 1 s pattern is 14.3 static deflections, of which 12.7 come from one harmonic through a factor 1/(2ζ). Rank the uncertainties in that number: the Fourier coefficients, the mode frequency, the damping ratio, the valve rise time. Which would you measure first?
(c) The thruster valve opens in 40 ms. Estimate the harmonic beyond which the real coefficients fall as 1/n² instead of 1/n, and say whether this changes the answer to the engineering question.
(d) The attitude controller changes the pattern every 20 s. Does the boom ever reach the steady-state amplitude of the table? Give the number that decides it.

Answer form: four short paragraphs, each with a number.

## 6. (Hard: the Gibbs constant and the convergence rate)

(a) Show that the N-term partial sum of the unit square wave (N odd, harmonics 1, 3, …, N) can be written S_N(t) = (4/π) Σ_{k odd ≤ N} sin(kω₀t)/k, and that its derivative is S_N′(t) = (4ω₀/π) Σ_{k odd ≤ N} cos(kω₀t) = (2ω₀/π) sin((N + 1)ω₀t)/sin(ω₀t). (Use the geometric-series formula of P02.07, problem 6.)
(b) Hence show that the first maximum of S_N after t = 0 is at t* = T/(2(N + 1)), and check it with `gibbs_overshoot`'s grid for N = 9, 29, 99.
(c) Integrate S_N′ from 0 to t* with the substitution u = (N + 1)ω₀t and the approximation sin(ω₀t) ≈ ω₀t for large N, and show S_N(t*) → (2/π)∫₀^π (sin u/u) du = (2/π) Si(π) = 1.17898. Evaluate Si(π) numerically (`scipy.special.sici` or `quad`).
(d) From Parseval's tail, derive the leading-order rms truncation error 2/(π√N) for the square wave and compare with `truncation_error_rms` at N = 21 and N = 81. Why does the Gibbs overshoot not contribute to this error in the limit?

Answer form: two derivations, three peak positions, one constant, two error comparisons with an explanation.
