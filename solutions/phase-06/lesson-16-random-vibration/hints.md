# P06.16  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. +6 dB/oct is G ∝ f²: over the 1.32 octaves from 20 to 50 Hz, G changes by (50/20)². On a log-log segment the trapezoid rule on linear axes replaces a curve that sags below the chord by the chord itself.
2. (b) Under r → 1/r, D(1/r) = D(r)/r⁴ and dr → −dr/r², and the limits 0 and ∞ swap. (c) Parseval: ∫|H(ω)|² dω over the whole axis equals 2π ∫₀^∞ h(t)² dt, and the impulse response of 1/(ωₙ² − ω² + 2iζωₙω) is e^{−ζωₙt} sin(ω_d t)/ω_d.
3. Compare G just below fₙ with G(fₙ) on a falling segment, and remember |H|² ≈ 1 over the whole plateau below the resonance.
4. Each normalised deviation is (ratio − 1)/u with u = 1/(2√(ζωₙT)). For the transfer function, fit |H| with the lesson's `transmissibility_sq` inside the model; give `curve_fit` a starting guess near 150 Hz.
5. On a flat input the response goes as √(Q G) and as √fₙ; turn 7.2/9.66 into a factor on G, on Q and on fₙ. Then ask which of those can change between a design estimate and a test.
6. A narrow-band process has one maximum per zero crossing. The absolute acceleration contains the input passed through at high frequency, which adds small wiggles; those are maxima too.
