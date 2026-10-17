# P02.03  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Compute ω_n and ζ first, then r = ω/ω_n; everything in (3) is a function of r and ζ alone. The static deflection F₀/k turns M into metres.
2. In (a) remember ẍ_p = −ω²x_p, which is why only two coefficients appear. In (b) multiply top and bottom of 1/(a + ib) by (a − ib). In (c) treat u = r² as the variable; the denominator is a quadratic in u. In (d) the t-terms cancel because ω_n² = k/m, so only the terms that came from differentiating t once survive. In (e) use sin ε ≈ ε and ω_n² − ω² = (ω_n − ω)(ω_n + ω).
3. An accelerometer wants M = 1 (quasi-static), an isolator wants M ≪ 1 (mass-controlled); both are statements about r. For (c) look at the sign of the r² coefficient in M⁻²: when it is negative M first rises above 1, and a small rise can stay inside the tolerance while delaying the fall.
4. `brentq(lambda f: magnification(2*np.pi*f, m, c, k)*F0/k - X_lim, f_a, f_b)` with one bracket below the peak and one above; check first that the peak exceeds the limit at all.
5. For (a) the amplitude is m_e e ω²|H|; divide top and bottom by mω² and let ω → ∞. For (b) compare 2 s for 20 Hz (10 Hz/s) with the lesson's sweep table. For (d) think about what a *phase* measurement shows that an amplitude one does not.
6. The instantaneous frequency is the derivative of the phase divided by 2π. The amplitude cannot build faster than 1 − e^{−t/τ}; the question is how long the drive lingers inside the half-power band at each rate.
