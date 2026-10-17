# P08.13  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. A 3.3 s record holds about 650 integral time scales; the relative scatter of its variance is √(2 × 0.866 L/U / T_rec). For the spectral plateau, ask where the first three bins of a 2¹²-point segment sit on the model spectrum.
2. (a) The average of a product of two fluctuations is not zero. (b) Write σ̂² − σ² as an integral of u′²(t) − σ² and use the Gaussian fourth-moment identity ⟨a²b²⟩ = ⟨a²⟩⟨b²⟩ + 2⟨ab⟩². (c) Solve for the exponents a, b in ν^a ε^b. (d) S(f) and R(τ) are a Fourier pair (Wiener–Khinchin).
3. Everything scales with U/L for the big eddies and with (ν³/ε)^(1/4) for the small ones; the sample rate must exceed twice the highest frequency you want, with margin.
4. A standard deviation estimated from n samples is itself uncertain by about 1/√(2(n − 1)). Compare the spectral plateau's frequencies with the corner of the model spectrum at 1.339kL = 1.
5. A moving average over w samples multiplies the spectrum by |sin(πfw/f_s)/(w sin(πf/f_s))|²; decimation folds the spectrum above the new Nyquist frequency back onto the band below it.
6. Think of the tails of the probability distribution, the derivative skewness, and the direction of energy transfer; for the database, think of what "simulation" means.
