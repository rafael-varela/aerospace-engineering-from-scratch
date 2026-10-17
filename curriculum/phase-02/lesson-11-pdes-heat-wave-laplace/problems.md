# P02.11  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a thicker tile by hand)

The tile of the lesson (α = 4.96 × 10⁻⁷ m²/s, T₀ = 300 K, T_s = 1300 K) is made 40 + s mm thick.

(a) Give the diffusion time L²/α and the Fourier number at 1200 s.
(b) Estimate the back-face temperature at 1200 s with the one-term formula, then with `tile_temperature`. Why is the one-term estimate low rather than high?
(c) At what time does the back face reach 450 K? (Use `brentq` on `tile_temperature`, or invert the one-term formula first.)

Answer form: four numbers with units and a one-sentence explanation.

## 2. (Derivation: the mixed-boundary slab and the harmonic term)

(a) Starting from the energy balance on a slab element, derive T_t = αT_xx, naming the assumption that makes the flux −k T_x and the one that lets you divide by dx.
(b) For θ(0) = 0, θ_x(L) = 0, separate variables and show that the eigenvalues are λ_n = (2n + 1)π/2L. What goes wrong if you try λ_n = nπ/L?
(c) Compute the coefficients of the constant initial condition θ₀ against sin λ_n x and show they are 4θ₀/((2n + 1)π).
(d) Verify directly that sin(kx) sinh(ky) satisfies u_xx + u_yy = 0, and that sin(kx) sin(ky) does not. Which sign of the separation constant does Laplace's equation force, and why does the heat equation force the other?

Answer form: four short derivations.

## 3. (Unfamiliar scenario: acoustic modes of a combustion chamber)

The pressure perturbation p′ in a gas-filled tube obeys the wave equation p′_tt = c² p′_xx with c the speed of sound. A rocket combustion chamber is modelled as a tube of length L = 0.8 + s/100 m closed at both ends (zero gas velocity, which means p′_x = 0 at each end), with c = 1100 m/s in the hot gas.

(a) Separate variables with Neumann conditions at both ends and show the modes are cos(nπx/L) with f_n = nc/2L, n = 0, 1, 2, …. What is the n = 0 "mode"?
(b) Give the first three longitudinal frequencies in hertz. These are the frequencies at which combustion instability ("chugging" at low frequency, "screech" at high) can lock on.
(c) If the injector end were instead open (p′ = 0 there), what would the modes and frequencies become?
(d) Use `string_mode_frequency` with a fictitious tension and mass per length chosen to give c = 1100 m/s, and explain why the string function can be reused for the closed–closed tube but not for the closed–open one.

Answer form: a derivation, six frequencies, and a short explanation.

## 4. (Numerical, uses your code: FTCS at the stability edge)

Use the 50 mm tile with the hot face at T_s = 1100 + 10s K, 51 grid points, and march to 600 s.

(a) Run `heat_ftcs` with r = 0.25 and r = 0.5. Report the maximum error against `tile_temperature` and the back-face temperature for each. Which is more accurate, and why is the r = 0.5 result worse although both are stable?
(b) Write a ten-line marcher of your own without the guard and run r = 0.51 for 100, 500, 1000 steps, recording max |T|. Compare the growth per step with the von Neumann amplification factor |1 − 4r| at the shortest wavelength.
(c) Show, with the amplification factor, that the growing mode alternates in sign from one grid point to the next.
(d) How many steps does r = 0.25 need to reach 600 s, and how does that count change with 101 points? What does this say about explicit schemes for fine grids?

Answer form: a table (r, max error, back face), three numbers for the growth, a short argument, two step counts.

## 5. (Numerical, uses your code: where you pluck matters)

The lesson's cable (L = 1.2 m, T = 200 N, μ = 0.010 kg/m, h = 5 mm) is plucked at a = L × s/40 from one end and released.

(a) Compute b₁ … b₈ with `plucked_string_coefficients` and tabulate b_n/b₁. Which modes are missing or nearly missing, and why?
(b) Plot (or tabulate at 21 points) the shape at t = 0, L/4c, L/2c with `dalembert` and with `wave_series` (60 terms). Where is the largest difference and what causes it?
(c) Run `wave_leapfrog` at s = 0.8 to t = 10 ms and report the energy drift relative to `string_energy` of the initial shape, and the maximum error against `dalembert`.
(d) Halve the number of grid points (keep s = 0.8). Does the energy drift double, quadruple, or stay? Explain from where the energy error comes.

Answer form: a table of eight ratios, three shapes, four numbers.

## 6. (Judgment: trusting the tile model)

The tile model holds the hot face at 1300 K, uses constant α, and insulates the back.

(a) Rank the three assumptions by how much each is likely to change the back-face temperature at 1200 s, with a direction (higher or lower) for each, and a one-line justification.
(b) Silica's conductivity roughly triples between 300 and 1300 K. Using `tile_temperature` with α × 2 as a crude bound, how much does the back face change, and what does that do to the required thickness?
(c) Propose one measurement on a single tile in a furnace that would check the constant-α assumption without needing to know α beforehand.
(d) A colleague says the 61.4 mm answer "is exact because it comes from the series". Reply in three sentences.

Answer form: four short paragraphs, each with at least one number.
