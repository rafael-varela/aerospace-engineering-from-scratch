# P02.03  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. ω_n = √(k/m), ζ = c/(2√(km)), r = 2πf/ω_n. 2. M = 1/√((1 − r²)² + (2ζr)²); X = M F₀/k; φ = atan2(2ζr, 1 − r²). 3. A = F₀(k − mω²)/D, B = F₀cω/D with D = (k − mω²)² + (cω)²; check √(A² + B²) = X and atan2(B, A) = φ. 4. 1/(2ζ); f_n√(1 − 2ζ²); Δf ≈ 2ζf_n (or `half_power_bandwidth`); Q = f_peak/Δf.

## 2.
(a) Substitute; collect cos ωt and sin ωt; two equations in A, B; solve by Cramer's rule. (b) Z = F₀/(k − mω² + icω); x_p = Re[Z e^{iωt}] = Re Z cos ωt − Im Z sin ωt; identify A and B. (c) Let u = r²: D(u) = (1 − u)² + 4ζ²u; dD/du = −2(1 − u) + 4ζ² = 0 ⟹ u = 1 − 2ζ²; substitute back; the peak exists only if u > 0. (d) With x_p = t(A cos + B sin): ẍ_p = −ω_n² x_p + 2ω_n(B cos − A sin); m ẍ_p + k x_p = 2mω_n(B cos − A sin); match to F₀ cos. (e) Replace sin((ω_n − ω)t/2) by (ω_n − ω)t/2 and ω_n² − ω² by 2ω_n(ω_n − ω); cancel.

## 3.
(a) M ≈ 1 needs r ≪ 1 (f_n ≫ f); M ≪ 1 needs r ≫ 1 (f_n ≪ f). (b) Scan r on a fine grid, find the first r where |M − 1| > 0.05, f_n ≥ f_max/r. (c) Repeat for each ζ; tabulate r_max; the r² coefficient 4ζ² − 2 is negative for ζ < 1/√2, so M overshoots 1 first; the overshoot stays under 5 % down to about ζ = 0.6, which pushes the 0.95 crossing to larger r. (d) Sensitivity: output ∝ 1/ω_n², so a tenfold f_n costs a hundredfold in signal.

## 4.
(a) `X = lambda f: magnification(2*pi*f, m, c, k)*F0/k`; confirm `peak_magnification*F0/k > X_lim`; `brentq(X - X_lim)` on [20, 55] and [55, 120]. (b) Same with c for ζ = 0.05 and 0.10; the limit is never exceeded once M_peak F₀/k < X_lim, i.e. ζ > F₀/(2kX_lim) approximately. (c) The half-power band is where M² > ½M_peak²; the tolerance band is where M > X_lim k/F₀; they coincide only if X_lim happens to be the −3 dB level. (d) `full_response(t, wn, F0, m, c, k, 0, 0)`, first t with |x| > X_lim; compare with −τ ln(1 − X_lim/X_ss).

## 5.
(a) X = m_e e ω²/√((k − mω²)² + (cω)²) → m_e e/m as ω → ∞; stiffening moves the resonance up but the high-frequency floor is set by the mass ratio, not the stiffness. (b) 10 Hz/s is between the lesson's 6.7 and 20 Hz/s rows (91 % and 78 % of the peak): the measured peak is low and wide, so Q is biased low; require ≲ ζ²f_n² ≈ 1 Hz/s. (c) Joint friction (Coulomb: equivalent ζ falls as amplitude rises); material/viscoelastic damping (roughly amplitude-independent); a two-level test shows Q changing for the first and constant for the second. (d) The phase of the response against frequency: a second mode produces a second 180° swing (or a dip and recovery) inside the band.

## 6.
(a) Phase Φ(t) = 2π(45t + 10t²/T); f_inst = Φ′/2π = 45 + 20t/T; rate = 20/T. (b) `rk4_forced(m, c, k, force, 0, 0, T, int(10000*T))`; `argmax(abs(x))`; `f_inst` at that time; divide by F₀/(cω_n). (c) Δt = 2ζf_n/rate; compare 1 − e^{−Δt/τ}; note it over-predicts slightly for fast sweeps because the band is only half traversed when the peak occurs. (d) The response lags the build-up, so on an upward sweep the peak appears above f_n; a downward sweep puts it below. Run the chirp with 65 − 20t/T.
