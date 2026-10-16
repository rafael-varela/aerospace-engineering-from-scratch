# P02.11  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for $s = 15$; numerical results are from `solution.py` and given to the figures shown.

## 1. A thicker tile by hand ($s = 15$: L = 55 mm)

**Derivation.** (a) L²/α = 0.055²/4.96 × 10⁻⁷ = **6098 s**; Fo = 1200/6098 = **0.1968**.

(b) One term: T_back = 1300 + (300 − 1300)(4/π)e^{−(π²/4)(0.1968)} = 1300 − 1000 × 1.2732 × 0.6154 = **516.5 K**. `tile_temperature` (50 terms): **521.8 K**. The one-term estimate is *low* because the n = 1 term, 4θ₀/(3π) sin(3π/2) e^{−9(π²/4)Fo} = −1000 × 0.4244 × (−1) × 0.0126 = +5.4 K at the back face, is positive there (sin λ₁L = −1 times the negative θ₀); at smaller Fourier numbers the higher terms matter more.

(c) One-term inversion: 0.85 = (4/π)e^{−(π²/4)Fo} gives Fo = 0.1638 and t = Fo L²/α = **999 s**; `brentq` on the full series gives **962 s**. The 55 mm tile buys 16 minutes, not 20.

**Check.** 50 mm at 1200 s gave 594.6 K (Fo = 0.238); a thicker tile at the same time has a smaller Fo and must be cooler: 521.8 < 594.6 ✓.

**Common mistakes.** Using L in millimetres in L²/α; dropping the sign alternation (−1)ⁿ at the back face; reading the 4/π factor as 4/(nπ) with n = 1 instead of (2n + 1) = 1.

## 2. The mixed-boundary slab and the harmonic term

(a) Element between x and x + dx, unit area. Stored energy ρc_pT dx; rate of change ρc_pT_t dx. Net inflow q(x) − q(x + dx) = −q_x dx. **Fourier's law** (assumption: heat flows down the gradient with constant k) gives q = −kT_x, so −q_x = kT_xx. Equate and **divide by dx** (assumption: the element is small enough that T_t and T_xx are uniform across it, i.e. the continuum limit): ρc_pT_t = kT_xx, T_t = αT_xx with α = k/(ρc_p).

(b) θ = X(x)Θ(t): Θ′/(αΘ) = X″/X = −λ². X = A sin λx + B cos λx; X(0) = 0 ⇒ B = 0; X′(L) = Aλ cos λL = 0 ⇒ cos λL = 0 ⇒ λL = (2n + 1)π/2. With λ_n = nπ/L instead, X′(L) = λ cos nπ = ±λ ≠ 0: those functions have the wrong slope at the insulated face and *do* satisfy X(L) = 0, which is the other problem (a face held at T_s), not this one.

(c) b_n = (2/L)∫₀ᴸ θ₀ sin λ_n x dx = (2θ₀/L)(1 − cos λ_n L)/λ_n = (2θ₀/L)(1/λ_n) = (2θ₀/L)(2L/((2n + 1)π)) = 4θ₀/((2n + 1)π). (The eigenfunctions are orthogonal on [0, L] with norm L/2, as for the full sines, because λ_n L is an odd multiple of π/2.)

(d) u = sin kx sinh ky: u_xx = −k² u, u_yy = +k² u (sinh″ = sinh), sum zero ✓. u = sin kx sin ky: u_xx = −k²u, u_yy = −k²u, sum −2k²u ≠ 0. Laplace needs X″/X = −Y″/Y: if one is −k² the other is +k², so one direction oscillates and the other grows or decays (sinh). The heat equation's Θ′ = αΘ·(X″/X) must decay for a physical solution, so X″/X < 0 and X oscillates; the time direction is then exponential decay, not growth.

**Common mistakes.** Writing the Neumann condition as T = 0; forgetting that the orthogonality norm is still L/2 for the quarter-wave set.

## 3. Acoustic modes of a combustion chamber ($s = 15$: L = 0.95 m, c = 1100 m/s)

(a) p′ = X(x)Θ(t), X″ + k²X = 0, X = A sin kx + B cos kx. X′(0) = Ak = 0 ⇒ A = 0; X′(L) = −Bk sin kL = 0 ⇒ kL = nπ, n = 0, 1, 2, …. Modes cos(nπx/L), f_n = nc/(2L). The n = 0 mode is a spatially uniform pressure with Θ″ = 0: a constant (or linear drift), the "bulk" mode; in a real chamber it is the Helmholtz/chugging mode set by the nozzle, not by (5).

(b) f₁ = 1100/(2 × 0.95) = **579 Hz**, f₂ = **1158 Hz**, f₃ = **1737 Hz**.

(c) Closed–open: X(0) = 0 at the open end, X′(L) = 0: sin kx with cos kL = 0, kL = (2n − 1)π/2, f_n = (2n − 1)c/(4L): **289, 868, 1447 Hz**. The fundamental halves (a quarter wave fits instead of a half wave) and only odd harmonics appear.

(d) `string_mode_frequency(n, 0.95, 1210.0, 0.001)` with T/μ = 1.21 × 10⁶ = 1100² gives 579, 1158, 1737 Hz: the function depends on T and μ only through c = √(T/μ), and the fixed–fixed string and the closed–closed tube share the same eigenvalue set nπ/L (sine vs cosine does not change the frequencies). The closed–open tube's eigenvalues (2n − 1)π/(2L) are those of the lesson's *tile* (3), not of the string; `string_mode_frequency` cannot produce them without a fictitious doubling of L and discarding the even harmonics.

**Check.** Dimensions: c/L = s⁻¹ ✓. Limiting case: L → ∞ makes all f_n → 0 ✓.

**Common mistakes.** Forgetting n = 0 for Neumann–Neumann; mixing up which end has the node.

## 4. FTCS at the stability edge ($s = 15$: T_s = 1250 K, 51 points, 600 s)

| r | steps | max error [K] | back face [K] |
|---|---|---|---|
| 0.25 | 1191 | 0.122 | 376.68 |
| 0.50 | 596 | 0.365 | 376.44 |

(a) r = 0.25 is three times more accurate. At r = 0.5 the amplification factor of the 2Δx zigzag is g = 1 − 4(0.5) = −1: that mode is neither damped nor grown, so any zigzag introduced by the initial step at the hot face persists forever, and the truncation error of the scheme (which is O(Δt) + O(Δx²) = O(Δx²) at fixed r, with a coefficient that vanishes at r = 1/6 and is largest at r = ½) is at its worst.

(b) Unguarded march at r = 0.51 (the back face by the ghost-point rule): max |T| = 1250 K at step 100 (still the boundary value: the growing mode starts from rounding), **8.1 × 10⁷ K at step 500, 1.9 × 10¹⁶ K at step 1000**. Growth over 500 steps: 2.4 × 10⁸; |1 − 4r|⁵⁰⁰ = 1.04⁵⁰⁰ = 3.3 × 10⁸. The factor per step is |g| = 1.04, as von Neumann predicts.

(c) At κΔx = π the Fourier mode is e^{iκjΔx} = e^{iπj} = (−1)^j: it alternates sign at every grid point. Its amplification factor is g = 1 − 4r sin²(π/2) = 1 − 4r, the most negative of all modes; it is the first to become unstable, and the instability looks like a sawtooth that grows by |1 − 4r| per step.

(d) At r = 0.25, Δt = rΔx²/α and 600 s takes 1191 steps with 51 points; with 101 points Δx halves, Δt quarters, and **4762 steps** are needed. Explicit schemes pay Δx⁻² in time steps for every refinement in space (Δx⁻³ work in all, in 1D), which is why fine-grid heat problems use implicit schemes (Phase 03).

**Common mistakes.** Taking r = 0.5 as "the best because the largest stable step"; counting steps with the unrounded Δt.

## 5. Where you pluck matters ($s = 15$: a = 0.45 m = 3L/8)

(a) b_n/b₁ = 1, 0.1913, −0.0460, −0.0676, −0.0166, 0.0213, 0.0204, **0** for n = 1 … 8. Mode 8 is exactly missing: sin(8π × 3/8) = sin 3π = 0. The pluck point is a node of mode 8 (and of 16, 24, …): a shape that is zero at the pluck point receives no coefficient from a triangle peaked there. Modes 3 and 5 are small because sin(3π × 3/8) and sin(5π × 3/8) are near ±0.38 and ±0.92 divided by n²; the 1/n² decay dominates.

(b) The two solutions differ by at most 0.035 mm at t = 0 (at x = 0.45 m, the kink), 0.017 mm at t = L/4c = 2.12 ms and at t = L/2c = 4.24 ms (at x = 0.15 m, where a kink of the travelling wave sits). The difference is the Gibbs error of the 60-term series at the kinks of the exact (d'Alembert) solution; away from the kinks the difference is an order of magnitude smaller.

(c) 241 points, s = 0.8, 354 steps to 10 ms: E₀ = 8.849 mJ, E = 8.857 mJ, drift **+0.084 %**; maximum error against d'Alembert **0.024 mm** (0.5 % of h). The exact triangle energy ½Th²(1/a + 1/(L − a)) = 8.889 mJ; the 0.4 % deficit of E₀ is the central difference across the kink.

(d) 121 points, 177 steps: drift **+0.17 %**, error 0.038 mm. The drift doubles when Δx doubles: the energy error is made in the one or two cells around each kink where the central difference mis-estimates u_x, so it is proportional to Δx (first order), not Δx² like the smooth-solution error.

**Check.** Energy is conserved exactly by the continuous equation; a drift that halves with the grid is a discretisation effect, not a bug.

**Common mistakes.** Normalising by b₁ with the sign of b₁ ignored; reporting energy with u_t from a one-sided difference (that gives a 5 % error that does not converge).

## 6. Trusting the tile model

(a) Ranking by effect on the 1200 s back-face temperature (50 mm tile, nominal 595 K): **constant cold-value α** is largest: doubling α (a crude stand-in for k tripling at temperature) raises the back face to 907 K, so the model is *low* by up to 300 K; **fixed hot face at 1300 K from t = 0** is second: a real surface heats up over the first minute or two and peaks later, so the model is *high* by tens of kelvin; **insulated back** is third: the aluminium skin and the bond line absorb heat, so the model is *high* by a smaller margin that depends on the skin thickness. The one that pushes the answer the wrong way (unsafe) is the α assumption.

(b) With 2α: back face **906.8 K** instead of 594.6 K, required thickness **86.9 mm** instead of 61.4 mm (×√2, as the Fourier scaling demands). The crude bound says the property uncertainty is worth 25 mm of tile; the real k(T) curve, not a constant, is what the design needs.

(c) Instrument one tile with thermocouples at x = L/2 and x = L, hold its face at a fixed temperature in a furnace, and record both histories. Equation (4) says θ(L/2, t)/θ(L, t) is a function of Fo alone (0.86, 0.75, 0.71 at Fo = 0.06, 0.12, 0.24): fit α to the back-face history, then the mid-depth history is a prediction with no free parameter. A systematic miss, growing as the tile heats, means α is not constant.

(d) "The series is the exact solution *of the model*: constant α, step-heated face, adiabatic back. Each of those is wrong by a known amount and the first alone moves the answer by 25 mm. The 61.4 mm is a model output to be reported as 60 mm ± the α uncertainty, with the k(T) test of (c) as the next step, not a design value."

## Check-yourself answers

1. T_t = αT_xx; α = k/(ρc_p) in m²/s. The Fourier number Fo = αt/L².
2. sin(nπx/L) with decay rate α(nπ/L)²; mode 3 decays 9 times faster.
3. λ_n = (2n + 1)π/(2L): the fixed face needs a node (sine) and the insulated face an antinode (zero slope), so an odd number of quarter-waves fits.
4. c = √(T/μ); f_n = nc/(2L). Plucking at L/4 puts the pluck at a node of mode 4: sin(4π/4) = 0, so b₄ = 0.
5. u = ½[F(x − ct) + F(x + ct)] with F the odd 2L-periodic extension of the initial shape; the odd extension makes the ends stay at zero, and periodicity makes reflected waves return inverted.
6. For heat, X″/X = Θ′/(αΘ) must be negative so that time decays, giving oscillatory X. For Laplace, X″/X = −Y″/Y, so one direction oscillates and the other is exponential (sinh): the solution is fixed by data on all four edges, not by initial data.
7. FTCS: r = αΔt/Δx² ≤ ½ (no mode may be amplified; the 2Δx zigzag is the first to go). Leapfrog: s = cΔt/Δx ≤ 1 (a wave may not cross more than one cell per step: the numerical domain of dependence must contain the physical one).
8. Fo = αt/L² = 4.96 × 10⁻⁷ × 4800/0.1² = 0.238, the same Fourier number as the 50 mm tile at 1200 s (four times the time, four times L²), so the back face is the same: 595 K.
