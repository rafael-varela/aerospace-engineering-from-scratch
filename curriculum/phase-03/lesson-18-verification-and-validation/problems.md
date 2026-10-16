# P03.18  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer, give every error with the resolution it belongs to, and label every data set you use as measured or synthetic.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: grid convergence by hand)

A drag coefficient is computed on three grids refined by r = 2: C_D = 0.0200 + 1 × 10⁻⁵ s on the fine grid, 0.0200 + 4 × 10⁻⁵ s on the medium grid and 0.0200 + 17 × 10⁻⁵ s on the coarse grid.

(a) Compute the observed order p, the Richardson-extrapolated C_D and the GCI of the fine grid (F_s = 1.25) by hand.
(b) Recompute the extrapolated value and the GCI with the formal order 2 and compare; compute the asymptotic-range ratio with p = 2 and say whether the grids are in the asymptotic range.
(c) By what factor must the fine grid be refined for a GCI of 0.1 %?

Answer form: p, two extrapolated values, two GCIs in %, one ratio, one factor.

## 2. (Derivation: truncation error, exact discrete solution, Richardson, planted bug)

(a) Derive the truncation error (2) of the θ-method with central differences by Taylor expansion about (x_i, t_{n+½}). Show that the O(Δt) term vanishes only for θ = ½.
(b) Show that u = 2t x(L − x) with q = 2x(L − x) + 4αt satisfies the Crank–Nicolson discrete equations (1) exactly for any Δx and Δt with zero Dirichlet ends.
(c) From f_h = f + Chᵖ on grids h, rh, r²h derive p and the Richardson extrapolation (5). Show that the asymptotic ratio GCI₃₂/(rᵖGCI₂₁) equals f₁/f₂ when p is the observed order.
(d) Show that evaluating the source at t_n instead of the θ-weighted level adds (½ − 0)Δt q_t + O(Δt²) to the truncation error for θ = ½, and hence that the observed order with Δt ∝ Δx becomes 1.

Answer form: four derivations; (b) ends with an identity, (c) with formulas (5).

## 3. (Unfamiliar scenario: verifying an ODE integrator by manufactured solution)

A probe cools by radiation in vacuum: dT/dt = −k(T⁴ − T_a⁴), k = 2 × 10⁻¹¹ K⁻³ s⁻¹, T_a = 250 K, T(0) = (350 + s) K. You have written Heun's method (RK2: k₁ = f(T, t), k₂ = f(T + hk₁, t + h), T ← T + ½h(k₁ + k₂)) and want to verify it without an exact solution of the nonlinear equation.

(a) Manufacture the solution T_m(t) = T_a + (T(0) − T_a)e^{−t/τ}, τ = 500 s: write the source q_m(t) that makes dT/dt = −k(T⁴ − T_a⁴) + q_m(t) have T_m as its solution, and evaluate q_m(0) by hand.
(b) Implement the modified ODE and Heun's method, integrate to t = 1000 s with n = 10, 20, 40, 80, 160 steps, tabulate |T_n − T_m| and the observed orders.
(c) Plant a bug: evaluate k₂ at time t instead of t + h. Repeat (b). Explain the result with the truncation error of the stage.
(d) Why does the test in (b) also verify the integrator for the *unmodified* cooling equation?

Answer form: one formula and one number, two tables, two short paragraphs.

## 4. (Numerical, uses your code: your own manufactured solution)

Build the manufactured solution u_m = (20 + s) sin(πx/L)(1 − e^{−t/τ}) + s·x(L − x)/L²·sin(2πt/τ), τ = 200 s, with a loss coefficient β = s × 10⁻⁴ s⁻¹ and the lesson's α and L.

(a) Let `manufactured_source` build q_m; check one value of q_m by hand (at x = 0.1 m, t = 50 s).
(b) Run the Crank–Nicolson convergence study on 25–400 cells with Δt = 300 s/nx; report the errors and observed orders. Repeat with `bug="source_time"`.
(c) Run the heater problem with your β on 25, 50, 100 cells, report the temperature at the heater at 600 s, the observed order, the extrapolated value and the GCI, check against an 800-cell run, and compute the refinement needed for a GCI of 0.01 %.

Answer form: one hand value, two convergence tables, one Richardson table with a check.

## 5. (Numerical, uses your code: validation with your own data set)

Generate synthetic data with `synthetic_measurements(HeatSetup(), seed=s)` (label it synthetic in everything you write).

(a) Evaluate the lossless model on 200 cells with Δt = 1 s at the sensors and report the validation metrics (points within 2u_val, reduced χ², largest |E|/u_val, E at t = 600 s). Take u_num from the lesson's GCI.
(b) Fit β on the first five sampling times only (hold out t = 600 s). Report β̂ and the hold-out row of E/u_val. Then fit on all six and report β̂ with its standard uncertainty from the curvature of χ² (σ_β = √(2/χ²'')). How many standard uncertainties is the truth away?
(c) Convert β̂ to a convective coefficient h for a 10 mm square aluminium section (ρ = 2700 kg/m³, c = 900 J/kg K). Then predict the heater temperature at 1200 s with β̂, with β_true and with β = 0. Which of the three numbers would you put in a design report, with what uncertainty, and what sentence must accompany it?

Answer form: one metrics table, two fits with uncertainties, three predictions and two sentences.

## 6. (Judgment: a vendor says "validated")

A CFD vendor states that its code is "validated for internal flows" and shows a plot in which its pressure-drop prediction for a duct lies within 3 % of a textbook correlation. Your project needs the heat-transfer coefficient in a cooling channel with a 90° bend to ±10 % for a thermal margin that decides whether a battery pack needs a second fan.

(a) Using the NASA-STD-7009 vocabulary, say what the vendor's plot does and does not demonstrate, and list the three pieces of verification evidence you would ask for first.
(b) Write three verifiable requirements (Appendix C style: one statement, a measurable criterion, a method) for the simulation you would run.
(c) The vendor's grid study on your geometry gives heat-transfer coefficients of 412, 395 and 390 W/m²K on grids refined by r = 1.5. Compute the observed order and the GCI; is the fine grid adequate for ±10 %? For ±2 %?
(d) The computed coefficient agrees with a single test point to 4 %, but the test's own uncertainty is 8 %. Write the two-sentence statement for the design review, including what the agreement proves and what it does not.

Answer form: four short paragraphs, each with at least one number.
