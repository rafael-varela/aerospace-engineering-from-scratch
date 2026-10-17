# P03.16  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for $s = 15$; numerical results are from `solution.py` and given to the figures shown.

## 1. One Gauss–Seidel sweep by hand ($s = 15$: b = 11.5 m, S = 16.8 m²)

**Derivation.** With W_wing = 0, L = 25 000 N. Parasite drag q(f_rest + c_dw S) = 2908.8 × (0.25 + 0.1344) = 1118 N; induced drag L²/(qπeb²) = 6.25 × 10⁸/(2908.8 × 2.670 × 132.25) = 608.5 N; **D = 1726.6 N**, L/D = 14.48. φ = D/L = 0.06906, κφ = 0.1778, **W_fuel** = 25 000 × (e^{0.1778} − 1) = **4864 N**, **W_TO = 29 864 N**, **M_root** = 3.75 × 29 864 × 11.5/(3π) = **136.65 kN m**. Discipline 2: h = 0.12 × 16.8/11.5 = 0.1753 m, A = 136 650/(0.1753 × 120 × 10⁶) = 64.96 cm², W_caps = 2700 × 9.807 × (2/3) × 6.496 × 10⁻³ × 11.5 = 1318.6 N, **W_wing = 1318.6 + 60 × 16.8 = 2326.6 N**. `aero_mission` and `structures` give the same to four figures.

(b) At the converged point (W_wing = 2449.9 N, M_root = 149.42 kN m) the partials give **k = 0.0506**. The error after the first sweep is |2326.6 − 2449.9| = 123.2 N; the target 10⁻⁶ × 2450 = 2.4 × 10⁻³ N needs n = ln(2.4 × 10⁻³/123.2)/ln 0.0506 = **3.6, so 4 more sweeps**. `mda_gauss_seidel` reaches a 10⁻¹² relative residual in 10 sweeps, with residuals 1.21 × 10⁴, 614, 31.1, 1.57, … each 0.0506 times the previous.

(c) No: with k = 0.05 each sweep gains 1.3 decades, four sweeps settle the design, and Newton would need the coupling derivatives for no visible gain. Newton pays off when k approaches 1 (part 4 of the lesson) or when the disciplines are expensive and every evaluation counts.

**Check.** Dimensions: M_root = N × m ✓; A = N m/(m × Pa) = m² ✓. The first-sweep W_wing (2327 N) is below the converged 2450 N because it was sized for a wing of zero weight.

**Common mistakes.** Using the 1 g moment instead of the ultimate (n = 3.75); forgetting the skin weight term; taking k from the first-sweep point rather than the converged one (the difference is small here but not in general).

## 2. Coupling, adjoint, deflection

(a) W_TO = L e^{κφ} with φ = f₀/L + L/(qπeb²), f₀ = q(f_rest + c_dw S). ∂φ/∂L = −f₀/L² + 1/(qπeb²), so ∂W_TO/∂L = e^{κφ}(1 + Lκ ∂φ/∂L) and **∂M_root/∂W_wing = (n b/(3π)) e^{κφ}(1 + κ(−f₀/L + L/(qπeb²)))**. From (4), **∂W_wing/∂M_root = (2/3)ρ_m g b²/(τSσ_a)**, a pure geometric factor. Their product k carries W₀ only through L inside e^{κφ} and the bracket: at the baseline 5.418 × 0.01051 = 0.0569; doubling W₀ changes it to 0.0642.

(b) Let y₂* be the fixed point. One sweep maps y₂ = y₂* + e to G₂(G₁(y₂* + e)) = G₂(G₁(y₂*) + G₁'e + O(e²)) = y₂* + G₂'G₁' e + O(e²) = y₂* + k e + O(e²). |e| shrinks every sweep iff |k| < 1: the contraction condition.

(c) dF/dx = ∂F/∂x + (∂F/∂y)(dy/dx) and dy/dx = −(∂R/∂y)⁻¹(∂R/∂x), so dF/dx = ∂F/∂x − (∂F/∂y)(∂R/∂y)⁻¹(∂R/∂x). Group the first two factors: ψᵀ := (∂F/∂y)(∂R/∂y)⁻¹, equivalently **(∂R/∂y)ᵀψ = (∂F/∂y)ᵀ**, and dF/dx = ∂F/∂x − ψᵀ ∂R/∂x. Direct: one solve with (∂R/∂y) per design variable, n_x solves. Adjoint: one solve with (∂R/∂y)ᵀ per function, n_F solves. Break-even at n_x = n_F; for the lesson's 2 and 2 it is a tie, for 200 variables and 3 functions the adjoint is 67 times cheaper.

(d) With the outer fibre at σ_a everywhere, M(y)/(EI(y)) = σ_a/(E h/2) = 2σ_a/(Eh), a constant κ_b. Integrating w'' = κ_b from the clamped root (w = w' = 0) gives w(y) = ½κ_b y², so at ultimate load w(b/2) = σ_a (b/2)²/(Eh); at limit load divide by 1.5: δ = σ_a (b/2)²/(1.5 E h). Divide by b/2 and substitute h = τS/b: **δ/(b/2) = σ_a b²/(3Eτ S)**. The limit δ_max gives AR ≤ 3Eτδ_max/σ_a: aluminium **12.6**; carbon (E = 140 GPa, σ_a = 400 MPa) 3 × 140 × 10⁹ × 0.12 × 0.06/(400 × 10⁶) = **7.56**. A stronger material used to its strength deflects *more*: stiffness-limited wings are designed below their strength, which is why the lesson's "fully stressed" spar is the assumption to doubt first.

**Common mistakes.** In (a), forgetting that W_fuel also depends on L (the e^{κφ} factor); in (c), transposing the wrong matrix (ψ solves the transposed system; the result formula uses the untransposed ∂R/∂x); in (d), using the ultimate-load deflection as the limit-load one.

## 3. An electric aircraft ($s = 15$: R = 350 km)

**Assumptions.** Constant L/D and η over the cruise, battery energy fully usable, no reserve, W_b proportional to total weight.

(a) W_b = gE/e = gRW/(η(L/D)e) = kW with **k = gR/(η(L/D)e)** = 9.807 × 3.5 × 10⁵/(0.75 × 15 × 9 × 10⁵) = **0.3390** (e = 250 Wh/kg = 9.0 × 10⁵ J/kg). W = W₀ + kW gives **W = W₀/(1 − k) = 12 103 N**, **W_b = 4103 N** (418 kg of battery for an 816 kg airframe). k = 1 means the battery needed to carry the aircraft's weight over the range equals the aircraft's weight: the fixed point runs away and no design closes.

(b) From W = W₀ the first error is kW/(1 − k)·(1 − k)… more simply e₀ = W − W₀ = kW; e_n = kⁿ e₀; relative error kⁿ·k/(1 − k) ≤ 10⁻⁶ gives n ≥ ln(10⁻⁶ (1 − k)/k)/ln k = **12.2, so 13 sweeps**. Newton: the residual R(W) = W − W₀ − kW is linear, so **one step** lands exactly on the answer.

(c) Total: dW/dR = W₀ (dk/dR)/(1 − k)² = W k/(R(1 − k)) = 12 103 × 0.339/(350 × 0.661) = **17.73 N/km**. Partial with W frozen: d(kW)/dR = kW/R = **11.72 N/km**. Ratio 1.513 = **1/(1 − k)**: every newton of battery added for range needs k newtons more battery to carry it, and so on, a geometric series summing to 1/(1 − k).

(d) k = 1 at **R = η(L/D)e/g = 1032 km**. The frozen-weight analysis would report a finite battery kW₀ = 8 kN at that range and a flyable aircraft of 16 kN; the coupled analysis says the aircraft does not exist.

**Check.** At R → 0, k → 0 and W → W₀ ✓. Units: gR/(η(L/D)e) = (m/s²)(m)/(J/kg) = (m²/s²)/(m²/s²) ✓ dimensionless.

**Common mistakes.** Using Wh/kg without converting to J/kg (a factor 3600); counting the first evaluation as a sweep; concluding from the partial that "range costs 11.7 N/km".

## 4. Your own aircraft ($s = 15$: W₀ = 27.5 kN, R = 1750 km)

(a) At (12 m, 16.8 m²): Newton converges to **M_root = 175.63 kN m, W_wing = 2853.3 N** in 2 iterations (residuals 5 × 10⁴, 7, 10⁻⁷: order 2.0); Gauss–Seidel needs 10 sweeps to 10⁻¹¹ with successive ratios 0.05902, 0.05902, 0.05902 against **k = 0.05902**. The baseline burns W_fuel = 6430 N (656 kg) at L/D = 15.63, W_TO = 36 783 N.

(b) Adjoint: dW_fuel/db = **−417.79 N/m**, dW_fuel/dS = **+76.14 N/m²**; finite differences (h = 10⁻⁶) agree to 5 × 10⁻¹¹ and 1 × 10⁻⁹ relative. The frozen partial ∂W_fuel/∂b = −499.64 N/m overstates the span benefit by **19.6 %**.

(c) MDF (adjoint, 11 SLSQP iterations, 45 MDA solves) and IDF (12 iterations) both give **b = 15.629 m, S = 19.386 m², AR = 12.60, W_fuel = 5631.3 N = 574.2 kg**, a saving of **81.4 kg (12.4 %)**. Active: the deflection (AR) limit to 10⁻¹⁴; stall margin 12.3 %, gate margin 21.9 %. W_wing = 4807 N, L/D = 18.69.

(d) λ_defl = **1316 N** per unit of the normalised constraint; 10 % more deflection allowance would save 131.6 N = **13.4 kg** of fuel per flight.

**Check.** Heavier aircraft and longer range than the lesson's: more fuel (574 vs 460 kg), larger wing (19.4 vs 17.5 m²), same aspect ratio because the limit is geometric ✓.

**Common mistakes.** Reporting k from the first sweep; comparing the FD derivative with the frozen partial and concluding the adjoint is wrong.

## 5. The finite-difference floor (baseline, errors in N/m)

| h (relative) | Newton 10⁻¹³ | GS 10⁻¹² | GS 10⁻⁶ | GS 10⁻³ |
|---|---|---|---|---|
| 10⁻¹ | 6.78 | 6.78 | 6.78 | 6.82 |
| 10⁻² | 6.66 × 10⁻² | 6.66 × 10⁻² | 6.59 × 10⁻² | 0.11 |
| 10⁻³ | 6.66 × 10⁻⁴ | 6.66 × 10⁻⁴ | 9.1 × 10⁻⁴ | 0.05 |
| 10⁻⁴ | 6.66 × 10⁻⁶ | 6.66 × 10⁻⁶ | 2.5 × 10⁻⁴ | 0.04 |
| 10⁻⁵ | 5.0 × 10⁻⁸ | 7.7 × 10⁻⁸ | 2.4 × 10⁻⁴ | 0.04 |
| 10⁻⁶ | 2.9 × 10⁻⁸ | 2.9 × 10⁻⁸ | 2.4 × 10⁻⁴ | 0.04 |
| 10⁻⁷ | 3.1 × 10⁻⁷ | 3.1 × 10⁻⁷ | 2.4 × 10⁻⁴ | 0.04 |

(a) Slope **2.00** over h = 10⁻¹ … 10⁻⁴ for the tight runs (the error falls by 100 per decade), until the rounding floor near 3 × 10⁻⁸ N/m (10⁻¹⁰ of the derivative, from cancellation in the difference of two fuel weights of 5 × 10³ N known to 10⁻¹⁶). The loose runs follow the same line until they hit floors of 2.4 × 10⁻⁴ (tolerance 10⁻⁶) and 0.04 N/m (10⁻³).

(b) Write the Gauss–Seidel answer as y_GS(x) = y(x) + ε(x) where ε is the stopping error, a smooth function of x of size about tol × |y| × k (the same number of sweeps is taken at x + h and x − h, so ε is nearly the same on both sides). The central difference of F(x, y_GS(x)) is dF/dx + (∂F/∂y)(dε/dx) + O(h²): the *value* of ε cancels, its *derivative* does not, and dε/dx is of order tol × |dy/dx|, independent of h. Hence a floor proportional to the tolerance (2.4 × 10⁻⁴/0.04 ≈ 10⁻⁶/10⁻³ × 6 ✓), not to 1/h. A floor growing like 1/h appears only when the two sides take different numbers of sweeps, which the h = 10⁻³ row for GS 10⁻⁶ begins to show.

(c) At h ≈ 10⁻⁷ SLSQP would get the gradient to 3 × 10⁻⁷/291 = 10⁻⁹ relative with either tight MDA, to 10⁻⁶ relative with GS 10⁻⁶ (still fine), and to 1.4 × 10⁻⁴ relative with GS 10⁻³, also still usable for descent but not for a tight KKT test: with a 10⁻³ loop the optimiser would wander in the last digits, stop early on "positive directional derivative", or report convergence to a point whose constraints are off by 10⁻³. With h = 10⁻³ or coarser steps the truncation error (6.7 × 10⁻⁴ N/m) would dominate everything else.

**Common mistakes.** Expecting the loose floor to grow as 1/h (that is the formula for *random* noise, not for a smooth systematic error); fitting the slope through the floor.

## 6. Scaling up

(a) **MDA.** Start with block Gauss–Seidel, because the commercial codes give no Jacobian; measure k from the ratio of the first three sweep residuals before anything else. k = 0.3 means 12 sweeps for 10⁻⁶; k = 0.9 means about 130 sweeps, each a panel-code run plus a finite-element solve, which is unaffordable. At k = 0.9 add Aitken acceleration or a Broyden (quasi-Newton) update of ∂R/∂y built from the sweep history, which converts the linear rate to superlinear without any discipline derivatives, or reformulate the coupling variables (for example pass the lift distribution rather than the total lift) to shrink k.

(b) **Derivatives.** Finite differences: 2 × 150 = 300 MDAs per gradient, each of several sweeps: hundreds of panel and FE runs per optimiser iteration. Direct: 150 coupled linear solves with ∂R/∂y, but ∂R/∂y and ∂R/∂x need the partials of each discipline. Adjoint: 13 solves (12 constraints + objective), the cheapest by a factor of 12 over direct and vastly over FD, but it needs partials from codes that do not provide them. The team can actually use only FD (with the floor problem of problem 5 multiplied by 150 variables), unless they replace the panel code by one with an adjoint (several open-source VLM/panel codes have one) and get stiffness-matrix derivatives from the FE model, which many FE codes do export; then the adjoint is available and the MDO becomes feasible at 150 variables.

(c) **Architecture.** MDF keeps every iterate consistent, so a run can be stopped and its last design used, and it needs no extra constraints; its cost is one MDA per evaluation, which with k = 0.3 and two expensive codes is 12 code runs per evaluation. IDF removes the inner loop and lets aerodynamics and structures run in parallel, at the price of adding the coupling variables (a lift distribution of perhaps 40 numbers, a deflection field of 40) and 80 equality constraints to the optimiser, and of inconsistent intermediate designs. With no derivatives and expensive codes, the number of code runs decides: IDF, because it runs each code once per evaluation and in parallel.

(d) **Checks before believing 7 %.** (1) The wing-weight law at the gate-limited span: the result sits on the span limit, so the fuel saving comes from induced drag bought with structure, and the weight model (a statistical equation? a sized FE model with buckling?) must be checked against an independent sizing of the final wing: ±20 % on wing weight is ±2–3 % on fuel. (2) The drag polar at the optimum's lift coefficient and aspect ratio, since the panel code's Oswald factor and profile drag at the new C_L drive the gain: compare with CFD or wind-tunnel data for a similar wing; a 0.05 error in e is 1.5 % of fuel here. (3) The mission model's fuel-fraction sensitivity (reserves, climb, descent credit): a Breguet-only cruise overstates the fuel benefit of L/D by the non-cruise fraction, typically 10–20 % of block fuel. Evidence: a sized structural model, measured or computed polars, and a flight-profile model, in that order; the 7 % is credible only after the first two.

## Check-yourself answers

1. x = (b, S); y = (M_root, W_wing); G₁(x, y₂): aerodynamics + mission, takes b, S, W_wing and returns M_root; G₂(x, y₁): structures, takes b, S, M_root and returns W_wing; R = y − G(x, y) = 0 defines a consistent design.
2. k = (∂M_root/∂W_wing)(∂W_wing/∂M_root), the product of the two off-diagonal entries of ∂G/∂y at the solution; k > 1 means Gauss–Seidel diverges and physically that each newton of wing weight generates more than a newton of additional wing weight through its own loads: the design does not close.
3. dy/dx = −(∂R/∂y)⁻¹ ∂R/∂x. Freezing W_wing drops the term (∂F/∂y)(dy/dx); since span raises the wing weight and weight costs fuel, the frozen derivative overstates the fuel benefit of span (−347 against −291 N/m at the baseline, 19 %).
4. (∂R/∂y)ᵀψ_i = (∂F_i/∂y)ᵀ, then dF_i/dx = ∂F_i/∂x − ψ_iᵀ ∂R/∂x. Direct: 200 solves; adjoint: 3.
5. Because the converged loop is only accurate to its tolerance; the stopping error is a smooth function of x, so its derivative survives the central difference as a floor of about (tolerance × derivative), while rounding adds a floor that grows as ε_machine F/h.
6. MDF converges the MDA inside every function evaluation and optimises over x only; IDF adds y to the design variables and imposes R(x, y) = 0 as equality constraints, converging the disciplines only at the optimum. Both solve the same KKT system at the solution, so they must agree.
7. Because sizing every section to σ_a makes the curvature 2σ_a/(Eh) regardless of the moment that produced it; the limit then constrains geometry only, b²/S ≤ 3Eτδ_max/σ_a, an aspect-ratio limit.
8. At the aspect-ratio limit more area permits more span (b = √(12.6 S)), and the induced-drag benefit of that span outweighs the profile drag and skin weight of the extra area.
