# P03.12  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact, which are numerical, and what the numerical error is.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: sizing an advection run)

A fuel line carries flow at a = 8 + s/5 m/s. The model domain is L = 3 m on n = 150 cells, and a Gaussian slug of standard deviation σ₀ = 5 cm must travel 2 m.

(a) Give Δx, the time step for ν = 0.8, the travel time and the number of steps.
(b) Give the upwind numerical diffusivity D at ν = 0.8 and the slug's standard deviation on arrival. Repeat for ν = 0.4.
(c) How many cells would keep the upwind spreading below 10 % of σ₀ at ν = 0.8?

Answer form: seven numbers with units, one sentence on what (c) implies for first-order schemes.

## 2. (Derivation: the integral law, the fluxes, and the modified equation)

(a) Integrate u_t + (au)_x = 0 over cell i and over one time step and derive equation (1) of the lesson. State precisely which quantity in (1) is approximated and which are exact.
(b) Derive the Lax–Wendroff update from the Taylor expansion in time (use u_tt = a²u_xx) and show it is (1) with F = ½a(ū_i + ū_{i+1}) − ½aν(ū_{i+1} − ū_i).
(c) Derive the modified equation (2) of the upwind scheme by Taylor expansion, including the diffusivity D = ½aΔx(1 − ν). Show that D = 0 at ν = 1 and explain, from the stencil, why the scheme is then exact.
(d) Derive |g|² for the upwind scheme and show that 0 ≤ ν ≤ 1 is necessary and sufficient for |g| ≤ 1 at every θ.

Answer form: four derivations with every step shown.

## 3. (Unfamiliar scenario: a line whose speed changes)

Downstream of the sensor the feed line narrows so that the bulk speed rises from a₁ = 10 m/s to a₂ = 20 m/s over a short converging section. Two models are proposed for a scalar carried by the flow: the advective form u_t + a(x)u_x = 0 and the conservation form u_t + (a(x)u)_x = 0.

(a) Show that for the conservation form the quantity q = a(x)u satisfies q_t + a q_x = 0, so that q is constant along characteristics. What happens to the peak and the length of a slug that passes from a₁ to a₂ under each model?
(b) The slug's mixture ratio is an intensive property (it does not change when the fluid accelerates); the mass of excess oxidiser per metre of line is extensive. Which model is right for which quantity, and what is the relation between them?
(c) Write the upwind flux for the conservation form on a grid where a is known at the interfaces, F_{i+1/2} = a_{i+1/2}ū_i, and show that the update conserves Σū_iΔx but not Σ(aū)_iΔx. Which of the two is physically conserved in the narrowing line?
(d) A colleague runs the advective form with the lesson's periodic code and a(x) varying. Explain in two sentences why that run cannot conserve anything and what boundary condition the problem actually needs.

Answer form: a derivation and two numbers for (a), a short paragraph each for (b), (c), (d).

## 4. (Numerical, uses your code: the sensor reading)

Change the square slug's width to 0.2 + s/100 m (all else as in `feed_line`), use ν = ½, and run `sensor_history` for the four schemes to t = 0.2 s.

(a) On n = 120 cells, tabulate the 10 % and 50 % arrival times, the peak and the minimum read by the sensor for each scheme, together with the exact front-arrival time at the sensor cell's centre.
(b) Repeat on n = 480 cells. Which readings converged, which did not, and which got worse?
(c) The controller trips when the reading exceeds 0.9. Which scheme–grid combinations trip, and which agree with the exact solution?
(d) Explain the early 10 % arrival of the first-order schemes from the modified equation, and the late 50 % crossing of Lax–Wendroff from its phase error.

Answer form: two 4-row tables, one list, one short paragraph.

## 5. (Numerical, uses your code: the Courant number as a knob)

For the Gaussian slug on n = 240, 480, 960, 1920 cells, run `convergence_study("upwind", …)` at ν = 0.25, 0.5 and 0.9.

(a) Tabulate the errors and observed orders. Predict from D = ½aΔx(1 − ν) the ratios e(0.25)/e(0.9) and e(0.5)/e(0.9) and compare with the finest grid.
(b) Run ν = 1 and report the error. Then, in a copy of `fv_step` without the CFL check, run ν = 1.05 for 300 steps from the Gaussian plus 10⁻⁶(−1)ⁱ and report max|u| after 100, 200 and 300 steps against |g(π)|ᵏ × 10⁻⁶.
(c) Why does the ν = 1.05 run look harmless for the first hundred steps even though it is unstable from the first?
(d) In production codes ν is usually set near 0.5–0.9 rather than 1. Give two reasons from this lesson and one from `P03.05`.

Answer form: a 3-row table with the predicted ratios, three numbers for (b), two short paragraphs.

## 6. (Judgment: choosing the scheme for the controller model)

The feed-line model will drive a controller simulation; the mixture fraction it reports must stay between 0 and 1, the front arrival must be within 2 ms, and the peak within 5 %. The real line has a turbulent eddy diffusivity of about 10⁻³ m²/s.

(a) At ν = ½, what grid would make the upwind scheme's numerical diffusivity smaller than the physical one, and how many steps would 0.1 s take? Repeat for ν = 0.9.
(b) Rank the four schemes against the three requirements on the 120-cell grid, using the lesson's sensor table. Which single requirement eliminates Lax–Wendroff regardless of grid?
(c) The minmod scheme's peak on the Gaussian slug is 0.875 on 120 cells (12 % low) while on the square slug it is 0.981. Explain the difference and say which one matters for the controller.
(d) Someone proposes to "fix" the Lax–Wendroff undershoot by clipping negative values to zero after each step. What does that do to conservation, and why is a limiter the honest version of the same idea?

Answer form: four short paragraphs, each with at least one number.
