# P08.14  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer. Label every velocity profile you use as **synthetic** or as **published data** (with its source and terms of use), and every skin-friction formula as a **correlation**.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

The light aeroplane at 3000 m, 60 m/s: ρ = 0.9091 kg/m³, ν = 1.8630 × 10⁻⁵ m²/s, q = 1636.4 Pa.

## 1. (Warm-up: the fuselage)

Treat the light aeroplane's fuselage as a flat plate of length L_f = 5.5 + s/15 m, turbulent from the nose.

(a) Compute Re_L and C_F from the four integrated correlations in `CF_plate`. How far apart are they?
(b) At the tail, compute the local c_f (White), the friction velocity u_τ, the wall shear stress, the viscous length ν/u_τ, δ from the 1/7th law and δ⁺ = δu_τ/ν.
(c) A CFD engineer wants the first cell at y⁺ ≈ 1 at the tail. How high is it? Where does the log layer (y⁺ = 30) start?
(d) By what factor does the turbulent C_F exceed the all-laminar Blasius value at the same Re_L?

Answer form: a Reynolds number, four coefficients and a spread; six numbers with units; two lengths; a factor.

## 2. (Derivation: the log law, the 1/7th law, and the mixed plate)

(a) In the overlap between the inner layer (where U depends on y, u_τ and ν) and the outer layer (where the velocity defect depends on y/δ and u_τ), both descriptions must hold. Show that this forces dU/dy = u_τ/(κy) and hence the log law u⁺ = ln(y⁺)/κ + B. Show that Spalding's formula reduces to u⁺ = y⁺ for small u⁺ and to the log law for large u⁺.
(b) With u/U = (y/δ)^(1/7), show that δ*/δ = 1/8, θ/δ = 7/72 and H = 9/7. With Blasius's pipe-flow wall shear τ_w = 0.0225ρU²(ν/(Uδ))^(1/4) and dθ/dx = c_f/2, derive δ/x = 0.3707Re_x^(−1/5) and c_f = 0.0577Re_x^(−1/5). Show that the average C_F over a plate is 5/4 of c_f at its end.
(c) Derive Prandtl's mixed-plate formula C_F = C_F,turb(Re_L) − (Re_tr/Re_L)[C_F,turb(Re_tr) − C_F,lam(Re_tr)] and state the assumption about the turbulent layer that it makes.

Answer form: three derivations; boxed results.

## 3. (Unfamiliar scenario: a rake on the wing)

A boundary-layer rake measures 25 points log-spaced from y = 0.05 mm to δ at x = 1.2 m on the light aeroplane's wing. Since no flight data are available, make the profile **synthetic**: `synthetic_profile(y, u_tau, nu, delta, noise=0.02, seed=s)`, with u_τ from White's local correlation at x = 1.2 m and δ from `power_law_bl`.

(a) Run `clauser_fit` with the default window (y⁺ ≥ 30, y/δ ≤ 0.2). How many points does it use, and how close is u_τ to the value that built the profile? Convert to τ_w.
(b) Repeat with y/δ ≤ 0.1 and ≤ 0.3, with y⁺ ≥ 50, and with κ = 0.384, B = 4.17. Tabulate u_τ and its error.
(c) Which choice moves the answer most, and why does the choice of log-law constants matter so little here? Would it matter as little for a profile measured only at y⁺ > 500?

Answer form: a count, two velocities, a percentage and a stress; a table of five rows; two sentences and a yes/no with a reason.

## 4. (Numerical, uses your code: correlations and the mixed plate)

Let Re_L = 3 × 10⁶ + s × 10⁵.

(a) Compare C_F from the integrated correlations with the integral of each local correlation (`integrate_cf`). What is the spread, and is it larger than the difference between "local integrated" and "integrated" for the same author?
(b) For Re_tr = 3 × 10⁵, 5 × 10⁵, 10⁶ and 2 × 10⁶, compute the mixed-plate C_F by Prandtl's subtraction (with the Prandtl–Schlichting and with the power-law correlation) and by the momentum-thickness match. Report the virtual origin x₀/L.
(c) Explain why the two methods drift apart as Re_tr grows. Which would you trust for a laminar-flow wing, and why?

Answer form: seven coefficients and a sentence; a table and four origins; a paragraph.

## 5. (Numerical, uses your code: a friction law with no correlation in it)

(a) Compute c_f from `log_wake_plate` at Re_x = 10⁶, 4.83 × 10⁶ and 10⁷ for Coles's wake parameter Π = 0.4, 0.55, 0.7, and with κ = 0.384, B = 4.17 (Π = 0.55). Compare each with White's correlation.
(b) Tabulate `log_wake_state` at δ⁺ = 10³, 10⁴, 10⁵: Ue⁺, c_f, Re_θ and H. How does H change with Reynolds number, and what does that say about the turbulent profile's fullness?
(c) The log-wake law and the correlations differ by a few per cent. Which inputs to the log-wake law would you want to measure to settle the difference, and which published data could you compare with (without redistributing them)?

Answer form: a table of 12 numbers; a table and two sentences; a short paragraph.

## 6. (Judgment: what do you tell the designer?)

The lesson's estimates of the wing's friction drag at cruise are 29.7 N (all laminar), 83–95 N (transition at the e⁹ point of `P08.12`), 150–153 N (Re_tr = 5 × 10⁵) and 163–168 N (turbulent from the leading edge, four correlations).

(a) Which number would you give for a drag estimate of the production aeroplane, and with what range? Justify each end of the range.
(b) Name two effects that none of the flat-plate numbers include and say which way each moves the drag.
(c) What single flight or tunnel measurement would narrow the range the most?

Answer form: a number with a range and two justifications; two effects with directions; one measurement.
