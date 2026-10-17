# P08.14  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. L_f from s; Re_L; `CF_plate(Re_L, kind)` for the four kinds; max/min. 2. `cf_local(Re_L, "white")`; u_τ = U√(c_f/2); τ_w = ρu_τ²; ν/u_τ; δ from `power_law_bl`; δ⁺. 3. y = y⁺ν/u_τ for y⁺ = 1 and 30. 4. Divide by `laminar_CF(Re_L)`.

## 2.
(a) Inner: U = u_τ f(y⁺) ⇒ y dU/dy = u_τ y⁺f′(y⁺); outer: (U_e − U)/u_τ = g(y/δ) ⇒ y dU/dy = −u_τ (y/δ)g′; equal for all y⁺ and y/δ ⇒ constant 1/κ; integrate. Spalding: expand the exponential; for large u⁺ the exponential dominates. (b) Integrals; θ = 7δ/72; substitute τ_w; integrate δ^(1/4)dδ; c_f from τ_w with that δ; average x^(−1/5). (c) Write drags per unit span as ½ρU²L C_F and subtract.

## 3.
1. c_f = `cf_local(U x/ν, "white")`, u_τ, δ. 2. y = `np.geomspace(0.05e-3, δ, 25)`; `synthetic_profile`. 3. `clauser_fit` default and the four variants; errors against the true u_τ; τ_w = ρu_τ². 4. Solve log_law(y⁺) = log_law(y⁺, 0.384, 4.17).

## 4.
1. `CF_plate` for four kinds and `integrate_cf` for three; spread. 2. Loop over Re_tr: `mixed_plate_prandtl` (two kinds) and `mixed_plate_momentum`. 3. Compare θ at transition with the turbulent θ that Prandtl's subtraction implies there.

## 5.
1. Loop over Π and constants: `log_wake_plate`; ratio to `cf_local(..., "white")`. 2. `log_wake_state(δ⁺)`. 3. κ, B, Π from a DNS or experimental profile; the Lee & Moser channel and Schlatter & Örlü boundary-layer profiles (downloaded by the learner under their citation terms).

## 6.
(a) Base on the e⁹ estimate, bounded below by a longer laminar run and above by transition near the nose. (b) Pressure gradient and form drag, roughness and steps, wing thickness, interference. (c) Measure the transition line (infrared thermography) in flight.
