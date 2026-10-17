# P05.18  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
(a) I = bh³/12 in m⁴; EI with E = 69 GPa. (b) P = mg; `deflection_point_load` at x = a, 0.150 and 0.300 m (the last is on the straight tail). (c) σ = 6Pa/(bh²); P_max = F_allow bh²/(6a) with F_allow = 72.5 MPa; divide by g for grams. (d) Relative terms u_E/E, u_b/b, 3u_h/h, 3u_a/a in quadrature; multiply by δ. Check with `predict_deflection(P, a, E, u_E, b, u_b, h, u_h, u_a)`. (e) q = ρgbh; `deflection_self_weight(a, q, L, EI)`; cancel by differencing.

## 2.
(a) Integrate twice; apply w(0) = 0 and w′(0) = 0; evaluate w and w′ at a; for x > a, w″ = 0 gives a line through (a, w(a)) with slope w′(a). (b) ln δ = ln 4 + ln P + 3 ln a − ln E − ln b − 3 ln h; first-order propagation of independent terms. For a fixed mark: ln g = 2 ln x + ln(3a − x) − ln 6; differentiate with respect to ln a. (c) w(a) = Pa³/(3EI) + θ₀a with θ₀ = Pa/k_θ; divide by P; factor. E_apparent = a³/(3 s I). (d) Write the loaded reading as r₀ + w_q + w_P and the unloaded as r₀ + w_q; subtract. The assumption is linearity (superposition).

## 3.
(a) I, EI for the rule; P = mg; linear Pa³/(3EI); `elastica_deflection(P, a, EI)`; the ratio; σ = 6Pa/(bh²). (b) `max_safe_load(107.5e6, a, b, h)`; repeat. `brentq` on linear/elastica − 1.01 over a bracket of masses. (c) Quadrature of 5 %, 3u_h/h, u_b/b, 3u_a/a; compare with the nonlinearity. The scale-factor uncertainties (E, h, b) make every point too high or too low by the same factor; the nonlinearity bends the curve. Compare the departure at the top load with the noise on one difference. (d) Look at the safety box: what can a hardened spring-steel rule do if it snaps or slips?

## 4.
(a) `data = make_synthetic_cantilever_data(s)`; P = mass × g; d = loaded − unloaded; `fit_line(P, d, sqrt(2)*0.30e-3)`; `modulus_from_slope(...)`; `en_number(E, uE, 69e9, 2e9)`. (b) Profile rows: x, w = loaded − unloaded; `fit_profile(x, w, 0.7*g, a, u_d)`. (c) From `data["truth"]`: EI_true = E b h³/12; s_true = a³/(3EI_true) + a²/k_θ; E_apparent = a³/(3 s_true I_nominal). Compare with (a) by the difference over its uncertainty. (d) w = `deflection_point_load(x, 0.7*g, a, EI)` with the nominal EI; `fit_profile(x, w, 0.7*g, a, sqrt(2)*0.01e-3)`; compare u(θ₀) with θ₀ = 0.7 g a/k_θ.

## 5.
(a) Four terms and their quadrature sum. (b) Replace the slope term; the thickness term is now the largest; decide what u_h halves it, then what u_a brings the total under 1 %. (c) Compare the 2 % bias with the target. Remove it with a profile fit at one load using the dial gauge; the correction carries an uncertainty u(θ₀) a/δ; add it in quadrature. (d) At a fixed mark, the load position enters with sensitivity 3a/(3a − x) and the mark position with d ln g/d ln x; both have uncertainties.

## 6.
For each candidate write E_apparent/E_true: a error → (a_true/a_used)³; h error → (h_true/h_used)³; clamp → 1/(1 + 3EI/(k_θ a)). Set each equal to 56/69 and solve. Judge plausibility from the method described (where was a measured from? was h measured or taken from the label? how was the strip clamped?). Cheapest tests first: re-measure a from the clamp edge, measure h with a caliper, then a profile fit for θ₀. Then judge the conclusion against E_n.
