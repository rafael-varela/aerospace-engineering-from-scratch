# P10.14  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. `reynolds_number(V, 0.2)` for three V. 2. `section_cd_friction(Re, 0.12, Re_tr=1e12)`. 3. V = 5 + s/3, c = 0.05 + s/200; Re; compare with `natural_transition_re("0012", 0.445)`.

## 2.
(a) Substitute Ue = Cs^m in (1); Ue′ = mCs^{m−1}; λ = θ²Ue′/ν. λ = 0.075 at m = 1. Separation needs λ ≤ −0.09: solve 0.45m/(5m + 1) = −0.09 (m < 0). (b) θ² ∝ ν ⇒ λ independent of ν; θ_s ∝ √ν ∝ Re^{−1/2}. (c) Closure needs s_s + ℓ ≤ s_TE with ℓ = N_ℓν/Ue_s; in chord units ν = c V/Re. (d) Trapezoid: h(u₀⁵ + u₁⁵)/2; exact for u linear: ∫u⁵ ds = h ∫₀¹(u₀ + (u₁ − u₀)t)⁵ dt = h(u₁⁶ − u₀⁶)/(6(u₁ − u₀)).

## 3.
1. ν = 1.1 × 10⁻⁵/0.017; Re = 120 × 0.08/ν. 2. `section_cd_friction(Re, 0.02, Re_tr=1e12)`. 3. Compare Re with the Re_min values of the lesson; argue from §2.5.

## 4.
1. `bubble_analysis(code, 0.80, 1.643e5)` for three codes. 2. `surface_for_cl("4402", 0.8)` and plot Ue against x. 3. `natural_transition_re(code, cl)` for the six cases (NaN means no change in 10³–10⁹).

## 5.
1. `ladson_blocks(read_ladson())`, `clmax_vs_re` for "I" and "VII"; pair by Re; percentage change. 2. `cd_at_zero_lift` for "VII" and `section_cd_friction(Re, Re_tr=None)`. 3. Argue: at high Re the grit adds roughness drag and thickens the boundary layer that must survive the stall gradient; at low Re forcing transition removes a bubble.

## 6.
(a) Loop over N_ℓ with `bubble_analysis(..., re_bubble=N)`. (b) Free-stream turbulence, the bubble's displacement effect on Ue, three-dimensional and unsteady effects (propeller wash). (c) Howarth: Thwaites separates 2.6 % late; scale to the 0012's separation arc length; Re_min ∝ 1/(s_TE − s_s). (d) Oil-flow or tufts/surface-pressure survey at cruise.
