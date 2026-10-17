# P10.07  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. A1 = cl/π; A2 = A1 + 4cm/π; α_i = A2/3. 2. dz/dx = α_i + A1 cos θ + A2 cos 2θ; z = ½∫₀^θ dz/dx sin θ′ dθ′ term by term. 3. Maximum: evaluate on a fine cosine grid. 4. `thin_airfoil(slope, alpha_i)`: cl, cm_c4, A0, alpha_L0.

## 2.
(a) Substitute the sine series of Δcp/4 into γ/(2V∞), identify Aₙ by orthogonality of sin nθ on [0, π], read dz/dx from (1) with A0 = 0, α − A0 = α_i. (b) ∫₀¹ dz/dx dx = 0 with dx = ½ sin θ dθ; product-to-sum for cos nθ sin θ. (c) Aₙ from the integral of sin nθ; sum the cosine series with the given identity; convert θ to x; integrate x-wise with z(0) = 0. (d) Bound Σ_{n>N} Aₙ Jₙ with |Jₙ| ≲ 1/n.

## 3.
(a) `camber_from_cl_cm(0.3, 0.0)`; evaluate z on a fine grid; find the minimum near the trailing edge and the sign of the slope there. (b) `thicken(cam, slope, 0.12, 100)`; `panel_solve` at α_i and 4°; `naca4("2412")` at 0° and 4°. (c) `thin_airfoil(slope, 0)["alpha_L0"]`; compare cl at a fixed α via 2π(α − α_L0).

## 4.
(a) Build the roof-top function with your a, `inverse_design(roof, 200)`. (b) `thin_airfoil(d["slope"], d["alpha_i"], n_terms=200, m=20000)`, `glauert_load`, maximum difference on 2–98 %. (c) Clip the slope at the ends, `thicken(..., 0.10, 120)`, `panel_solve`, `panel_load` on 10–90 %.

## 5.
(a) Draw w; `decision_matrix(res["raw"], ..., w)`; winner and margin to second. (b) Rows without the 0012; `decision_matrix` with the nominal weights; compare each criterion's min and max with and without it. (c) `weight_flip_points(norm, w, 2)`.

## 6.
(a) Drag at the design cl (`P10.06`, TR-824 polars), cl_max and stall type (`P10.15`, TR-824), Reynolds-number sensitivity (`P10.14`). (b) Its viscous polar and its behaviour at the trainer's Reynolds number; a structural check of the thin trailing edge. (c) Fitting weights to a preferred answer removes the information the matrix was meant to provide; add the missing criteria and report the sensitivity instead.
