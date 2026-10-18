# P10.22  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
Sort each item by what it compares with: an exact property (code verification), another grid of the same code (solution verification), or a measurement (validation). Item (f) is about the data's own uncertainty.

## 2.
(a) f₃ − f₂ = C h₁ᵖ(r²ᵖ − rᵖ), f₂ − f₁ = C h₁ᵖ(rᵖ − 1); divide; monotone convergence needed (same sign). Richardson from f* = f₁ − C h₁ᵖ. (b) GCI as a 95 % band on the fine-grid error; divide by the coverage factor 2. (c) Substitute c′ = κc and a₀′ = 2π; S′ = κS, AR′ = AR/κ; CL′ = κ⁻¹ CL... collect.

## 3.
1. κ = a₀(2.0e6)/2π. 2. `grid_study(8.333, 1.0, 0.0, kappa=κ)`. 3. Sensitivity by central difference in κ on the fine lattice; u_input = |∂CL_α/∂κ| u_κ. 4. Combine u_num and u_input; U = 2u; per degree too. 5. Compare with the validation's mean signed error and scatter.

## 4.
1. Build the table with `validation_table(κ, u_κ)`; recompute u_val for each scatter with `validation_uncertainty(..., scatter_rel=...)`; refresh "within". 2. `validation_table(1.0, 0.0)` and `validation_metrics`. 3. Argue from what u_val contains.

## 5.
1. `vlm_cla(AR, taper, sweep, nc, ns_half, κ)` for (nc, ns_half) = (2, 6), (4, 12), (8, 24), (16, 48). 2. `observed_order_3`, `gci_fine`, `richardson` on the first three. 3. Compare with the 16 × 96 value.

## 6.
Argue from the data README's caveats, the model's assumptions (flat wake, linear), and the section-slope proxy; write the conclusion with numbers and scope.
