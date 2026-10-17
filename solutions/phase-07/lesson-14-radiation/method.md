# P07.14  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. F₁₃ = `vf_coaxial_disks(r, r, h)`. 2. F₁₁ = 0 ⇒ F₁₂ = 1 − F₁₃. 3. Reciprocity: F₂₁ = A₁F₁₂/A₂ with A₂ = 2πrh. 4. Symmetry: F₂₃ = F₂₁. 5. Summation: F₂₂ = 1 − F₂₁ − F₂₃. 6. Row 3 mirrors row 1. Check rows and A_iF_ij.

## 2.
(a) Write A₁F₁₂ as the double integral; swapping 1 and 2 leaves the integrand unchanged. ΣQ_i = Σ_i A_iJ_i − Σ_iΣ_j A_iF_ijJ_j; replace A_iF_ij by A_jF_ji, sum over i using Σ_i F_ji = 1. (b) G = J − Q/A from the definition; substitute into J = εE_b + (1 − ε)G; for known T use G = ΣF_ijJ_j instead. (c) Q₁ = (E_b1 − J₁)ε₁A₁/(1 − ε₁), Q₁ = A₁F₁₂(J₁ − J₂), Q₁ = (J₂ − E_b2)ε₂A₂/(1 − ε₂); add the three potential differences. Limits: F = 1, A₁ = A₂; then A₁/A₂ → 0 with F₁₂ = 1. (d) Each gap is a parallel-plate pair; the same Q passes through all gaps in series; add their resistances.

## 3.
(a) h(T_g − T_b) = εσ(T_b⁴ − T_w⁴); `brentq` on [T_w, T_g]. (b) Change ε or h and repeat. (c) Shield: 2h_s(T_g − T_s) = ε_sσ(T_s⁴ − T_w⁴) (two convecting faces, one radiating face toward the wall); then the bead with T_s in place of T_w. (d) Compare errors and practical problems (coating durability, h depends on flow, shield fouling and its own lag).

## 4.
1. Radiator: o = 0, u = (0, 1, 0), v = (0.6, 0, 0). Array: o = 0, u = (0, 1, 0), v = 1.0 (cos φ, 0, sin φ). 2. `vf_numerical` both ways at each n; compare A F. 3. F_∞ = 2F₄₀ − F₂₀; repeat for φ = 90° against the closed form. 4. Build F for radiator, array, space (large area, reciprocity for its row), call `solve_enclosure` with Q = 200 W on the radiator, T on the others.

## 5.
1. `box_view_factors(1, 1, 1)`; make face 5 black at 0 K, faces 0–4 at 400 K with ε_w. 2. `solve_enclosure`; Q_out = −Q₅. 3. ε_eff = Q_out/(A₅σT⁴). 4. Repeat for other ε_w and for `box_view_factors(1, 1, 3)` (the open face is a 1 × 1 end). 5. Explain with multiple reflections.

## 6.
(a) Compare the linearised and exact fluxes at 20 °C. (b) Compare the absorbed solar flux with α = 0.93 and 0.15 against the emitted flux. (c) `radiator_temperature` of the 1.079 m² panel with α_s = 0.30 in normal sun. (d) The lesson's array result (F 0.269–0.294 and its temperature rise). (e) Ideal shields against real blankets; compute the flux with an effective emittance typical of real blankets (label it illustrative).
