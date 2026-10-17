# P07.12  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. t_cu, h_c from s. 2. `wall_heat_flux([(0.001, 11.6), (t_cu, 398)], 3000, 300, 3000, h_c)`; the resistances and shares by hand. 3. Repeat with copper only; compare the shares.

## 2.
(a) Four neighbour conductances k t Δy/Δx and k t Δx/Δy; divide by V = ΔxΔy. (b) Σ_i (K T)_i over all nodes: conduction terms cancel in pairs; free rows are zero; what is left is the fixed-node residuals, face and edge terms. (c) Energy balance of a strip of width W and length dx. (d) ξ = 1 − 4 Fo sin²(kΔx/2); |ξ| ≤ 1 ⇔ Fo ≤ ½.

## 3.
1. n = 121 nodes over 12 mm (interface at node 100); k per interval; ρc per node with the interface averaged; `conduction_system_1d(..., ("flux", 1e5), ("adiabatic",))`; `march(theta=0.5)` with Δt = 0.05 s. 2. e = √(kρc) of silica; formula; compare with √(αt) against 10 mm. 3. `E_stored`, `Q_in` against q″ t.

## 4.
1. P_chip, T_wall from s; three grids. 2. `brentq` on `t_frame` for T_max = 373.15 K. 3. `heat_balance`; the overlap-area generation.

## 5.
1. Reference march; four Δt per scheme; `observed_order`. 2. Fo = αΔt/Δx²; amplification factors of the two schemes for stiff modes. 3. Accuracy needed in time against cost.

## 6.
1. `time_to_surface_limit` for two h; ratio. 2. Lower k lowers effusivity; time ∝ e² roughly for a semi-infinite convective wall at fixed βt-group. 3. The thickness sweep of the lesson. 4. Rank by sensitivity.
