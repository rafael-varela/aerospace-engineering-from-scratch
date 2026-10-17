# P04.02  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. h = 2000 + 100s, V = 120 + s; h_E, K = ½mV², U = mgh. 2. `brentq(lambda ht: ht - zoom_altitude(h, V, 1.2 * stall_speed(ac, ht)), h, h + 5000)`; gain = h_t − h. 3. `ps = specific_excess_power(ac, h, V)`; `climb_rate_with_acceleration(ps, V, 0.02)`.

## 2.
(a) Differentiate, use F = T + D + L + W, L ⊥ v, W conservative → (T − D)V; divide by W. (b) D(q) = qSC_D0 + kW²/(qS); dD/dq = 0 → q²S²C_D0 = kW² → C_D0 = kC_L²; then C_D = 2C_D0 and L/D = C_L/(2C_D0) = 1/(2√(kC_D0)). (c) P = DV with q ∝ V²: P = aV³ + b/V, dP/dV = 0 → 3aV⁴ = b → C_D0 = kC_L²/3; C_L ∝ 1/V² → V_mp/V_md = (1/3)^{1/4}; L/D = C_L/(4C_D0) with C_L = √(3C_D0/k). (d) V₀δV = −gδh (energy); δL = W(2δV/V₀ − βδh) (lift at constant C_L, linear in the perturbations); m δḧ = δL (vertical, γ small, cos γ ≈ 1); combine → δḧ = −(2g²/V₀² + gβ)δh. Evaluate with `density_scale_gradient(3000)` and `phugoid_period_with_gradient`; simulate with `simulate_flight(..., thrust_on=False, drag_on=False)` and `oscillation_period`.

## 3.
(a) ḣ = P/W; Δh = ḣ × 21 600 s. P_s varies with altitude only through the aerodynamic power needed to stay level, which the 150 W "left over" already accounts for at the start; say how it changes. (b) E = PΔt; Wh = J/3600; mass = Wh/250. (c) Build `uav`; at 1000 m and at the top: `min_power_speed`, `stall_speed`, P_r from `power_required_available`, sink = P_r/W; `scipy.integrate.quad(lambda h: 1/sink(h), 1000, top)`. (d) Compare the mass, the energy, and what else the aircraft must carry or do (night length, power to fly level, solar cell area).

## 4.
(a) `grid = np.linspace(energy_height(0, 150), 9000 + 100s, 161)`; `V, h, ps = rutowski_path(ac, grid)`; `time_to_energy(grid, ps)`; read V[−1], h[−1]. (b) `climb_time_accelerate_then_constant(ac, 150, Vc, hE1)` for Vc = 130, 150, 160. (c) `hs = np.linspace(0, h[-1], 161)`; `ps_max = [best_climb(ac, x)[1] for x in hs]`; `np.trapezoid(1/ps_max, hs)`. Missing: the kinetic energy at the end differs from that at the start. (d) `replace(ac, mass=1.1*ac.mass)`; same steps.

## 5.
(a) Loop n; bank = degrees(arccos(1/n)). (b) `brentq` on n of P_s; then g√(n² − 1)/V in rad/s → deg/s, and V²/(g√(n² − 1)). (c) Define n_max(V) by the same `brentq`; `minimize_scalar(lambda V: -n_max(V), bounds=(110, 0.75*speed_of_sound(5000)), method="bounded")`. (d) q at 5 km, n = qSC_Lmax/W; Δh_E = P_s(n = 4) × 5 s; at constant h, ΔV ≈ gΔh_E/V.

## 6.
Check where m appears: h_E (nowhere), P_s (through C_L in the induced drag and through the division by W), L/D (nowhere, at given C_L), V at given C_L (through W/S). Recall that the energy method assumes instant exchange along h_E lines. Compare a 10⁻⁸ m energy error with the other uncertainties in the model.
