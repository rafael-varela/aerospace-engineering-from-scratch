# P05.07  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. A = π(25.4² − 22.2²)/4. 2. σ = −P/A; δ = `axial_deformation(-P, 600, A, E)`; |σ|/Fty. 3. `constrained_thermal_stress(E, alpha, 30)`. 4. ΔT = Fty/(Eα).

## 2.
(a) δ = (P/E)∫dx/A(x) = (PL/(E(A₁ − A₀))) ln(A₁/A₀); limit by l'Hôpital or ln(1 + x) ≈ x. (b) Elongations u and −u; forces N₁ = k₁(u − α₁L₁ΔT), N₂ = k₂(−u − α₂L₂ΔT); junction N₂ − N₁ + P = 0; solve for u; the result is linear in (P, ΔT). (c) Compliance (1/(πE))∫ds/((u − a)(u + a)) with du = tan α ds; partial fractions; evaluate between s = 0 and t. (d) Solve (8), (9) for F_b, F_m; set F_m = 0; for the thermal case, write the free misfit as the sum of the bolt's extra stretch and the members' extra squeeze with a common force change ΔF.

## 3.
(a) ΔT = −203 K; σ = −EαΔT; force σA; compare with Fty. (b) Free shrinkage δ = α L |ΔT|; line stiffness k_l = EA/L; force F = δ/(1/k_l + 1/k_bel); require F ≤ (Fty/3)A, solve for k_bel; bellows stroke = δ − F/k_l. (c) Nodes at 0, L, L + 1; elements (line, EA = k_l L) and (bellows, EA = k_bel × 1); fixed nodes 0 and 2; eps0 = [αΔT, 0]. (d) α of aluminium falls at cryogenic temperature; E rises; the yield strength rises.

## 4.
(a) For each T: F_i = `preload_from_torque(T, K, 8)`; add ΔF_cold and ΔF_hot from `LESSON_JOINT.preload_at`; fraction with `separation_load(...) < 18e3`; fraction with `joint_loads(..., 12e3, ...)[0] / A_t > 880`. (b) F_i,min = (1 − C)·18 kN − ΔF_cold; F_i,max = 880·A_t − ΔF_hot − C·12 kN. (c) Same with sd 0.01; `np.percentile(F_i, [1, 99])`.

## 5.
(a) N(r) = ρAω²((r₀ + L)² − r²)/2; at r = r₀ this is ρAω²(r₀L + L²/2); divide by A. δ = ∫N/(EA) dr. (b) Nodes `np.linspace(0, 100, n + 1)`; q at each node; assemble consistent loads; fixed node 0; tip u[-1]; root stress N[0]/A. Compare tip with (a); fit the order of the root-stress error. (c) αΔT L; add to the centrifugal stretch.

## 6.
(a) Separation load ∝ F_i; the bolt sees C times the alternating load only while clamped; preload loss by embedding and temperature eats a margin. (b) `bolt_stiffness(8, 1.25, 20.0, 0.0, E_Ti)`, then C. Higher C: the bolt sees more of the alternating load (worse fatigue) and the series stiffness in (12) rises (more cold loss). (c) Near an edge the cone is truncated; k_m falls; C rises. (d) Stretch measurement gives F_i = k_b × stretch; with the applied torque, K = T/(F_i d).
