# P07.04  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
(a) `heat_transfer_sgen(Q, 340, 290)`. (b) `bodies_equilibrate(1, 900, 300 + 2s, 1, 900, 300)`; C ΔT²/(4T_m²) with C = 900 J/K. (c) `throttle_sgen(p1, 1e5)`; multiply by T₀.

## 2.
(a) Path 1→2 actual, 2→1 reversible; ∮δQ/T ≤ 0 gives ∫δQ/T ≤ S₂ − S₁; name the difference S_gen; for a CV add the entropy carried by mass and set the time derivative to zero for steady flow. (b) Reversible process: δq = T ds, δw = p dv; first law; h = u + pv, dh = du + p dv + v dp; ideal gas substitutions. (c) S_gen = C₁ ln(T_f/T₁) + C₂ ln(T_f/T₂) = −Σ C_i ln(T_i/T_f) ≥ −Σ C_i (T_i/T_f − 1) = 0 because Σ C_i T_i = (Σ C_i) T_f. (d) S_gen = Q_C/T_C − Q_H/T_H ≥ 0 with Q_C = Q_H(1 − η).

## 3.
(a) W/350; Q(1/290 − 1/350); Q(1/3 − 1/290); add. (b) Telescoping sum. (c) Compare magnitudes; consider whether a radiator could reject heat at 3 K; consider which T₀ is meaningful in orbit.

## 4.
(a) `ds_ideal_gas(300, 1e5, T2, p2)`. (b) Three `ds_path_integral` calls with n = 20001. (c) Errors for four n; log₂ of successive ratios. (d) cp ln(T₂/T₁) − R ln(p₂/p₁) with each cp; compare with (a).

## 5.
(a) Edited case dict; `ranked()`; T₀ × each entry. (b) Valve + duct together are R ln(1.0/0.28) per kilogram whatever p_prv is. (c) `heat_exchanger(1.0, 700, T_out, mdot_fan, 320)` for three flows; watch for a `ValueError`. (d) e.g. raise p_prv to the duct need and take bleed from a lower stage (smaller p₁), or replace the valve by a turbine; recompute.

## 6.
(a), (b) `check_cycle` with the given Q and T. (c) Think of streams at very different temperatures or compositions. (d) Separate S_gen (independent of T₀) from T₀ S_gen. (e) Look for a liquid or a two-phase fluid being throttled.
