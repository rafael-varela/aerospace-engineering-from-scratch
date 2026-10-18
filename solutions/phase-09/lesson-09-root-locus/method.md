# P09.09  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
T = 0.8 + s/50; `num, den = altitude_plant(t_path=T)`. (a) Poles 0, −1/T, −2; segments by counting; `asymptotes`. (b) `breakaway_points`; K = −D/N there. (c) den = Tτ s³ + (T + τ)s² + s; Routh s¹ entry vanishes at 50K = (T + τ)/(Tτ); ω² = 1/(Tτ). (d) `gain_for_damping(num, den, 0.5)`; `np.roots(char_poly(num, den, K))`; `step_response(K*num, char_poly(num, den, K), t)` and `overshoot_and_settling`; predicted 16.3 % and 4/(ζω_n).

## 2.
(a) 1 + KG = 0 ⇔ G = −1/K; for K > 0 real, −1/K has angle 180° and magnitude 1/K. (b) Angle of (s − r) for real s and real r is 0 if r < s, 180° if r > s; a conjugate pair contributes θ + (−θ) = 0. (c) D/N = s^{n−m} − (Σp − Σz)s^{n−m−1} + … = (s − σ_a)^{n−m} + O(s^{n−m−2}); then (s − σ_a)^{n−m} = −K with K → ∞. (d) On the real axis between two poles K(σ) = −D/N has an extremum: dK/dσ = −(D′N − DN′)/N² = 0. Root sum: −(coefficient of s^{n−1})/(leading coefficient) of D + KN; KN has degree m ≤ n − 2, so it does not touch s^{n−1}.

## 3.
den = (s² − 4)(0.1s + 1) = 0.1s³ + s² − 0.4s − 4; char 0.1s³ + s² + (K − 0.4)s + (K − 4). Segments by counting (poles 2, −2, −10; zero −1). Routh: all coefficients positive and 1·(K − 0.4) > 0.1(K − 4). `jw_crossings`. `gain_for_damping(num, den, 0.5)`; roots.

## 4.
`pd_design_table((0.5 + s/10, 1.0, 3.0))`; compare the third pole with −1/T_d and with the pair's real part.

## 5.
Rows 1 and 2 from alternating coefficients; fill; apply ε for a lone zero in the first column; for a zero row, differentiate the auxiliary polynomial of the row above. Count sign changes; compare with `routh_rhp_count` and `np.roots`.

## 6.
(a) `altitude_plant(tau=…)` and `altitude_plant(v=…)` with `np.roots(char_poly(num, den, 0.01037))`; ζ = −Re/\|p\| of the complex pair; `jw_crossings`. (b) Noise on ḣ is multiplied by K T_d and fed to the pitch command. (c) 5 s settling needs ζω_n ≥ 0.8; 10 % overshoot needs ζ ≥ 0.59; check whether the P-control locus has any point satisfying both (compare the real part of the pair along the locus with −0.8).
