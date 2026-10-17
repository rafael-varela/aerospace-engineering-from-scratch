# P05.12  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. I = π(D⁴ − d⁴)/64, A = π(D² − d²)/4, EI. 2. P = −6 × 0.3 g; q = −6 ρ A g. 3. `cantilever_point(L, L, L, P, EI)`, `cantilever_udl(L, L, q, EI)`, sum; slope PL²/(2EI) + qL³/(6EI). 4. M_root = PL + qL²/2; σ = |M| (D/2)/I. Compare the deflection with the probe's size and the stress with F_ty.

## 2.
(a) M = q(L − x)²/2; integrate twice; C₁, C₂ from w(0) = w′(0) = 0. Second route: integrate q four times, four constants, conditions w(0), w′(0), M(L) = EIw″(L), V(L) = EIw‴(L). (b) Superposition: q L⁴/(8EI) = R L³/(3EI). Integration: M(x) = q(L − x)²/2 − R(L − x); integrate twice with C₁, C₂; impose w(0) = 0, w′(0) = 0, w(L) = 0; solve for R. The superposition used w(L) = 0 (compatibility at the prop). (c) Tip effect of a point load at η; multiply by w(η)dη and integrate; substitute t = η/s; use the two integrals; prove the second with u = 1 − t².

## 3.
(a) Unknowns: root force, root moment, strut force (3); equations: vertical force and moment (2); one redundant. (b) R from `cantilever_udl(a, s, q, EI) = R a³/(3EI)`. (c) w(x) = `cantilever_udl` + `cantilever_point(x, s, a, −R, EI)`; root moment q s²/2 − R a; moment at the strut q(s − a)²/2; compare with q s²/2 and q s⁴/(8EI). (d) `fd_beam(s, EI, q, n, point_loads=[(a, −R)])` with a on a node; w at that node. (e) Pinned root: moments about the root: R a = q s²/2.

## 4.
1. I from `section_properties(i_section(h, b_f, 0.010, 0.003))` or the closed form; EI = E I. 2. `elliptic_tip(L, s, EI)`; `bending_moment` + `deflection_by_integration`; `fd_beam(s, EI, elliptic, 400)`. 3. Repeat with `tapered_EI(EI, 0.6, s)` for integration and FD. 4. Factor = w_tip(tapered)/(0.05 s).

## 5.
(a) Integrate four times; four conditions; closed form; then `fd_beam(..., bc="clamped-clamped")` for n = 20–160 and `observed_order` of the max error. (b) Two loops over n with a = 0.3 and a = 0.3333; errors against `simply_supported_point`; orders. (c) Two `fd_beam` calls with unit loads; compare the two cross deflections.

## 6.
(a) δ_s ≈ M₀/(G A_web); compare with the bending deflection; it adds deflection. (b) Error of order w′² at the measured slope. (c) Factor = predicted/measured (subtract the shear part first). (d) Subtract the root contribution M₀s/k_θ as well, then recompute the factor. (e) Measure the root rotation directly; measure deflection at several stations to separate root rotation (linear in y) from bending (cubic-like).
