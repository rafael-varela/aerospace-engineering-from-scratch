# P03.07  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for $s = 15$; numerical results are from `solution.py` and given to the figures shown.

## 1. Verlet on the oscillator ($s = 15$: hω = 0.75)

**Derivation.** With ω = 1 and h = 0.75, x = (hω)² = 0.5625.

(a) M = [[1 − x/2, h], [−hω²(1 − x/4), 1 − x/2]] = **[[0.71875, 0.75], [−0.64453, 0.71875]]**; det M = 0.71875² + 0.75 × 0.64453 = 0.51660 + 0.48340 = **1.00000** (exactly: (1 − x/2)² + x(1 − x/4) = 1).

(b) Coefficient 1 − x/4 = **0.859375**; the true energy swings by x/4 = **14.06 %** of itself (peak to peak, relative), while H̃ is constant.

(c) ω_h/ω = (2/0.75) asin(0.375) = 2.6667 × 0.38440 = **1.02506**, a frequency error of **+2.51 %**; the leading term x/24 = **2.34 %** (the next term, 3x²/640, adds 0.15 %). The step is 2π/0.75 = **8.4 steps per period**; the stability limit hω = 2 is **π = 3.14 steps per period**.

**Check.** The eigenvalues of M are e^{±iω_h h}: cos(ω_h h) = 1 − x/2 = 0.71875 = cos(0.7688), and 0.7688/0.75 = 1.02506 ✓.

**Common mistakes.** Writing the modified-energy coefficient with a plus sign; quoting the energy swing as an error in the *solution* (the trajectory is bounded and its frequency is wrong by 2.5 %, which is the real error).

## 2. Three properties of kick-drift-kick

(a) J = [[0, I], [−I, 0]]. For the kick, M = [[I, 0], [hA, I]] with A = ∂a/∂q: MᵀJM = [[I, hAᵀ], [0, I]] [[0, I], [−I, 0]] [[I, 0], [hA, I]] = [[I, hAᵀ], [0, I]] [[hA, I], [−I, 0]] = [[hA − hAᵀ, I], [−I, 0]], which equals J iff **A = Aᵀ**. For a = −∇V, A = −∂²V/∂q∂q is symmetric (equality of mixed partials). For the drift, M = [[I, hI], [0, I]]: MᵀJM = [[I, 0], [hI, I]] [[0, I], [−I, 0]] [[I, hI], [0, I]] = [[I, 0], [hI, I]] [[0, I], [−I, −hI]] = [[0, I], [−I, hI − hI]] = J, with no condition. If MᵀJM = J and NᵀJN = J then (MN)ᵀJ(MN) = Nᵀ(MᵀJM)N = NᵀJN = J, so any composition of kicks and drifts, with any step sizes including negative ones (Yoshida's w₀), is symplectic.

(b) For a(q) = f(|q|) q, a kick gives q × (p + h f q) = q × p + h f (q × q) = q × p; a drift gives (q + h p) × p = q × p + h (p × p) = q × p. Each sub-step conserves L exactly, so every composition does, at every h. (RK4 evaluates the force at intermediate points that are not the stored q, and its p-update is not parallel to the stored q: no such argument exists, and the test shows 10⁻⁷ drift.)

(c) Verlet V_h = K_{h/2} D_h K_{h/2}. Its inverse is K_{−h/2} D_{−h} K_{−h/2} (reverse the order, negate each step), which is V_{−h}: applying the step with −h to the output returns the input exactly. Symplectic Euler SE_h = D_h K_h has inverse K_{−h} D_{−h} ≠ D_{−h} K_{−h} = SE_{−h}, so it is not symmetric. For N steps, SE_h^N = D_h K_h D_h K_h ⋯ D_h K_h. Insert D_{h/2} D_{−h/2} = I at the front and regroup each D_h in the interior as D_{h/2} D_{h/2}: SE_h^N = D_{h/2} (D_{h/2} K_h D_{h/2}) (D_{h/2} K_h D_{h/2}) ⋯ (D_{h/2} K_h D_{h/2}) D_{−h/2} = **D_{h/2} Ṽ_h^N D_{−h/2}**, with Ṽ_h the drift–kick–drift Verlet step. The O(h) error of symplectic Euler is the fixed map D_{h/2} (a shift of q by hp/2) applied once at the end and undone at the start; after a whole period Ṽ_h^N is close to the identity, the two shifts nearly cancel, and only Ṽ's O(h²) error remains. The code measures order 2.00 at t = P and 0.94 at t = 1.5P.

**Common mistakes.** Using Ṽ = K D K in (c) and concluding the conjugation is by a kick (it is by a drift, because SE_h ends with a drift); forgetting that negative sub-steps are still symplectic (needed for Yoshida).

## 3. Molniya orbit with a fixed step

**Assumptions.** Two-body, μ = 3.986 × 10¹⁴ m³/s²; r_p = 6 916 km, r_a = 46 284 km, P = 43 175 s = 11.99 h.

**Derivation.** (a) v_p = √(μ(1 + e)/r_p) = **10 014 m/s**; v_a = v_p r_p/r_a = **1 496 m/s**; ω_p = v_p/r_p = **1.448 × 10⁻³ rad/s**; ω_a = v_a/r_a = **3.233 × 10⁻⁵ rad/s**; ratio **44.79**. By conservation of angular momentum v_p r_p = v_a r_a, so ω_p/ω_a = (v_p/r_p)(r_a/v_a) = (r_a/r_p)² = ((1 + e)/(1 − e))² = (1.74/0.26)² = 44.79 ✓.

(b) N_p = P ω_p/0.05 = **1 250 steps per orbit**; N_a = P ω_a/0.05 = **27.9**; the fraction wasted on the slow part is 1 − 27.9/1 250 = **97.8 %**.

(c) Five orbits:

| steps/orbit | hω_p | Verlet max \|ΔE/E\| | RK4 max \|ΔE/E\| |
|---|---|---|---|
| 2 000 | 0.031 | 4.70 × 10⁻⁴ | 2.87 × 10⁻⁸ |
| 8 000 | 0.0078 | 2.94 × 10⁻⁵ | 7.4 × 10⁻¹¹ |

Verlet's swing scales as (hω_p)² (ratio 16.0), RK4's error as (hω_p)⁴ (ratio 390). RK4's drift at 2 000 steps is about −2.9 × 10⁻⁹ per orbit, so after 1 000 orbits its energy error is about 3 × 10⁻⁶, still 150 times below Verlet's bounded 4.7 × 10⁻⁴. For this satellite at a fixed step, **RK4** (or better, Yoshida at 2 000 steps, which has both the fourth power and the bound). The answer differs from the e = 0.1 case because the fixed step is set by the perigee rate, which is 45 times the apogee rate: hω_p is large wherever the cost is affordable, and the second-order method's h² bound is then worse than the fourth-order method's secular h⁴ t drift for any propagation shorter than N* ≈ 10⁵ orbits.

(d) Use the Sundman transformation dt = r dτ and integrate in τ: equal steps in τ are steps of physical length h r, short at perigee and long at apogee, with no error estimate and no feedback from the solution. For a symplectic method in τ to keep its properties the transformed equations must still be Hamiltonian, which they are for the extended Hamiltonian Γ = r (H − E₀) with E₀ the initial energy.

**Check.** P from a: 2π√(a³/μ) = 2π√(1.882 × 10²³/3.986 × 10¹⁴) = 43 175 s ✓.

**Common mistakes.** Setting the step from the mean motion n = 2π/P (which gives hω_p = 0.5 at 200 steps per orbit and a 5 % energy swing); concluding "symplectic is always better".

## 4. Your own eccentricity ($s = 15$: e = 0.25)

**Derivation.** (a) a = 6 300/0.75 = **8 400 km**; P = 2π√(a³/μ) = **7 661.8 s** (127.7 min); v_p = √(μ × 1.25/6.3 × 10⁶) = **8 893 m/s**; ω_p = **1.4116 × 10⁻³ rad/s**; hω_p = ω_p P/200 = **0.0541**.

(b) 50 orbits at 200 steps per orbit:

| method | max \|ΔE/E\| | final ΔE/E | slope per orbit | max \|ΔL/L\| | position error |
|---|---|---|---|---|---|
| Verlet | 3.73 × 10⁻⁴ | +2.2 × 10⁻⁵ | +1.8 × 10⁻⁸ | 1 × 10⁻¹⁴ | 2 430 km |
| RK4 | 1.38 × 10⁻⁶ | −1.35 × 10⁻⁶ | −2.70 × 10⁻⁸ | 4.3 × 10⁻⁷ | 3.7 km |
| Yoshida | 1.10 × 10⁻⁶ | +2 × 10⁻¹⁴ | −1.4 × 10⁻¹⁴ | 1.8 × 10⁻¹⁴ | 6.6 km |

(c) (hω_p)² ratio: (0.0541/0.0386)² = **1.96**; measured 3.73 × 10⁻⁴/1.06 × 10⁻⁴ = **3.5**. The energy swing of Verlet is governed by the h² term of the modified Hamiltonian, which involves the second derivative of the force along the path, not only ω_p²; at higher e the perigee passage is sharper (the force changes faster relative to the step) and the h² term grows faster than ω_p² does. The (hω_p)² rule is a lower bound, good for scaling at fixed e.

(d) RK4's final error, 1.35 × 10⁻⁶, is smaller than Verlet's swing, 3.7 × 10⁻⁴, by 270 times. But Verlet's slope is noise (1.8 × 10⁻⁸ per orbit on an oscillation of 3.7 × 10⁻⁴, with a final value within the swing), while RK4's final value *equals* slope × 50 = −1.35 × 10⁻⁶: a secular drift. At 5 000 orbits RK4 would be at −1.35 × 10⁻⁴, comparable to Verlet's bound and still growing; Yoshida would still be at 10⁻⁶. The number that supports the choice is the slope: −2.7 × 10⁻⁸ per orbit for RK4 against 10⁻¹⁴ for Yoshida.

**Check.** The same 50-orbit run with the perigee at 7 000 km × 0.75 instead of 6 300 km gives identical relative numbers: at fixed e and fixed steps per orbit the two-body problem is scale-free.

**Common mistakes.** Letting a = 7 000 km at e = 0.25, which puts the perigee at 5 250 km, inside the Earth; reading Verlet's fitted slope as a drift.

## 5. Matched cost and the crossover

(a) Verlet at 800 steps per orbit, 20 orbits: swing **6.59 × 10⁻⁶**. RK4 at 200 steps per orbit, 100 orbits: slope **−7.28 × 10⁻⁹ per orbit** (final −7.28 × 10⁻⁷: linear). N* = 6.59 × 10⁻⁶/7.28 × 10⁻⁹ = **905 orbits**.

(b) Yoshida at 200 steps per orbit (600 evaluations): swing **2.28 × 10⁻⁷**, bounded. It sits below both from the first orbit, at three quarters of the cost.

(c) The position error is the phase along the orbit, which grows linearly with time for every method (period error × number of orbits); at 200 steps RK4's per-orbit phase error (17 m) is smaller than Yoshida's (42 m) because its error constant is smaller on this problem. The energy error is the size and period of the ellipse. RK4's phase error will *also* become worse than Yoshida's, once its energy drift has changed the period enough: with ΔE/E = 7.3 × 10⁻⁹ N, the period error is 1.5 × 7.3 × 10⁻⁹ N and the accumulated phase grows as N², overtaking Yoshida's linear 42 m per orbit at N ≈ 2 × 42/(1.5 × 7.3 × 10⁻⁹ × 4.4 × 10⁷ m) ≈ 170 orbits.

(d) At 10⁴ orbits and 800 evaluations per orbit: Verlet (800 steps) ends at |ΔE/E| ≤ **6.6 × 10⁻⁶**; RK4 (200 steps) at about **7.3 × 10⁻⁵** and growing; Yoshida (200 steps, 600 evaluations, or 267 steps for the full budget) at **2.3 × 10⁻⁷** bounded (8 × 10⁻⁸ at 267 steps). Run Yoshida.

**Common mistakes.** Comparing at equal steps instead of equal cost; forgetting that RK4's position error inherits its energy drift quadratically.

## 6. When symplectic is the wrong tool

(a) **Adaptive RK** (DOP853 with events): burns are non-Hamiltonian, the flyby needs steps that vary by 100× across the trajectory, and 5 days is a few orbits.
(b) **Symplectic** (Verlet or Yoshida, with J₂ and lunisolar perturbations, all conservative): 30 years is 1.5 × 10⁵ orbits, where a secular energy drift of RK4 would masquerade as orbital decay.
(c) **Adaptive RK**: 3 minutes, drag-dominated, nothing conserved, a sharp peak that needs step control.
(d) **Symplectic** (a Wisdom–Holman mixed-variable method): 10⁸ orbits; only the bounded-error property keeps the question answerable at all.
(e) **Either**, and preferably adaptive: drag removes the conservation argument; the error budget is the atmosphere model (tens of percent), not the integrator.
(f) **Verlet**: a Hamiltonian system where a bounded energy *is* the physical requirement (the ensemble temperature), the force evaluation is the whole cost, and 10⁶ steps of RK4 would both cost four times more and heat or cool the system artificially.

## Check-yourself answers

1. p_{n+½} = p_n + (h/2)a(q_n); q_{n+1} = q_n + h p_{n+½}; p_{n+1} = p_{n+½} + (h/2)a(q_{n+1}). The velocity (momentum) is evaluated at the half step.
2. det M = 1: the map preserves area in the phase plane, so a region of initial conditions keeps its area; explicit Euler expands it (det = 1 + h²ω²) and RK4 contracts it (det = |R(iωh)|² < 1).
3. H̃ = p²/2 + (ω²/2)(1 − h²ω²/4)q², constant to round-off; the true energy oscillates with relative amplitude h²ω²/4.
4. A symplectic method follows exactly the flow of a modified Hamiltonian H̃ = H + O(h^p); H̃ is conserved and H − H̃ is bounded, so H's error is bounded. RK4 has no such H̃, and its per-step energy error accumulates linearly in time.
5. A kick changes p along q (central force) and a drift changes q along p; neither changes q × p, so the composition conserves it exactly.
6. w₁ = 1/(2 − 2^{1/3}) = 1.3512 for the outer sub-steps, w₀ = 1 − 2w₁ = −1.7024 for the middle one; the h³ error term of Verlet can only be cancelled by a combination w₀³ + 2w₁³ = 0 with w₀ + 2w₁ = 1, which forces one weight negative.
7. N symplectic Euler steps equal N drift-kick-drift Verlet steps conjugated by a half-step drift; after a whole period the Verlet part returns near the start and the two half-drifts cancel, leaving Verlet's O(h²) error.
8. A re-entry or a thrust arc: the force is not conservative (drag, thrust), so there is no energy to conserve, the time scales vary sharply, and an adaptive Runge–Kutta method with step control is both more accurate and cheaper.
