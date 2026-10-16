# P03.07  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. x = (s/20)². M = [[1 − x/2, h], [−hω²(1 − x/4), 1 − x/2]] (write with ω = 1, h = s/20). 2. Coefficient 1 − x/4; relative swing x/4. 3. ω_h/ω = (2/√x) asin(√x/2); error ω_h/ω − 1 against x/24; steps per period 2π/(hω); limit π.

## 2.
(a) J = [[0, I], [−I, 0]]. Kick M = [[I, 0], [hA, I]]: MᵀJM = [[hAᵀ − hA, I], [−I, 0]] = J iff A = Aᵀ. Drift M = [[I, hI], [0, I]]: MᵀJM = J identically. A = ∂a/∂q = −∂²V/∂q² is symmetric. Products of symplectic matrices are symplectic (MᵀJM = J and NᵀJN = J give (MN)ᵀJ(MN) = Nᵀ(MᵀJM)N = J). (b) Kick: q × (p + h f q) = q × p + h f q × q = q × p. Drift: (q + hp) × p = q × p. (c) Verlet = K_{h/2} ∘ D_h ∘ K_{h/2}; the inverse of K_h is K_{−h} and of D_h is D_{−h}; the composition with −h is K_{−h/2} D_{−h} K_{−h/2}, which is the inverse of the original in reverse order: symmetric. Symplectic Euler = D_h ∘ K_h; with −h it is D_{−h} K_{−h} ≠ (D_h K_h)⁻¹ = K_{−h} D_{−h}: not symmetric. N steps: D_h K_h D_h K_h ⋯ D_h K_h; insert D_{h/2} D_{−h/2} = I at the start and D_{−h/2} D_{h/2} between every K_h and the next D_h... regroup as D_{h/2} (D_{h/2} K_h D_{h/2})^N D_{−h/2}.

## 3.
(a) v_p = √(μ(1 + e)/r_p), v_a = v_p r_p/r_a; ω = v/r; ratio = (r_a/r_p)² = ((1 + e)/(1 − e))². (b) N_p = P ω_p/0.05, N_a = P ω_a/0.05; waste = 1 − N_a/N_p. (c) `energy_drift(stormer_verlet_step, q0, p0, 5, 2000)` etc. with `perigee_state(26600e3, 0.74)`. Compare swings; apply the (hω_p)² and (hω_p)⁴ scalings to 1 000 orbits and to RK4's drift. (d) dt = r dτ (Sundman); the transformed system must remain Hamiltonian (it does, with the extended Hamiltonian r(H − E₀)) for a symplectic method in τ to keep its properties.

## 4.
(a) a = r_p/(1 − e); P = 2π√(a³/μ); v_p = √(μ(1 + e)/r_p); ω_p = v_p/r_p; hω_p = ω_p P/200. (b) `energy_drift(..., 50, 200)`, `drift_slope`, `angular_momentum` along `propagate`, `kepler_exact(50 P)`. (c) Ratio of (hω_p)² values; the measured ratio exceeds it because the perigee passage sharpens with e (the force's third derivative along the path grows faster than ω_p²). (d) Compare max |ΔE/E| and the sign of the slope; extrapolate RK4's slope to 5 000 orbits.

## 5.
(a) `energy_drift(stormer_verlet_step, q0, p0, 20, 800)` → swing; `energy_drift(rk4_step, q0, p0, 100, 200)` → `drift_slope`; N* = swing/|slope|. (b) `energy_drift(yoshida4_step, q0, p0, 20, 200)`. (c) Position error = phase along the ellipse (grows ∝ N for both); energy error = size/period of the ellipse (bounded vs secular). (d) Verlet 800: swing stays; RK4 200: |slope| × 10⁴; Yoshida 200 (600 evaluations, under budget): swing stays; pick by the two numbers and the position requirement.

## 6.
(a) Adaptive RK: burns and a flyby are events with widely varying time scales; not Hamiltonian during burns. (b) Symplectic (with J₂ and lunisolar terms, all conservative): 10⁵ orbits, secular drift would be mistaken for physics. (c) Adaptive RK: drag dominates, nothing is conserved, 3 minutes long. (d) Symplectic (Wisdom–Holman class): 10⁸ orbits. (e) Either, with the drag model dominating the error budget; a symplectic method gains nothing, an adaptive one handles the slow decay well. (f) Verlet: 10⁶ steps of a Hamiltonian system where a bounded energy is the thermodynamic requirement and force evaluations are the whole cost.
