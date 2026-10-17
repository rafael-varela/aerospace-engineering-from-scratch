# P08.12  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. c from s; Re_c = Uc/ν. 2. `flat_plate_transition(30, ν, 9.0)`: x0, x_tr; compare x_tr with c. 3. N_crit = `mack_ncrit(0.003)`; recompute x_tr/c. 4. Flight is quieter than most tunnels.

## 2.
(a) u = U + u′, v = v′, p = P + p′; drop products of primes; eliminate p′; substitute ψ′; divide by iα; non-dimensionalise with U_max and h. (b) Rayleigh: (U − c)(φ″ − α²φ) − U″φ = 0; divide by U − c, multiply by φ*, integrate by parts using φ = 0 at the walls; imaginary part: c_i ∫U″|φ|²/|U − c|² dy = 0. (c) e^N; amplitude ratio × initial amplitude.

## 3.
1. h from s; ν = μ/ρ. 2. U_c = 5772.22 ν/h; U_b = 2U_c/3. 3. Re_b = U_b 2h/ν and U_b 4h/ν. 4. Compare with the qualitative observation and explain subcritical (bypass) transition.

## 4.
(a) Loop over N with `neutral_reynolds(α₁, N)`. (b) Loop over α; `critical_point(60)`. (c) `os_eigenvalues(60, 1, 1e4)`; count c_i > 0; repeat with N = 100 and match eigenvalues.

## 5.
(a) `mack_ncrit`; `flat_plate_transition(60, ν, N_crit)`. (b) `transition_reynolds(H, 0.664115, 9)` for each H; x_tr = Re_x,tr ν/U. (c) `thwaites_wedge(m)` → (k, H) → `transition_reynolds(H, k, 9)` with Ue ≈ 60 m/s.

## 6.
(a) Bypass transition (high Tu, wakes), roughness (insects, rivets, steps and gaps), crossflow instability, attachment-line contamination. (b) Sweep: crossflow profiles are three-dimensional. (c) Infrared thermography or hot films in flight; transition well ahead of 60 % contradicts it.
