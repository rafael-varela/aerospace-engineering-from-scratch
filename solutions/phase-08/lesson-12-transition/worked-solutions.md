# P08.12  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and probe scripts that call it, given to the figures shown. Transition points come from the supplied **approximate-fit** chart (Drela and Giles 1987) and Mack's **correlation**; they are estimates, not measurements.

## 1. A sailplane wing (s = 15: c = 0.90 m)

**Assumptions.** Flat plate (zero pressure gradient), Blasius layer (H = 2.591, Re_θ = 0.664√Re_x), N_crit = 9 for quiet air, smooth surface.

(a) Re_c = 30 × 0.90/1.4607 × 10⁻⁵ = **1.85 × 10⁶**. The chart gives Re_θ0 = 242, so instability starts at Re_x0 = (242/0.664)² = 1.33 × 10⁵, **x₀ = 64.6 mm**; transition at Re_θ = 1108, Re_x = 2.79 × 10⁶, **x_tr = 1.36 m**, beyond the 0.90 m chord: **yes, laminar to the trailing edge** (in this model).

(b) Mack: N_crit = −8.43 − 2.4 ln(0.003) = **5.51**; Re_θ,tr = 242 + 5.51/0.01039 = 772; x_tr = **0.659 m = 0.73c**.

(c) (a) for the aircraft: the atmosphere in smooth gliding conditions is much quieter than most tunnels, which is why N_crit ≈ 9 is used for flight. It matters for a sailplane because its performance (glide ratio) depends heavily on how much of the wing stays laminar; a tunnel with 0.3 % turbulence would under-predict the laminar run by a quarter of the chord and mislead a design comparison unless N_crit is matched to the tunnel.

**Check.** x_tr scales with ν/U at fixed N_crit: halving the speed doubles x_tr.

**Common mistakes.** Using Re_c instead of Re_x; forgetting that a real sailplane section has a pressure gradient that changes H (Problem 5c).

## 2. Derivations

(a) Write u = U(y) + u′, v = v′, p = P(x) + p′ and subtract the base flow (which satisfies the equations with a constant dP/dx balancing νU″). Dropping products of primes (small disturbances):
u′_t + U u′_x + v′U′ = −p′_x/ρ + ν∇²u′,  v′_t + U v′_x = −p′_y/ρ + ν∇²v′,  u′_x + v′_y = 0.
Take ∂/∂y of the first minus ∂/∂x of the second to remove p′; with ψ′ = φ(y)e^(iα(x − ct)), u′ = φ′e^(…), v′ = −iαφe^(…), and ∇² → D² − α²:
iα[(U − c)(φ″ − α²φ) − U″φ] = ν(φ⁗ − 2α²φ″ + α⁴φ);
non-dimensionalising with the centreline speed and the half-height (Re = U_c h/ν), **(U − c)(φ″ − α²φ) − U″φ = (φ⁗ − 2α²φ″ + α⁴φ)/(iαRe)**. No slip and no penetration give φ = φ′ = 0 at a wall. Neglected: the products of disturbances (nonlinear terms) and, for a boundary layer, the slow growth of the layer (V ≠ 0 and dU/dx ≠ 0 are O(1/√Re_x)); the **parallel-flow assumption** is exact for channel flow but only an approximation for a growing boundary layer.

(b) With Re → ∞: (U − c)(φ″ − α²φ) − U″φ = 0. Divide by U − c, multiply by φ* and integrate between the walls; integrating φ*φ″ by parts with φ = 0 at the walls gives −∫(|φ′|² + α²|φ|²)dy − ∫U″|φ|²/(U − c) dy = 0. The first integral is real. With 1/(U − c) = (U − c*)/|U − c|², the imaginary part is **c_i ∫U″|φ|²/|U − c|² dy = 0**, so if c_i ≠ 0, U″ must change sign: an inflection point is necessary (Rayleigh's criterion). Plane Poiseuille flow has U″ = −2 everywhere: **inviscidly stable**. Blasius has U″ ∝ f‴ = −½ff″ ≤ 0 (f, f″ ≥ 0): **no inflection, inviscidly stable** too. What destabilises both is **viscosity** itself: the viscous terms near the wall and at the critical layer (where U = c_r) shift the phase between u′ and v′ so that the disturbance extracts energy from the mean shear. These are Tollmien–Schlichting waves, which exist only at finite Re (hence a critical Re and a neutral curve that closes at large Re).

(c) **e⁹ = 8103, e⁷ = 1097.** An initial amplitude of 10⁻⁵U becomes about **0.08U** at N = 9: a velocity fluctuation of 8 % of the free stream, where the linear theory stops (nonlinear terms are no longer small) and secondary instabilities break the wave down to turbulence within a short distance. The e^N method does not model that breakdown; it uses N_crit as a calibrated proxy for "amplitude large enough".

**Common mistakes.** Forgetting the U″φ term (it comes from the v′U′ term after cross-differentiation); claiming Poiseuille flow is linearly stable at every Re (inviscidly, yes; viscously, no: Re_c = 5772).

## 3. A rocket engine's cooling channel (s = 15: h = 2.5 mm)

**Assumptions.** Fully developed plane Poiseuille flow (channel wide compared with its height), water at 20 °C, ν = 1.0038 × 10⁻⁶ m²/s, smooth walls.

(a) U_c = 5772.22 ν/h = **2.32 m/s**; U_b = 2U_c/3 = **1.55 m/s**.

(b) U_b·2h/ν = **7.70 × 10³**; on the hydraulic diameter 4h of a wide channel, **1.54 × 10⁴**.

(c) Transition in channel and pipe flow is **subcritical**: finite-amplitude disturbances (from the inlet, the manifold, wall roughness) trigger turbulence through nonlinear and non-modal growth at Reynolds numbers where every linear eigenmode decays. Linear theory gives the Reynolds number above which *some* infinitesimal wave grows, which is an upper bound for this flow, not a prediction. The e^N method assumes a quiet stream in which slowly amplified waves are the only route; a cooling channel fed by a pump through bends is anything but quiet, and the designer should assume turbulent flow (and size the heat transfer and pressure drop for it, as in `P08.07`).

**Common mistakes.** Mixing the centreline and bulk velocities; using the pipe value of `P08.07` (2300, a different flow and a different definition).

## 4. The neutral curve (s = 15: α₁ = 0.98)

(a) Re_neutral(0.98) = **5934.38** (N = 30), **5926.9966** (N = 40), **5926.9948** (N = 50), **5926.9948** (N = 60). N = 30 has two correct digits, N = 40 six, N = 50 and 60 eight or more: spectral convergence (each ten points gain several digits) until round-off.

(b) Lower branch at α = 0.90, 0.95, 1.00, 1.05, 1.08: **6965, 6208, 5815, 5890, 6478**. The minimum lies between α = 1.00 and 1.05; `critical_point(60)` finds **Re_c = 5772.22 at α_c = 1.02054**, matching the published 5772.22 and 1.02056 to the digits given (Re_c to 3 × 10⁻⁷ relative).

(c) **One** mode grows (c = 0.23752649 + 0.00373967i, Orszag's published eigenvalue). Comparing N = 60 with N = 100, only **7 of the 59** eigenvalues agree to 10⁻⁶; all of them have c_i ≥ −0.183, and the unmatched ones include near-degenerate pairs around c_r ≈ 0.9 with c_i ≈ −0.09. The heavily damped part of the spectrum (and nearly coincident pairs, which are sensitive to round-off) is unreliable at this N; you tell by refining N and keeping only the eigenvalues that do not move. The least-stable mode, the one that matters, converged long before.

**Common mistakes.** Trusting the whole computed spectrum because the top eigenvalue converged; bracketing `neutral_reynolds` where Im c does not change sign.

## 5. What the transition point depends on (s = 15: Tu = 0.20 %)

(a) N_crit = −8.43 − 2.4 ln(0.002) = **6.49**; x_tr = **0.528 m (35 % of the chord)**, against 0.865 m at N_crit = 9.

(b) N_crit = 9, Re_θ = 0.664√Re_x:

| H | Re_θ0 | Re_x,tr | x_tr (m) |
|---|---|---|---|
| 2.55 | 369 | 4.38 × 10⁶ | 1.361 |
| 2.591 | 242 | 2.79 × 10⁶ | 0.865 |
| 2.61 | 206 | 2.32 × 10⁶ | 0.720 |
| 2.65 | 157 | 1.66 × 10⁶ | 0.514 |

Changing H from Blasius's 2.591 to Thwaites's 2.61 (0.7 %) moves x_tr by **−17 %**. The chart is very steep in H because the critical Re_θ depends exponentially on 1/(H − 1): an integral method's small error in H is a large error in transition.

(c) With Thwaites's θ and H for wedge flows and Ue ≈ 60 m/s: **0.195 m (m = −0.05), 0.706 m (m = 0), 1.79 m (m = +0.05)**. A mild pressure gradient moves transition by a factor of 2.5–3.6 either way, comparable to the effect of turbulence between Tu = 0.07 % and 0.3 %. For a real wing, whose pressure gradient is far from zero, the pressure distribution matters at least as much as the atmosphere's turbulence, and it is the one the designer controls.

**Common mistakes.** Using Blasius's k = 0.664 with Thwaites's H (inconsistent; either is defensible if stated); taking Tu in per cent instead of as a fraction in Mack's formula.

## 6. When not to trust e^N

(a) (1) **Bypass transition**: strong free-stream turbulence or the wakes of upstream bodies (a propeller slipstream over the light aeroplane's wing root, a turbine blade row) put large disturbances directly into the layer and skip the slow wave growth. (2) **Roughness**: insect debris on the leading edge, rivet heads, steps and gaps at panel joints trip the layer at the roughness itself (the classic loss of laminar flow on a bug-spattered wing). (3) **Crossflow instability and attachment-line contamination** on swept wings: three-dimensional profiles with inflection points amplify stationary vortices much faster than Tollmien–Schlichting waves, and turbulence from the fuselage can run along the leading edge.

(b) Sweep: the boundary-layer profile on a swept wing has a crossflow component normal to the outer streamline, inflectional and inviscidly unstable; the two-dimensional Falkner–Skan fit knows nothing of it.

(c) A flight measurement of the transition line: infrared thermography of the wing in cruise (laminar and turbulent regions have different surface temperatures because their heat transfer differs), or a chordwise row of hot-film gauges. Transition observed well ahead of 60 % chord on a clean wing in smooth air would contradict the claim; transition that moves forward when the leading edge is lightly roughened would show that the margin is small.

## Check-yourself answers

1. A normal mode is a disturbance of the form φ(y)e^(iα(x − ct)): a wave with wavenumber α travelling at phase speed c_r and growing in time at rate αc_i. It grows if c_i > 0.
2. Re_c = 5772.22 (based on centreline speed and half-height, at α_c = 1.02056; published, Orszag 1971): the lowest Reynolds number at which any two-dimensional infinitesimal wave grows in plane Poiseuille flow. Below it every such wave decays.
3. The eigenfunctions are smooth, so a global polynomial interpolant converges faster than any power of 1/N (spectral accuracy: eight digits at N = 50), and the Chebyshev points cluster near the walls where the viscous layers of the eigenfunction are thin. A low-order finite-difference grid would need far more points.
4. N = ln(A/A₀), the logarithm of the amplification of the most amplified wave since it started to grow. N ≈ 9 (an amplification of 8100) was found by calibrating against observed transition in low-disturbance flight and wind tunnels.
5. The instability point is where waves *start* to grow (Re_x0 ≈ 1.3 × 10⁵ on the e^N chart: 41 mm from the wing's leading edge; the Orr–Sommerfeld critical point of the Blasius profile, Re_δ* ≈ 520, is nearer, at about 28 mm); transition needs them to grow by e⁹, which takes until Re_x ≈ 2.8 × 10⁶ (0.865 m, 58 % chord).
6. Higher free-stream turbulence means larger initial disturbances and therefore a smaller amplification needed: Mack's correlation lowers N_crit from 9.0 at Tu = 0.07 % to 2.6 at 1 %, and moves the wing's x_tr from 58 % to 11 % of the chord.
7. A favourable gradient makes the profile fuller (lower H), which raises the critical Re_θ and lowers the growth rate; the layer is also thinner. Both delay transition (m = +0.05: x_tr 1.79 m against 0.71 m).
8. Bypass transition, roughness-induced transition, crossflow instability or attachment-line contamination (any one).
