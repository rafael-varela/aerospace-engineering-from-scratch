# P10.18  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. rot = Rotor(mass=1000 + 10s). 2. For each ρ: `hover_momentum(rot.W, rot.A, rho)` and `hover_power(rot, rho=rho)`. 3. Compare P_induced and P_profile between the altitudes. 4. C_T/σ from the `hover_power` output; relate to the mean lift coefficient 6C_T/σ.

## 2.
(a) dL = B ½ρ(Ωr)²c a(θ − λ/r) dr; divide by ρπR²(ΩR)², use r = y/R and σ = Bc/(πR); integrate θ₀r², θ_tw r³ and λr from 0 to 1. (b) ∫λ dC_T = λC_T with λ = √(C_T/2). (c) Expand (r + μ sin ψ)³, average over ψ (odd powers of sin vanish), integrate in r; then argue with (u_T² + u_R²)^{3/2} ≥ |u_T|³. (d) Put μ = 0; then take √(μ² + λ²) ≈ μ and multiply by ΩR.

## 3.
1. For each tip speed: rot = Rotor(tip_speed=...); `hover_power`, `hover_collective`. 2. Induced: C_T^{3/2}(ΩR)³ = (T/ρA)^{3/2}. Profile: ∝ (ΩR)³. 3. μ = 60/ΩR; the retreating blade sees (1 − μ)ΩR at the tip and must carry its share of the thrust at a lower speed, so its cl rises. 4. `best_speeds(rot)`.

## 4.
1. rot = Rotor(f=0.8 + s/50); `best_speeds`, `autorotation`. 2. P_av = 0.9 × 160 kW; rate of climb (P_av − P(V))/W on a fine grid of V; maximum at V_mp. 3. `brentq` on P(V) − P_av for V_max. 4. Check: `forward_flight_power(V_mp, rot, climb=ROC)` should equal P_av. 5. Plot, tangent of slope min(P/V).

## 5.
1. Loop over μ; `profile_power_integral(mu, radial_flow, n=400)`; K = (I − 1)/μ². 2. Locate the reverse-flow region (r < −μ sin ψ). 3. Loop over K with Rotor(K=...); `best_speeds`; P at 60 m/s.

## 6.
(a) Rotor(kappa=...), Rotor(cd0=...); `hover_power`, `autorotation`. (b) μ from `forward_flight_power`; reverse-flow radius r = μ at ψ = 270°; retreating tip (1 − μ)ΩR; advancing tip Mach (1 + μ)ΩR/a. (c) Momentum theory needs a slipstream leaving the disk in one direction; in slow vertical descent the rotor descends into its own wake. (d) Tail-rotor and transmission power; the inflow change and the non-ideal state in descent.
