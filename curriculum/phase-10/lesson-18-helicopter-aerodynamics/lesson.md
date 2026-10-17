# P10.18  Helicopter Aerodynamics

> Hovering is the most expensive thing a helicopter does; flying forward is cheaper until the fuselage's drag catches up, and the bottom of that bucket sets endurance, range and how gently it can autorotate.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the hover results of momentum theory (induced velocity, ideal power, figure of merit) and Glauert's forward-flight inflow equation with its two limits. 2. Derive the blade-element thrust with uniform inflow and linear twist, C_T = (σa/2)(θ₀/3 + θ_tw/4 − λ/2), and the profile power σC_d0/8 (1 + 3μ²), and explain the 4.65 used instead of 3. 3. Implement and test a hover and forward-flight performance model: collective pitch, power parts, Newton solution of the inflow, the power curve, minimum-power and maximum-range speeds, and an autorotation estimate. 4. Verify it by independent routes: blade element against momentum theory, strip integration against the closed form, disk integration of the profile power, and the energy balance of the trimmed rotor. 5. Decide which modelling input (κ, C_d0, f, K) most changes hover power, best speeds and autorotative descent rate. |
| **Prerequisites** | `P10.17` (actuator disk, blade element momentum, figure of merit), `P10.04` (section lift slope), `P03.02` (Newton's method), `P03.04` (numerical integration), `P00.14`. From Phase 08 (recapped in §2.1): `P08.03` (control-volume momentum balance). Skill check: you can integrate (r + μ sin ψ)³ over r ∈ [0, 1] and ψ ∈ [0, 2π]. |
| **Engineering question** | For the reference light helicopter (1100 kg, R = 5.0 m, two blades of 0.25 m chord, tip speed 210 m/s, NACA 0012), how much main-rotor power does it need to hover at sea level, at what airspeed is the power least and the range greatest, and how fast does it descend in steady autorotation? The answers, for the stated κ = 1.15, C_d0 = 0.010 and f = 1.0 m²: 128.4 kW in hover (figure of merit 0.63, collective 7.9° at 0.75R); minimum power 74.9 kW at 23.7 m/s (46 kn); maximum range at 37.9 m/s (74 kn); least autorotative descent rate 6.9 m/s (1366 ft/min), with a glide ratio of 4.4 at 37.9 m/s. |
| **Assumptions and limits** | Uniform (or κ-corrected) induced inflow; Glauert's momentum model in forward flight; constant chord and section C_d0; small angles; level, trimmed flight with the rotor thrust balancing weight and fuselage drag; no tail-rotor, transmission or accessory power; no compressibility (advancing tip Mach 0.81 at 70 m/s), no retreating-blade stall or reverse-flow correction beyond the profile-power factor; sea-level ISA. Momentum theory fails in the vortex-ring and turbulent-wake states of steep descent, so vertical autorotation is not modelled. |
| **Analytical work** | Hover momentum theory and figure of merit; the blade-element thrust and the collective pitch; the induced-power factor; the profile power with and without radial flow; Glauert's inflow and its limits μ → 0 and μ ≫ λ; the trim tilt and the energy balance C_Tλ = C_Pi + C_Pp; the power curve, the tangent construction for maximum range, and the energy estimate of autorotation. |
| **Computational deliverable** | `starter/main.py`: `thrust_coefficient`, `hover_momentum`, `figure_of_merit`, `be_thrust_coefficient`, `be_thrust_strips`, `hover_collective`, `hover_power`, `glauert_inflow`, `profile_power_integral`, `forward_flight_power`, `best_speeds`, `autorotation` (`Rotor`, `power_from_coefficient` and `main` are given). Tests prove: momentum, energy and continuity in hover; λ_h = v_h/ΩR; FM = 1 for an ideal rotor and rising with loading; the blade-element closed form against strips; the collective reproduces C_T and its induced power equals momentum theory; Glauert's inflow satisfies its equation, tends to λ_h at μ → 0 and to C_T/(2μ) at high μ; the profile factor is 3 without radial flow and 4.4–4.9 with it; the energy balance holds to 10⁻¹⁰; the tangent condition at the maximum-range speed; climb adds T·V_c; autorotation is W V_d = P. |
| **Practical activity** | **computational experiment**: momentum against blade element in hover; figure of merit against loading; Glauert's inflow against its limits; the profile-power factor by disk integration; the power curve and its sensitivity to κ, C_d0, f and K; the autorotation estimate. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. No rotor is run; the course never asks for a rotor test at home. |
| **Evidence** | Completed `main.py`, passing tests, `helicopter_power.png`, the derivations of §2, the experiment tables of §4 and a half-page performance note for the reference helicopter stating the hover power, the best speeds, the autorotative descent rate and the input that most affects them. |
| **Assessment** | `assessments/phase-10/grading-rubric.md`; competencies 1–4. Pass: the hover and forward-flight results are derived; all public tests pass; the blade-element/momentum and energy-balance checks are reported; the note ranks the inputs by their effect. |
| **References** | `@ntrs-nasa-cr-3082-rotary-wing-aerodynamics` (Stepniewski, *Rotary-Wing Aerodynamics* Vol. I, https://ntrs.nasa.gov/citations/19790013868) Chapter II: 3.5 Vertical Descent Rates; 3.6 Induced Velocity and Thrust in Nonaxial Translation; 3.7 Power Required in Nonaxial Translation; 3.9 Induced Power in Horizontal Flight; 5.2 Nonuniform Downwash and Tip Losses; Chapter III: 2.8 Rotor Power, and Aerodynamic and Overall Efficiencies in Hover; 3.1 Velocities; 3.4 Blade Profile Drag Contribution to Rotor Power and Drag; 3.8 Propulsive Thrust and Power Required in Horizontal Flight (required-free). `@leishman-intro-aerospace-flight-vehicles` "Helicopters" (https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/helicopters/): Rotor Flow Environment; Hovering Flight Analysis; Axial Climbing & Descending Flight; Forward Flight Analysis; Blade Element Analysis in Forward Flight; Autorotation; Limitations of the Helicopter (required-free, link only). Chapter VI of CR-3082, 9.1 Hover and 9.2 Forward Flight (supplementary-free: how the section's drag and lift feed into rotor performance). Optional textbooks for deeper study: J. G. Leishman, *Principles of Helicopter Aerodynamics* (the full treatment of non-uniform inflow, the induced-velocity curve in descent, profile-power corrections and trim; optional, paid; a library copy serves). |
| **Extensions** | ★ Add a tail rotor (radius 0.8 m, same tip speed) that balances the main-rotor torque at a 6 m tail arm, and recompute hover and the power curve. ★ Replace the uniform-inflow hover blade element by blade element momentum with Prandtl's tip loss (`P10.17`) and compute κ instead of assuming it. ★ Add altitude: repeat the power curve at 2000 m and 3000 m and find the hover ceiling for an available power of 160 kW. |
| **Time** | unmeasured |

## 1. The engineering problem

A helicopter's rotor does three jobs: it holds the aircraft up, it tilts to pull it forward, and in an engine failure it windmills to let it down gently. All three are questions of power. The reference light helicopter has a gross mass of 1100 kg, one two-bladed main rotor of 5.0 m radius with 0.25 m chord blades of NACA 0012 section, turning at 401 rpm for a tip speed of 210 m/s. Its engine and transmission set the power available. To know whether it can hover out of ground effect, how fast to fly for the longest endurance or the longest range, and how quickly it descends with the engine stopped, you need the power required against airspeed.

The model hierarchy is the one of `P10.17`. **Momentum theory** gives the least possible hover power for this disk. **Blade element theory** with uniform inflow adds the blades: their pitch, their drag, and the induced-power factor κ that says how far real inflow is from uniform. In forward flight, **Glauert's inflow model** extends momentum theory to an edgewise disk, and the **power curve** adds three parts: induced, profile and parasite. Each level is checked against its neighbour.

Success criteria, written first: blade-element thrust with uniform inflow matches the strip integration to 10⁻⁵ and its induced power equals momentum theory; Glauert's inflow tends to the hover value as μ → 0 and to C_T/(2μ) at high speed; the profile-power factor computed by integrating over the disk is 3 without radial flow and about 4.6 with it; the trimmed rotor's power satisfies the energy balance to 10⁻¹⁰; the maximum-range speed satisfies the tangent condition dP/dV = P/V.

## 2. Understand and derive

### 2.1 Recap: momentum theory in hover (`P08.03`, `P10.17`)

For an actuator disk of area A in still air producing thrust T, mass, momentum and energy conservation give a far-wake velocity twice the disk's induced velocity v_h, and

**v_h = √(T/(2ρA)),  P_ideal = T v_h = T^{3/2}/√(2ρA).**   (1)

Rotor coefficients use the tip speed ΩR: C_T = T/(ρA(ΩR)²), C_P = P/(ρA(ΩR)³), λ = inflow/ΩR, so λ_h = √(C_T/2) and C_P,ideal = C_T^{3/2}/√2. The **figure of merit** compares the ideal with the actual power, **FM = P_ideal/P**.

For the reference rotor: A = 78.54 m², disk loading T/A = 137.4 N/m², **v_h = 7.49 m/s**, λ_h = 0.03566, C_T = 0.002543, and the **ideal hover power is 80.81 kW**.

### 2.2 Blade element in hover

At radius r (nondimensional r = y/R), the blade section sees the tangential speed Ωr and the inflow λΩR. With the small-angle inflow angle φ = λ/r, pitch θ(r) and lift slope a, the local lift coefficient is a(θ − λ/r). Summing B blades of chord c over the disk, and writing the **solidity** σ = Bc/(πR) (blade area over disk area, 0.03183 here):

dC_T = (σa/2)(θr² − λr) dr.   (2)

With **linear twist** θ = θ₀ + θ_tw r and uniform λ,

**C_T = (σa/2)(θ₀/3 + θ_tw/4 − λ/2).**   (3)

(The θ/3 says that thrust comes mostly from the outer blade, where the dynamic pressure grows as r².) In hover λ = λ_h from (1), so (3) gives the collective. For the reference rotor (a = 5.73 /rad, θ_tw = −8°, both stated assumptions): **θ₀ = 13.86°, θ_0.75 = 7.86°**; strip integration of (2) gives the same C_T, 0.0025433. The mean lift coefficient of the blade is 6C_T/σ = 0.479, comfortably below stall.

**Induced power.** The blade element's induced power is ∫λ dC_T = λ C_T: with uniform inflow from (1), exactly C_T^{3/2}/√2, the momentum result. Real inflow is not uniform and the tip loses lift (`P10.17`), so the induced power is κ times the ideal with **κ ≈ 1.15** for a well-twisted rotor (a stated assumption; CR-3082 Chapter II 5.2 discusses its origin).

**Profile power.** Each section's drag ½ρ(Ωr)²c C_d0 at speed Ωr absorbs ½ρc C_d0 Ω³r³; summing over B blades and r,

**C_P0 = σC_d0/8.**   (4)

With C_d0 = 0.010 (assumed: the smooth NACA 0012 section has cd ≈ 0.006–0.008 near zero lift at these Reynolds numbers, Ladson's data in `curriculum/phase-10/data/`, and real blades with erosion and roughness more), the **hover power** is

C_P = κ C_T^{3/2}/√2 + σC_d0/8:  **92.93 kW induced + 35.45 kW profile = 128.38 kW, FM = 0.629.**   (5)

### 2.3 Forward flight: Glauert's inflow

In forward flight at speed V the disk is nearly edgewise; it tilts forward by α (the tip-path-plane angle) to supply a propulsive force. Glauert assumed the momentum relation with the total velocity through the disk in place of the axial one:

T = 2ρA v √((V cos α)² + (V sin α + v)²).

In rotor coefficients, with the **advance ratio** μ = V cos α/(ΩR) and λ = (V sin α + v)/(ΩR),

**λ = μ tan α + C_T/(2√(μ² + λ²)),  λ_i = C_T/(2√(μ² + λ²)).**   (6)

*Limits.* At μ = 0, λ_i = λ_h. At μ ≫ λ, λ_i → C_T/(2μ), i.e. v = T/(2ρAV): the rotor acts like a wing of span 2R deflecting the air that passes through a circle of its diameter, which is why the induced power falls rapidly with speed. Equation (6) is a quartic in λ; Newton's method from λ_h converges in four iterations at every μ tested.

**Trim to first order.** The fuselage's drag D = ½ρV²f (f is the **equivalent flat-plate area**, assumed 1.0 m², typical of a light helicopter) must be balanced by the forward tilt of the thrust: tan α = D/W, T = √(W² + D²).

### 2.4 The power curve

The rotor's shaft power is T(V sin α + v) + profile power, i.e. C_P = C_T λ + C_P0. With κ for the non-ideal inflow,

**C_P = κ C_T λ_i + (σC_d0/8)(1 + Kμ²) + DV/(ρA(ΩR)³).**   (7)

*Energy balance (the check).* With κ = 1 the rotor's induced-plus-propulsive power is T(V sin α + v) = Tv + TV sin α. Trim gives T sin α = D, so TV sin α = DV: the propulsive part of the rotor's work is exactly the fuselage's drag times speed, and C_T λ = C_Pi + C_Pp. The tests check it to 10⁻¹⁰.

*The profile factor K.* A blade element at radius r and azimuth ψ sees u_T = r + μ sin ψ (tangential) and u_R = μ cos ψ (radial), in units of ΩR. If drag acts only on the tangential component, power ∝ u_T³, and averaging over the disk,

(1/2π)∫∫(r + μ sin ψ)³ dr dψ = 1/4 + (3/4)μ²  ⇒  **K = 3.**   (8)

If the section's drag acts on the whole in-plane velocity, the power is ∝ (u_T² + u_R²)^{3/2}, and the integral (computed by `profile_power_integral`) gives K = 4.53 at μ = 0.1, 4.60 at 0.2, 4.71 at 0.3 and 4.84 at 0.4. **K = 4.65** is a commonly used constant that represents this extra drag of the radial flow over the useful speed range; this lesson uses it, and Experiment 4 shows that the choice between 3 and 4.65 moves the maximum-range power by 1.3 %.

### 2.5 Best speeds and autorotation

The **minimum-power speed** V_mp minimises (7): best endurance and best rate of climb, since the excess power (P_available − P) is largest there. The **maximum-range speed** minimises P/V, the energy per metre: it is where a line from the origin is tangent to the power curve (dP/dV = P/V).

In steady **autorotation** the engine delivers nothing; the descent supplies the power, W V_d = P_required. The least descent rate is therefore P_min/W at V_mp, and the flattest glide is at the maximum-range speed. This energy estimate neglects the change of inflow in descent and the transmission's own losses (the rotor drives no engine but still drives the tail rotor and accessories), so real descent rates are higher.

### 2.6 Where the model stops

- **Descent through the rotor's own wake** (vortex-ring and turbulent-wake states) breaks momentum theory; vertical autorotation is therefore not modelled.
- **High μ.** The retreating blade sees (r − μ)ΩR at ψ = 270°: at 70 m/s (μ = 0.32), the inner 32 % of the retreating blade is in reverse flow and the retreating tip moves at only 143 m/s, so it must work at high cl and stalls; the advancing tip reaches Mach 0.82. Neither is modelled.
- **κ, C_d0 and f** are assumptions; Experiment 5 ranks them.

## 3. Build the computational model

Open `starter/main.py`. `Rotor` (with properties A, σ, W, Ω), `power_from_coefficient` and `main` are given.

### Step 1: hover

`thrust_coefficient`, `hover_momentum` (1) with far-wake velocity 2v_h and wake area A/2, `figure_of_merit`, `be_thrust_coefficient` (3), `be_thrust_strips` (midpoint rule on (2)), `hover_collective` (solve (3) for θ₀ with λ = λ_h) and `hover_power` (5).

### Step 2: Glauert's inflow by Newton's method

```python
lam = np.sqrt(CT / 2) + mu * np.tan(alpha)
for it in range(1, max_iter + 1):
    s = np.sqrt(mu**2 + lam**2)
    f = lam - mu * np.tan(alpha) - CT / (2 * s)
    fp = 1 + CT * lam / (2 * s**3)
    step = f / fp
    lam -= step
    if abs(step) < tol: break
```

f′ > 1 everywhere, so the iteration is well behaved.

### Step 3: the profile factor and the power curve

`profile_power_integral(mu, radial_flow)` averages |u_T|³ or (u_T² + u_R²)^{3/2} over a midpoint grid in r and ψ and divides by 1/4. `forward_flight_power(V, rotor, climb)` trims (tan α = D/W, T = √(W² + D²)), solves (6), and returns the parts of (7) (plus C_T V_c/ΩR in a climb). `best_speeds` uses `minimize_scalar` on P and on P/V; `autorotation` divides by W.

Running `main()` prints (abridged):

```
rotor: W = 10791.0 N, A = 78.540 m^2, disk loading 137.4 N/m^2, sigma = 0.03183, Omega = 42.0 rad/s (401 rpm)
momentum hover: v_h = 7.489 m/s, lambda_h = 0.03566, ideal power 80.81 kW, C_T = 0.002543, C_T/sigma = 0.0799
hover power: induced 92.93 kW (kappa = 1.15), profile 35.45 kW (sigma Cd0/8, Cd0 = 0.01), total 128.38 kW, FM = 0.629
blade element, uniform inflow: theta0 = 13.858 deg, theta_75 = 7.858 deg (twist -8 deg); C_T closed form 0.0025433, strips 0.0025433
minimum power 74.86 kW at V = 23.66 m/s (46 kn); maximum range at V = 37.88 m/s (74 kn), P = 92.49 kW
  V =    0 m/s: mu = 0.0000, induced  92.93 kW, profile  35.45 kW, parasite   0.00 kW, total 128.38 kW
  V =   20 m/s: mu = 0.0952, induced  34.38 kW, profile  36.95 kW, parasite   4.90 kW, total  76.23 kW
  V =   40 m/s: mu = 0.1897, induced  17.48 kW, profile  41.38 kW, parasite  39.20 kW, total  98.06 kW
  V =   60 m/s: mu = 0.2799, induced  12.04 kW, profile  48.37 kW, parasite 132.30 kW, total 192.71 kW
  V =   70 m/s: mu = 0.3211, induced  10.68 kW, profile  52.45 kW, parasite 210.09 kW, total 273.22 kW
autorotation (energy balance): least descent rate 6.94 m/s (1366 ft/min) at 23.7 m/s; best glide 37.9 m/s, descent 8.57 m/s, glide ratio 4.42
energy balance at 50 m/s (kappa = 1): C_T lambda = 9.9737e-05, C_Pi + C_Pp = 9.9737e-05
```

(Speeds in knots and descent rates in ft/min are given because helicopter flight manuals use them.) It also writes `helicopter_power.png`: the power curve with its three parts and the tangent for maximum range, and λ_i/λ_h against μ with the high-speed limit.

### Verify

Run `python -m pytest curriculum/phase-10/lesson-18-helicopter-aerodynamics`. The tests cover momentum, energy and continuity in hover; λ_h = v_h/ΩR; FM = 1 for an ideal rotor, below 1 and rising with loading otherwise; (3) against strips for three twist and inflow cases, its linearity in θ₀ and fall with λ; the collective reproducing C_T and the blade-element induced power equalling C_T^{3/2}/√2 to 10⁻⁵; the hover power parts; (6) satisfied to 10⁻¹³ in under 30 Newton steps, tending to λ_h and to C_T/(2μ) with a shrinking error; K between 2.95 and 3.1 without radial flow (the small excess is the reverse-flow region, where |u_T|³ ≠ u_T³) and between 4.4 and 4.9 with it; the zero-speed power equal to hover; the energy balance and the parasite power DV; the bucket and the tangent condition; heavier is costlier and a climb adds TV_c exactly; autorotation's W V_d = P.

### Compare

Blade element with uniform inflow and momentum theory agree exactly, as §2.2 says they must: induced power λC_T = C_T^{3/2}/√2. Glauert's inflow (6) agrees with its high-speed limit to 0.8 % at μ = 0.1, 0.05 % at μ = 0.2 and 0.01 % at μ = 0.3. The analytic K = 3 is reproduced by disk integration (3.002 at μ = 0.1). The energy balance closes to 10⁻¹⁰.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** No rotor is run. Write each prediction before running it.

1. **Figure of merit against loading.** Hover the same rotor at 700, 900, 1100 and 1300 kg. Predict the trend of FM. *Finding:* C_T/σ = 0.0508, 0.0654, 0.0799, 0.0944; P = 82.6, 104.2, 128.4, 154.8 kW; **FM = 0.496, 0.574, 0.629, 0.670**. The profile power is fixed while the induced power grows as T^{3/2}, so a more heavily loaded rotor is a more "efficient" one by this measure: FM compares rotors fairly only at the same disk loading.

2. **Glauert's inflow.** λ_i/λ_h at μ = 0–0.3 against C_T/(2μ). Predict the speed at which the limit is good to 1 %. *Finding:* λ_i/λ_h = 1.000, 0.925, 0.648, 0.354, 0.178, 0.119 at μ = 0, 0.02, 0.05, 0.1, 0.2, 0.3; the high-speed limit is within 1 % by μ = 0.1 (V ≈ 21 m/s). Below μ ≈ 0.05 the two differ greatly: the transition from hover is where induced power collapses.

3. **The power curve.** P(V) and its parts. Predict V_mp. *Finding:* 128.4 kW in hover, 76.2 kW at 20 m/s, **74.9 kW minimum at 23.7 m/s**, 98.1 kW at 40 m/s, 192.7 kW at 60 m/s, 273.2 kW at 70 m/s. At V_mp the induced power has fallen by two thirds and the parasite power (∝ V³) has barely started; above 40 m/s the parasite power dominates.

4. **Profile factor.** K by integration, and the speeds with K = 3, 4.65, 5. Predict the sensitivity. *Finding:* K = 4.51 at μ = 0.05 rising to 4.84 at μ = 0.4 with radial flow. With K = 3, 4.65, 5: V_mp = 24.0, 23.7, 23.6 m/s; P_min = 74.1, 74.9, 75.0 kW; V_max range = 38.2, 37.9, 37.8 m/s; P at 60 m/s = 188.1, 192.7, 193.7 kW. The factor matters at high speed (2.4 % at 60 m/s) and hardly at the bucket.

5. **Which input matters most?** κ = 1.10/1.15/1.25 and C_d0 = 0.008/0.010/0.012: hover power, FM and least descent rate. Predict the ranking. *Finding:* κ: 124.3/128.4/136.5 kW, FM 0.650/0.629/0.592, V_d 6.82/6.94/7.17 m/s. C_d0: 121.3/128.4/135.5 kW, FM 0.666/0.629/0.596, V_d 6.24/6.94/7.63 m/s. A ±20 % change of C_d0 moves hover power by ±5.5 % and the autorotative descent rate by ±10 %: the profile drag matters more in autorotation (at V_mp the profile power is half the total) than in hover. In forward flight, f dominates above 40 m/s (Problem 4).

Raw data format: CSV tables `mass_kg, CT_over_sigma, P_kW, FM` (1); `mu, lam_i_over_lam_h, CT_over_2mu_over_lam_h` (2); `V_ms, P_induced_kW, P_profile_kW, P_parasite_kW, P_kW` (3); `mu, K_noradial, K_radial` and `K, V_mp, P_min_kW, V_mr, P_60_kW` (4); `parameter, value, P_hover_kW, FM, Vd_ms` (5).

## 5. Validate and improve

| quantity | route 1 | route 2 | difference | cause |
|---|---|---|---|---|
| hover C_T at the collective | closed form (3): 0.0025433 | strips (200): 0.0025433 | < 10⁻⁷ | midpoint rule on a quadratic |
| induced power, uniform inflow | blade element λC_T | momentum C_T^{3/2}/√2 | 10⁻⁵ (2000 strips) | none: the same physics |
| λ_i at μ = 0.3 | Glauert (6): 0.00424 | C_T/(2μ): 0.00424 | 0.01 % | λ ≪ μ |
| profile factor | analytic 3 | integration 3.002 at μ = 0.1 | 0.002 | reverse-flow region |
| energy at 50 m/s, κ = 1 | C_Tλ = 9.9737 × 10⁻⁵ | C_Pi + C_Pp = 9.9737 × 10⁻⁵ | < 10⁻¹⁰ | none |

The checks verify the code, not the inputs. The engineering answer: **hover 128.4 kW (FM 0.63); least power 74.9 kW at 23.7 m/s (46 kn); best range at 37.9 m/s (74 kn); least autorotative descent about 6.9 m/s (1370 ft/min), glide ratio 4.4.** Treat the descent rate as a lower bound: the energy estimate leaves out the tail rotor and transmission, which the rotor still drives in autorotation, and the descent's effect on the inflow; a flight-manual value for a real helicopter of this size is the right comparison, and it will be higher. To improve: replace κ by a blade-element-momentum hover calculation with tip loss (Extension 2), add the tail rotor (Extension 1), and use measured C_d0 for the actual blade (Ladson's 0012 data near zero lift give 0.0062–0.0066 with free transition and 0.0084–0.0110 with a wraparound leading-edge grit strip).

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (1), (3)–(8); `helicopter_power.png`; the five experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) giving hover power, the best speeds and the autorotative descent rate, with the input that most affects each.

### Adversarial self-checks

- **What if the mass rises?** Hover power grows faster than mass (T^{3/2} for the induced part): 1300 kg needs 154.8 kW, 21 % more for 18 % more mass. Altitude does the same through ρ (Problem 1).
- **What if V → 0?** (7) returns exactly the hover power (a test); the forward-flight model is continuous with hover, but the real transition region (μ < 0.05) has strongly non-uniform inflow that (6) does not capture.
- **Which assumption fails first in the real world?** In hover, κ and C_d0 (Experiment 5); at high speed, retreating-blade stall and advancing-tip compressibility (§2.6); in autorotation, the uniform-inflow picture itself.
- **What measurement could falsify the model?** Flight-test power against speed (engine torque × rotor speed) at a known weight: the bucket's speed and depth test κ, C_d0 and f separately because they dominate different parts of the curve.
- **Does the result survive a change of method or grid?** The strips, the Newton tolerance and the disk-integration grid change the numbers by less than 10⁻⁵; the modelling inputs change them by 5–10 %.
- **Which decision changes if the uncertainty doubles?** A hover-ceiling or payload decision: if the power required were uncertain by ±11 % instead of ±5.5 %, a margin of 10 % between available and required hover power would no longer guarantee hover out of ground effect.

## Check yourself

1. Define C_T, C_P, λ, μ and σ for a rotor. What is λ in hover according to momentum theory?
2. What is the figure of merit, and why does it rise with disk loading for a given rotor?
3. Why is the blade-element C_T proportional to θ/3 rather than θ/2 for an untwisted blade?
4. State Glauert's inflow equation and its two limits.
5. Name the three parts of the forward-flight power and the speed range in which each dominates.
6. Why is the profile-power factor 3 when only the tangential velocity is counted, and why do engineers use about 4.65?
7. How do you find the maximum-range speed on a power curve, and why is it faster than the minimum-power speed?
8. Why does the energy estimate of the autorotative descent rate come out lower than flight-manual values?
