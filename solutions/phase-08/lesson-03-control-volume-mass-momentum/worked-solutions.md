# P08.03  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe script that calls it, given to the figures shown. Engine and vehicle data are illustrative.

## 1. A rocket engine at sea level and in vacuum (s = 15: A_e = 0.065 m²)

**Assumptions.** Steady, uniform exit flow, axial exit velocity, no air enters; the stress on every crossed face is the pressure alone.

(a) Control volume round the engine, cut through the mount. Mass enters only from the propellant lines (with no x-momentum) and leaves through the exit. Momentum: the mount force F_mount plus the pressure forces equals the momentum outflow ṁV_e. The ambient pressure acts on the whole closed surface except the exit, where p_e acts; the net pressure force on the contents is −(p_e − p_a)A_e in the exhaust direction. So the mount must push the engine back (towards the exhaust) with ṁV_e + (p_e − p_a)A_e: the engine pushes the vehicle forward with **F = ṁV_e + (p_e − p_a)A_e**.

(b) Sea level: 12 × 2800 + (40 000 − 101 325) × 0.065 = **29 614 N**; at 11 km: **34 729 N**; vacuum: **36 200 N**. c_sl = **2468 m/s**, I_sp,sl = **251.6 s**; c_vac = **3017 m/s**, I_sp,vac = **307.6 s**.

(c) p_ISA(h) = 40 kPa at **h = 7185 m**. Below it the nozzle is **overexpanded** (p_e < p_a): the pressure thrust is negative, and if p_e is far enough below p_a the boundary layer separates inside the nozzle and the exhaust detaches from the wall, with side loads that can damage it. Upper-stage nozzles are designed for low p_e for exactly this reason and are never fired at sea level.

**Check.** The vacuum and sea-level thrusts differ by p_a A_e = 6586 N (test `test_rocket_thrust_vacuum_minus_sea_level_is_pa_Ae`).

**Common mistakes.** Using absolute rather than gauge pressure for the pressure-thrust term; quoting I_sp in m/s.

## 2. The box around a jet engine

(a) Steady: d/dt ∫ρb dV = 0. With B = m (b = 1): **Σ ρ u_n A = 0** (mass in = mass out). With B = mu_x (b = u_x), and dB_sys/dt = ΣF_x by Newton's second law for the system: **ΣF_x = Σ ρ u_n u_x A**, where ΣF_x is the pressure on the control surface plus every force through solids cut by it.

(b) Engine frame, air arriving at V₀ along +x. Inlet: −ρ₀V₀A_box (in). Jet: +ṁ_e = (1 + f)ṁ_a. Outer exit: +ρ₀V₀(A_box − A_e). Fuel: −fṁ_a. Side by conservation: ṁ_side = ρ₀V₀A_box + fṁ_a − ṁ_e − ρ₀V₀(A_box − A_e) = **ρ₀V₀A_e − ṁ_a**. It is the difference between the free-stream air that would have crossed the jet's exit area and the air the engine actually swallowed. At the lesson's cruise point ρ₀V₀A_e = 8.71 kg/s < ṁ_a = 10 kg/s, so ṁ_side = **−1.29 kg/s**: air flows *in* through the side wall, drawn towards the engine to supply its intake.

(c) Momentum outflow plus pressure terms, face by face: inlet −ρ₀V₀²A_box − (−p_a A_box); jet ṁ_eV_e + p_eA_e; outer exit ρ₀V₀²(A_box − A_e) + p_a(A_box − A_e); side ṁ_sideV₀ (n_x = 0, no pressure term); fuel 0. Sum: ṁ_eV_e + V₀(−ρ₀V₀A_e + ṁ_side) + (p_e − p_a)A_e = ṁ_eV_e + V₀(−ṁ_a) + (p_e − p_a)A_e. So **F = ṁ_eV_e − ṁ_aV₀ + (p_e − p_a)A_e**, and A_box has cancelled. The side face's momentum ṁ_sideV₀ combined with the deficit of the outer exit to leave exactly the ram drag ṁ_aV₀: without the side face the box would give a size-dependent answer.

(d) Every u_x becomes u_x + W, so Σρu_nu_xA gains WΣρu_nA, and the pressures are unchanged. The force is unchanged exactly when **Σρu_nA = 0**: Galilean invariance of the thrust is a consequence of conservation of mass. The test `test_thrust_is_galilean_invariant_because_mass_balances` checks both halves: the full box gives the same force in both frames, and the box with the side face left out (mass residual of order 1 kg/s) gives a force that changes with the frame.

**Common mistakes.** Leaving out the side face; using p n_x on the side wall (its normal has no x-component); counting the fuel's momentum (it enters with zero x-velocity).

## 3. A quadcopter in hover (s = 15: m = 3.0 kg)

**Assumptions.** Actuator disc: uniform, steady, inviscid, no swirl; ambient pressure on the far boundaries; the slipstream contracts to its final speed w far below.

(a) Control volume enclosing the stream tube from far above (air at rest) to far below (speed w): the pressure is ambient everywhere on it, so the only force is the rotor's. Mass: ṁ through every section. Momentum: **T = ṁw − 0**. Energy (no losses, no heat): the power supplied equals the kinetic-energy flux leaving, **P = ½ṁw²**. At the disc the rotor does work at the rate P = Tv_i = ṁwv_i, so ṁwv_i = ½ṁw², **v_i = w/2**: half the final wake speed is gained above the disc, half below. With ṁ = ρAv_i, T = 2ρAv_i².

(b) A = π(0.25)²/4 = 0.04909 m², T = 3.0 × 9.80665/4 = **7.355 N** per rotor. Sea level: v_i = √(T/(2ρA)) = **7.82 m/s**, w = **15.64 m/s**, P = Tv_i = **57.5 W** per rotor, **230.1 W** in total.

(c) At 3000 m: v_i = **9.08 m/s**, P = **66.8 W** per rotor, **267.1 W** total, **16.1 % more** (the ratio √(ρ₀/ρ) = 1.161). The rotor must make the same thrust ṁw with less mass flow, so it must accelerate each kilogram more, and kinetic energy goes as w²: thin air costs power.

(d) Profile drag of the blades (viscous losses); tip losses and non-uniform inflow; swirl in the wake; motor, propeller and speed-controller efficiency. A measured-to-ideal ratio near 2 is typical of small rotors.

**Check.** Dimensions: T/(ρA) has units of (m/s)² ✓. Limit: doubling A at the same T lowers v_i by √2 and the power by √2, the reason helicopters have large rotors.

**Common mistakes.** Writing T = ṁv_i (the momentum change is to the far wake, w); assuming the pressure is ambient at the disc (it jumps there).

## 4. Your own turbojet in cruise (s = 15: ṁ_a = 18 kg/s)

(a) A_e = **0.1824 m²**; gross (1.02 × 18 × 650) = **11 934 N**; ram drag 18 × 236 = **4248 N**; pressure (40 000 − 22 632) × 0.1824 = **3168 N**; net **10 854 N**.

(b)

| A_box (m²) | mass residual (kg/s) | box force − formula (N) |
|---|---|---|
| 0.5 | −5.6 × 10⁻¹⁶ | 0.0 |
| 5 | −5.6 × 10⁻¹⁶ | 1.5 × 10⁻¹¹ |
| 50 | −5.6 × 10⁻¹⁶ | 7.3 × 10⁻¹² |

The side face carries ρ₀V₀A_e − ṁ_a = −2.33 kg/s (inflow).

(c) V_e = √(650² + 2 × 1150 × 900 × [1 − (22 632/40 000)^(0.33/1.33)]) = **833.8 m/s**, T_e = 781.4 K, A_e (from continuity at p_a, with the code's gas constant) = **0.2182 m²**, net thrust **11 061 N**: **1.9 % more**. Holding V_e fixed while setting p_e = p_a would give 7686 N, but that engine is not the same engine: it throws away the exhaust's pressure excess instead of converting it into velocity. The right comparison expands the same gas further.

(d) Thrust per unit propellant is greatest when the expansion extracts all the velocity the exhaust's pressure can give at ambient pressure; the momentum gained by expanding from p_e to p_a, ∫A dp, exceeds the (p_e − p_a)A_e that underexpansion collects at the exit.

**Common mistakes.** Holding A_e fixed when the exit state changes; comparing at the same V_e.

## 5. Non-uniform jets and gravity losses

(a) 1/7th law: β = **1.02041** (quadrature 1.020407 at 200 points, closed form 1.020408); laminar: β = **4/3**. With ṁ_e = 18.36 kg/s and V̄ = 650 m/s, writing ṁV̄ underestimates the gross thrust by **244 N (2.2 % of net)** for the 1/7th law and **3978 N** for the laminar profile. A jet-engine exhaust is turbulent (Re of order 10⁶), so its profile is nearly flat and the 1/7th law is the realistic one; the laminar factor is for slow, viscous pipe flow.

(b)

| t_b (s) | F/W₀ | burnout speed, RK4 (m/s) | c ln(m₀/m_f) − g t_b (m/s) | gravity loss (m/s) | burnout altitude (m) |
|---|---|---|---|---|---|
| 4 | 30.6 | 1793.355 | 1793.355 | 39.2 | 3035 |
| 8 | 15.3 | 1754.128 | 1754.128 | 78.5 | 5912 |
| 16 | 7.6 | 1675.675 | 1675.675 | 156.9 | 11 197 |
| 32 | 3.8 | 1518.769 | 1518.769 | 313.8 | 19 884 |
| 60 | 2.0 | 1244.182 | 1244.182 | 588.4 | 29 045 |

RK4 (400 steps) and the closed form agree to 10⁻⁹ m/s.

(c) Every second of burn spends g × 1 s of velocity holding the rocket up, so the shorter the burn, the smaller the gravity loss. In practice the acceleration (30g for a 4 s burn) is limited by the structure and payload, and a fast rocket low in the atmosphere pays a large aerodynamic drag, which this model ignores; the real optimum burn time is a trade between gravity and drag losses.

**Common mistakes.** Using the burnout mass in place of the initial mass in the thrust-to-weight ratio; forgetting the gravity term in Tsiolkovsky for a vertical burn.

## 6. What does the test stand measure?

(a) A control volume round the engine, cut at the flexures and the load cell. The cut solids carry F = ṁ_eV_e − ṁ_aV₀ + (p_e − p_a)A_e. On the stand V₀ ≈ 0, so the ram drag is zero; the gross thrust and any pressure thrust (non-zero unless the nozzle is perfectly expanded at ground pressure) remain. The load cell reads the static gross thrust, minus the flexures' stiffness times the displacement if they are not zeroed.

(b) Ram drag ṁ_aV₀ in flight; a different ambient pressure and temperature (density changes ṁ_a; p_a changes the pressure thrust); installation effects: intake spillage drag, interference of the exhaust with the airframe, power and bleed-air extraction.

(c) Yes, slightly: air moving through the test cell at 5 m/s brings in a ram-drag term ṁ_a × 5 m/s (90 N for ṁ_a = 18 kg/s) and a static-pressure gradient along the cell that pushes on the engine's projected area. Test cells are calibrated for these effects.

(d) The momentum-flux factor β and the omitted pressure term: ṁV̄ is low by (β − 1)ṁV̄ (about 2 % for a turbulent exhaust) and misses (p_e − p_a)A_e entirely if the nozzle is choked. The load cell measures the sum directly; that is why stands measure force, not flow.

## Check-yourself answers

1. A control volume is a region of space chosen by the analyst, through whose boundary fluid flows; a system is a fixed collection of matter. Newton's laws apply to systems; the Reynolds transport theorem translates them to control volumes.
2. dB_sys/dt = d/dt ∫_CV ρb dV + ∮_CS ρb (u_rel·n) dA: the rate of change of the system's B equals the rate of change inside the control volume plus the net outflow through its surface.
3. F = ṁ_eV_e − ṁ_aV₀ + (p_e − p_a)A_e: gross momentum thrust, ram drag, pressure thrust.
4. Because the momentum balance holds for any control volume; the side face carries the mass the outer faces leave unbalanced, and its momentum cancels the box-size terms exactly.
5. The ratio J A/(ρQ²) of the true momentum flux to the momentum flux of a uniform profile with the same mean velocity; 4/3 for laminar pipe flow, about 1.02 for a turbulent 1/7th-law profile.
6. In a frame moving with the rocket, the propellant leaves with relative speed V_e and the vehicle's momentum changes by the thrust; integrating m dv/dt = F − mg with m = m₀ − ṁt gives Δv = c ln(m₀/m_f) − g t_b.
7. Steady in the vane's frame, the jet arrives at relative speed V − u with mass flow ρA(V − u); the force is ρA(V − u)²(1 − cos θ), the power Fu = ρA(V − u)²u(1 − cos θ) is greatest at u = V/3.
8. Because ram drag ṁ_aV₀ grows with speed while the jet speed V_e is set by the engine; when V₀ approaches V_e the net thrust vanishes. That is why fast aircraft need high exhaust speeds and efficient slow aircraft use propellers and fans that accelerate a large mass of air a little.
