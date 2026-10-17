# P08.03  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
(a) Control volume round the engine and its tanks' feed lines, cut at the mount; momentum outflow ṁV_e through the exit; pressure p_a on all faces except p_e on the exit; F = ṁV_e + (p_e − p_a)A_e. (b) `rocket_thrust` with p_a = 101 325, 22 632, 0; c = F/ṁ; I_sp = c/g₀. (c) Solve p_ISA(h) = p_e; over- and underexpanded.

## 2.
(a) Steady: the volume integral's time derivative vanishes; B = m gives Σρu_nA = 0, B = mu_x gives ΣF_x = Σρu_nu_xA. (b) Inlet ρ₀V₀A_box in; jet ṁ_e out; outer ρ₀V₀(A_box − A_e) out; fuel ṁ_f in; side by difference. (c) Sum ρu_nu_xA + p n_xA over faces; the side face carries u_x = V₀; collect V₀ terms using the mass balance. (d) ΣF_x changes by WΣρu_nA.

## 3.
(a) Control volume from far above to far below: T = ṁw (pressure ambient on its boundary); P = ½ṁw² (energy); at the disc P = Tv_i = ṁwv_i, so v_i = w/2. With ṁ = ρAv_i: T = 2ρAv_i². (b) T = mg/4, v_i = √(T/(2ρA)), P = Tv_i. (c) Same with ρ at 3000 m; ratio √(ρ₀/ρ). (d) Profile and tip losses, swirl, non-uniform inflow, motor and electronics efficiency.

## 4.
(a) `engine_case`. (b) Loop over A_box with `turbojet_control_volume`, `mass_outflow`, `solid_force_x`. (c) V_e,new = √(V_e² + 2c_pT_e[1 − (p_a/p_e)^((γ−1)/γ)]), T_e,new = T_e(p_a/p_e)^((γ−1)/γ), A_e from `exit_area` at p_a; `turbojet_thrust`. Compare with (a). (d) Energy argument: the pressure excess at the exit could have been turned into velocity.

## 5.
(a) `profile_fluxes` for (1 − r/R)^(1/7) and 1 − (r/R)²; error = (β − 1)ṁV̄. (b) For each t_b: ṁ = 60/t_b, F = ṁc, `rocket_burn(100, ṁ, F, t_b, 400)`, `tsiolkovsky_delta_v(c, 100, 40, t_b)`. (c) Gravity loss g t_b; structural and acceleration limits, drag in the dense lower atmosphere.

## 6.
(a) Control volume round the engine, cut at the flexures and the load cell; on the stand V₀ ≈ 0 so ram drag vanishes; pressure thrust remains if the nozzle is not perfectly expanded at ground pressure. (b) Ram drag, the intake and nozzle installation (spillage, interference with the airframe), the different ambient pressure and temperature. (c) A test-cell flow brings in a small ram-drag term and a pressure gradient along the cell. (d) The momentum-flux factor β and the exit pressure term.
