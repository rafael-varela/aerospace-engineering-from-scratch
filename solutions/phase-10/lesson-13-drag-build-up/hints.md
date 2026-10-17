# P10.13  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. `cf_mixed(Re, x_tr)`; compare the value at 0.5 with the value at 0.
2. (a) Over the laminar part the turbulent law would have produced friction C_f,t(Re_tr) Re_tr in "Re-weighted" units; replace it by the laminar one. (b) Maximise C_L/(C_D0 + kC_L²). (c) Insert C_L* into V = √(2W/(ρSC_L)).
3. Build each component with `component(name, L_ref, S_wet, FF, Q, V, x_tr)`, then multiply its C_D by 12/11; the wing's exposed area is the wing area minus the fuselage's width times the chord.
4. Pass `gear_dq=` to `cd0_build_up` and recompute e with `oswald_e`, because e depends on C_D0.
5. At 4° the upper surface's suction peak and adverse gradient move transition forward; a NaN means even fully turbulent friction does not reach the measured drag.
6. The section estimate's error applies only to the airframe part of C_D0, and the build-up chose the wing's transition location independently.
