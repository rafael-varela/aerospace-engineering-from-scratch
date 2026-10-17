# P09.11  Hints

One hint per problem. Read only the one you need.

1. Put the three terms over the common denominator s(Tf s + 1). At very high frequency the integral term vanishes and Kd s/(Tf s + 1) → Kd/Tf. At low frequency Ki/s is the largest term, so the phase is near −90°.

2. Write u = C_r r − C y and substitute y = P u. The characteristic polynomial is the numerator of 1 + PC, and C contains no b. For (b), divide numerator and denominator of T by C and let s → 0. For (c), num_C = (Kp Tf + Kd)s² + (Kp + Ki Tf)s + Ki: each coefficient is linear in the gains. For (e), set every derivative in the plant equations to zero: α̇ = 0 with q = 0 forces α = 0, and then q̇ = 0 fixes δ.

3. Change `w_servo` in the dictionary from `pitch_plant()` and rebuild the TF. Simulate with `simulate_pitch_loop(p, gains, np.radians(10), 8.0, 1e-3, b=...)` and pass the result to `step_metrics`. Decide with the peak elevator against the ±15° limit as well as the overshoot.

4. φ/δ_a = L_δ/(s(s − L_p)). `pid_place_three_poles(np.array([30.0]), np.array([1.0, -L_p, 0.0]), poles, 0.01)`. The pole you did not place comes from the derivative filter: look near −1/Tf. For U/R, the numerator is num_Cr · den_P and the denominator is the characteristic polynomial.

5. Back-calculation drains the integrator at the rate (u − v)/Tt. Small Tt means fast draining. Call `simulate_pitch_loop(..., u_limit=p['delta_max'], anti_windup='backcalc', Tt=...)` with a 10 s horizon.

6. The gain margin of L = PC is the "ultimate gain" of the loop PC: pass its numerator and denominator to `ultimate_gain`. For the phase margin, sweep ω, find |L| = 1 and read 180° + arg L with the phase unwrapped. Flight-test engineers measure margins without going unstable: think of injecting a small known signal into the closed loop.

7. While saturated, u = δ_max is constant and v = P + I + D. Substitute into İ = Ki e + (u − v)/Tt and collect the terms in I. The frozen-signal test is a first-order lag in I; its fitted time constant should be close to Tt.
