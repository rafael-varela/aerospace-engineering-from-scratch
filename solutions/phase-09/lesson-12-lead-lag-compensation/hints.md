# P09.12  Hints

One hint per problem. Read only the one you need.

1. T = 1/(ω_m√α). The zero is at 1/T, the pole at 1/(αT), and the gain at ω_m is 1/√α. In dB, 20 log₁₀. A lead depends on frequency only through ω/ω_m.

2. For (a), set d/dω[arctan(ωT)] = d/dω[arctan(αωT)] and solve for ω². Then use tan(a − b) = (tan a − tan b)/(1 + tan a tan b) and convert tan to sin. For (b), the phase at 10z is arctan 10 − arctan 10β; let β → ∞. For (d), solve |L(jω)| = 1 as a quadratic in ω², then PM = 90° − arctan(ω_c/(2ζω_n)).

3. `p = {'K': 20.0, 'tau': 0.05 + s/1000}`, `ng, dg = servo_tf(p)`, then `design_lead(ng, dg, 30, 55)` and `design_lead_lag(ng, dg, 30, 55, 200)`. Measure with `margins(*series(...))`. For the ramp error, simulate the closed loop with a midpoint-held ramp, as `main()` does.

4. Build the denominator with `np.polymul([1, 0], np.polymul([0.1, 1], [0.02, 1]))`. Phase −180° means arctan(0.1ω) + arctan(0.02ω) = 90°, i.e. 0.1ω × 0.02ω = 1. The lead adds phase at ω_c but also raises the gain at the phase crossover.

5. Noise at the motor command scales with the controller's high-frequency gain. The lag's pole sits at z/β, so its tail decays with time constant of the order of β/z. Sampling with a hold behaves like a delay of dt/2, whose phase cost is ω_c dt/2.

6. Total lead = PM − (180° + arg G(jω_c)). Per stage, use `lead_alpha(total/2)`. Build both stages with `lead_compensator(1.0, alpha, T)`, multiply them with `series`, then choose Kc so that |L(jω_c)| = 1.
