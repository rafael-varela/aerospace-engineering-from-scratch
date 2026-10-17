# P09.03  Hints

One hint per problem. Read only the one you need, then `method.md` if the hint is not enough.

1. Multiply numerators and denominators for the series; cross-multiply for the parallel; for the loop, the denominator is d₁d₂ + n₁n₂. The DC gain is the ratio of the constant coefficients.

2. (a) Write E = R − HY and Y = GE and eliminate E. (b) In n_G d_H/(d_G d_H + n_G n_H), where does d_H appear? Divide numerator and denominator by n_G and let the gain grow. (c) For the companion A, (sI − A)⁻¹B is (1, s)ᵀ/(s² + a₁s + a₀). (d) The signal at the pick-off point must be the same before and after the move.

3. The inner loop is the wheel and the body with K_d feedback of ω; the outer loop closes K_p around it and the angle integrator. With J and τ_w the denominator is J τ_w s³ + J s² + K_d s + K_p.

4. Loop over K_q with `dict(pitch_model(), K_q=..., K_theta=...)`; the settling time depends on the residues as well as the poles, so the two criteria need not agree.

5. (a) Draw the controller as 1 − 3/(s + 2): a lag block driven by r and a gain. The plant mode e^{t} is still there; what excites it? (b) A root of multiplicity n moves by about ε^{1/n} when a coefficient changes by ε. (c) Divide the elevator limit by the peak elevator per unit command.

6. h/θ = V(1/T_θ2)/(s(s + 1/T_θ2)). T(s) has the zero (s + 1/T_θ2): look at what multiplies what in the characteristic polynomial. For the unreduced simulation, add two blocks to `pitch_blocks` and extend the wiring of `pitch_connect`.
