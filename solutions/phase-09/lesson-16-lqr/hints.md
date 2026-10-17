# P09.16  Hints

One hint per problem. Read only the one you need.

1. The scalar CARE is a quadratic in P; take the root that makes a − b²P/r negative. For the limits, write the pole as −√(a² + b²q/r) and let q/r go to 0 and to ∞. "Expensive control" still has to stabilise an unstable mode, and the cheapest way to do so reflects the pole about the imaginary axis.

2. (a) Adjoin λ to the dynamics in the cost and set the derivative with respect to u to zero. (b) Write x = X₁c and λ/2 = X₂c for the same coefficient vector c. (c) Compute J₂HJ₂⁻¹ block by block. (d) Start from −(AᵀP + PA) = Q − PBR⁻¹BᵀP, add and subtract jωP, and multiply on the left by Bᵀ(−jωI − Aᵀ)⁻¹ and on the right by (jωI − A)⁻¹B; recognise RK = BᵀP.

3. A = [[0, 1], [0, L_p]], B = [[0], [L_δ]]. Bryson's rule needs radians, not degrees. The peak aileron for an initial bank error and zero rate is at t = 0.

4. `satellite_model(5 + s)`; the hand formula has J inside k₂ only. Time at the limit: count samples with |u| = u_max in the saturated simulation and multiply by the step.

5. Think about where the loop is broken for the theorem (plant input, single input, full state), and about what an estimator, a delay and a flexible mode each insert into the loop. The delay margin is PM in radians divided by ω_c; compare ωτ at ω_c with the phase margin.

6. (a) Take r → ∞ in the formula from Problem 1: P/r → 2a/b². The closed loop with gain kK is a − 2ka. (b) Call `cost_of_gain` seven times. (c) Change one thing at a time and look at whether J(0.5K) and J(2K) still agree.
