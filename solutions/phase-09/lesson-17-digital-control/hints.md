# P09.17  Hints

One hint per problem. Read only the one you need.

1. For one state, the augmented matrix is 2 × 2 and its exponential is e^{−T/τ} and 1 − e^{−T/τ}. For Tustin substitute s = (2/T)(z − 1)/(z + 1) and clear the fraction by multiplying top and bottom by (z + 1). Matched pole-zero: the servo has one pole and one zero at infinity. A pole at a negative real z makes the response alternate in sign every sample.

2. (a) Compute M², M³ for M = [[A, B], [0, 0]] and spot the pattern. (b) Integrate the variation-of-parameters formula over the two sub-intervals separately, and change variables so each integral starts at zero. (c) Put z = e^{jθ} into (2/T)(z − 1)/(z + 1) and simplify with half-angle identities. (d) The hold's impulse response is a pulse of width T; take its Laplace transform and factor out e^{−sT/2}.

3. A biquad (second order over second order) has no zeros at infinity, so matched pole-zero maps both zero pairs and pole pairs. Search |N_d(e^{j2πfT})| on a fine grid of f for its minimum. Plain Tustin moves the notch to (1/(πT)) arctan(πf_bT).

4. Call `sample_rate_study(np.array([10.0 + s]))` and `discrete_loop_margin(1/f, 1/f)` for f in a loop. `hybrid_simulation(T, tau, t_end, n_sub=20)` returns samples every T/20, so take every 20th point.

5. Write the loss with a full-sample delay as 1.5 ω_cT and solve for ω_s = 2π/T. For the anti-aliasing filter, evaluate the phase of 1/((s/ω_a)² + √2 s/ω_a + 1) at ω_c.

6. Bisect on k between a stable and an unstable value; double k until unstable to find the bracket. The critical pole is the eigenvalue of largest modulus at the critical k; its angle divided by T is its frequency. For the continuous loop, scan max Re eig(A − kBK).
