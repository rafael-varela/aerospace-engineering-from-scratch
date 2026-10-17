# P02.08  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Linearity: transform each term separately. For (b), s² + 4s + 13 = (s + 2)² + 9, so the numerator s + 2 is exactly the shifted s of the damped-cosine entry.
2. For (c), the three poles are 0, p, p̄; the residue at 0 is G(0) = 1 and the residue at p is ω_n²/[p(p − p̄)]; write r e^{pt} + r̄ e^{p̄t} = 2 Re(r e^{pt}) and expand. For (d), the derivative of (6) collapses to a single sine term because the cosine terms cancel.
3. P_ind = R/(s²(τs + 1)): a double pole at 0 plus a simple pole at −1/τ. Partial fractions: A/s² + B/s + C/(τs + 1); find A by covering s², then C by covering (τs + 1), then B by comparing the s¹ coefficient (or use SymPy's `apart`).
4. Every number is a function call; use `np.linspace(0, 1, 100001)` for the grid. The pulse is `step_response(...) - delayed_step_response(..., 0.05)`.
5. Rate limiting makes the response depend on amplitude, which is the definition of nonlinear. For (c), overshoot fixes a minimum pole angle (a wedge) and settling time fixes a maximum real part (a vertical line).
6. sin(ω_p t) + sin(ω_p (t − t_p)) u_{t_p}(t) = sin(ω_p t) − sin(ω_p t) u_{t_p}(t) when ω_p t_p = π: the delayed sine is the negative of the original, so the sum is zero after t_p. The transform follows from the shift rule.
