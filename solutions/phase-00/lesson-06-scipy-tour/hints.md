# P00.06  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. `quad` wants a function of one variable; `minimize_scalar` minimises, so hand it the negative thrust; `brentq` needs a bracket where the function changes sign, which is `[τ, t_end]` for a threshold below the peak.
2. In (a) choose `u = t` so that `du` is simple; in (b) divide the tail integral by `a τ²`; in (c) the single-interval error is `−(h³/12) F''(ξ)` and there are `T/h` intervals.
3. Piecewise functions have corners. `quad` may warn or lose accuracy if a corner lies inside its interval; pass the corner times as `points=` or integrate piece by piece.
4. The observed order between `n` and `2n` is `log₂(e_n / e_{2n})`; it should approach 2. Quote digits to where the error analysis says they stop being meaningful, not to where the printout stops.
5. Relative error propagates: `ρ ∝ p` at fixed `T`, and `D ∝ ρ`, so the relative error in `D` equals the relative error in `p`.
6. Match the *kind* of question (root, integral of an ODE, probability, least squares, eigenvalues) to the submodule, not the vocabulary of the question.
7. A constant offset in one reading is the vector `δr = [0, 0.5, 0]`; its effect on `x` is `B⁻¹ δr`, independent of the true load. A fourth, redundant measurement makes `B` tall; consistent readings lie in its column space and inconsistent ones do not.
