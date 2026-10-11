# P00.07  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. In level flight $L = W = mg$. Everything else is one formula each, in the order $q \to C_L \to C_D \to D$. Check the units of $q$ before you divide by it.
2. For (a) use the quotient rule and look only at the numerator. For (b) write the quantity to be minimised as $C_D\,C_L^{-3/2}$ and differentiate; the algebra is the same shape as (a) with a different constant in front of $k$.
3. Repeated doubling: multiply by 2, write down the integer part as the next bit, keep the fraction. For the ulp, find the power of two just below the number and subtract 52 from its exponent. `math.ulp` is the check, not the derivation.
4. A known case you can do by hand needs numbers that make the arithmetic trivial, for example $q = 1$ psf and $S = 1$ ft². Think about which test depends on the *scale* of the output and which only on its *shape*.
5. The step is `mdot * dt`; the quantity that matters for drift is whether that step is exactly representable in binary. Fit the slope with `np.polyfit` on `log10` values, skipping any zero drift.
6. Ask of each failure: did I *intend* the output to change? If yes, the reference moves and the commit says why. If no, the size of the difference tells you whether it is rounding (ulps) or a formula (percent).
7. Mutations 1 to 4 are caught; one of them by only a single test. Mutation 5 needs a case where `min(|a|, |b|)` and `max(|a|, |b|)` give different verdicts: choose $a$ and $b$ far apart and a large relative tolerance.
