# P00.14  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Three independent sources, three standard uncertainties, one root-sum-square. The Type A term uses n − 1 in s and then divides by √n; the two Type B terms are half-widths divided by √3.
2. Square the first-order expansion and let every cross term average to zero. For products and powers, divide both sides by y first; the logarithm of a product is a sum.
3. The factor 2 on the velocity term is the whole story. Write each relative uncertainty, square them, and look at which square is largest.
4. `fit.u_slope` is the yardstick for "consistent with 1": a difference under 2 u is unremarkable, over 3 u is real. For (d), read what `polyfit` does with the covariance when `cov=True` and when `cov="unscaled"`.
5. Think about what a scale is physically: a spring or a strain gauge, a reference voltage, a mechanism with friction, a plastic case that creeps. Each one has a signature shape.
6. For (a) the offset enters both readings with the same sign, so in the difference the sensitivity coefficients are +1 and −1 on the *same* quantity. For (b): an uncertainty is a spread around a value you believe; a factor 4.45 is not a spread.
7. The Student's t factor for 95 % and ν = n − 1 degrees of freedom replaces 2 when a budget is Type A dominated. For (c) remember that √h curves downward: the output distribution is skewed.
