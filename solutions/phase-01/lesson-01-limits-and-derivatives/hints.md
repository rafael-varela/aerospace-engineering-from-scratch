# P01.01  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Everything follows from u = a₀α/C_Lmax. The slope is a₀ sech²u = a₀(1 − tanh²u), so you only need tanh u once; the "fraction remaining" is sech²u itself.
2. Write both expansions one under the other and subtract or add before dividing. For (c), write out the five-point stencil with step h/2 and compare term by term with [4D(h/2) − D(h)]/3; they have the same four function values with the same weights.
3. M depends on h only through T, and T depends on t only through h: dM/dt = (dM/dT)(dT/dh)(dh/dt). dM/dT = −M/(2T) because a ∝ √T. Above 11 km one of the three factors is zero.
4. Fit the order only where the error is still falling with h. For (c), h_opt ∝ (f/f‴)^(1/3): the lift model at 16° has a much larger third derivative relative to its value than sin does. For (d), there is no subtraction in Im f(x + ih)/h.
5. The noise error of a central difference is at most (noise amplitude)/h, because two independent errors of ±0.005 can differ by 0.01 and you divide by 2h. Truncation grows with h, noise shrinks with h: set them equal.
6. Ask in each case: is the function smooth across the stencil, and is the function value known to more digits than the subtraction throws away? One of the five cases is not a derivative problem at all.
7. If α(C_L) is the inverse of C_L(α) then dα/dC_L = 1/(dC_L/dα) at corresponding points. For bisection, each step halves the bracket: count the halvings from 4° to 10⁻⁶ rad.
