# P02.09  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. The integrals in (2) run only over (0, dT) where the pulse is 1. Combine a_n² + b_n² with 1 − cos x = 2 sin²(x/2) to get the amplitude.
2. Products of sines and cosines become sums: sin α sin β = ½[cos(α − β) − cos(α + β)]. An integral of cos(kω₀t) over a whole period is zero unless k = 0.
3. The period of |sin| is half that of sin. Integrate |sin(2πf₀t)| cos(nω₀t) over one period of |sin| with ω₀ = 2π(2f₀); the product-to-sum formula leaves integrals of sines with integer multiples that evaluate to simple fractions.
4. Everything is a call to `response_coefficients` and `partial_sum`; peak = max |partial_sum| over one period. The dangerous periods are T = n/f_mode for odd n; a band of ± 0.2 Hz in f_mode is a band in T of width n × 0.2/f_mode² around each.
5. The overshoot is in the *truncated series*, not in the torque: the full series converges to the square wave at every point except the jump. A valve rise time t_r makes the corner smooth on a scale t_r, which affects harmonics with period shorter than about π t_r.
6. For (a), write cos(kω₀t) = Re e^{ikω₀t} and sum the geometric series with ratio e^{2iω₀t} over k = 1, 3, …, N. For (c), (N + 1)ω₀t* = π.
