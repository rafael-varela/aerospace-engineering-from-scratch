# P02.01  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Burn time is propellant mass over mass flow; the burnout velocity is Tsiolkovsky's logarithm minus g t_b, and the altitude formula (4) needs m(t_b) = m_dry, not m₀, in the logarithm's prefactor.
2. In (a) substitute u = m₀ − ṁτ so that dτ = −du/ṁ and the integrand becomes ln(m₀/u); integrate by parts with "1 du" as the part you integrate. In (b) the chain rule dh/dt = (dh/dv)(dv/dt) lets you swap the independent variable from t to v. In (c) use ln(1 + x) = x − x²/2 + … and ask where the series converges.
3. Divide by C and move the T term to the left: dT/dt + (G/C)T = (P + GT_a)/C. The integrating factor is e^(Gt/C). The steady state is where dT/dt = 0, and the time constant is the coefficient's reciprocal.
4. Copy `sounding_rocket()` into a dict and overwrite `"cd"`. The burnout index is where `ts` is closest to t_b; the apogee is `argmax` of the altitude column. For the frozen-density bracket you need k = ρ(h_b) C_D A/(2 m_dry).
5. The equilibrium of dv/dt = −g − kv|v| is approached with time constant 1/|f′(v*)| = 1/(2k|v*|). Compare that time with how long the rocket spends falling through one scale height (8.5 km) at its speed. Terminal velocity is only reached when the time constant is short compared with the rate at which the equilibrium itself moves.
6. For (a), write out y_N = y₀ + h Σ f(t_n) and compare with the trapezoid rule, whose error you know from P01.03. For (b), unroll the recursion as a geometric sum and use 1 + hL ≤ e^(hL). For (c), the bound assumes errors can *grow* at rate L; look at the sign of ∂f/∂v.
