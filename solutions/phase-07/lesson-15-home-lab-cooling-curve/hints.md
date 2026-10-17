# P07.15  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. τ = ρcV/(hA) is ρc/h times V/A; a longer bar has relatively less end area.
2. For the lag, try T_m = T∞ + a e^(−t/τ) + b e^(−t/τ_p) and fix a from the ODE and b from T_m(0) = T∞. For small Bi, expand J₀ and J₁ for small λ.
3. Without air, h is h_rad alone; separate variables in the nonlinear ODE and integrate dT/(T⁴ − T_w⁴) numerically (or, for a 0 K shroud, exactly).
4. Use the same weighting, window and T∞ as the lesson; a free T∞ and a curved record are a dangerous combination.
5. Bi uses L_c = V/A for the criterion; the water bath has a much larger h and steel a much smaller k.
6. Each candidate cause changes a different part of the record: the level of τ, its drift across windows, or the shape near the end.
