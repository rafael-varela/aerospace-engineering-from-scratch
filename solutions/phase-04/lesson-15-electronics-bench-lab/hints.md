# P04.15  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. `divider_design(9.0, 3.0, r_total)`, then `nearest_e_series(..., direction="up")` for the top and `"down"` for the bottom; current is V/(R_t + R_b). One LSB at the ADC is 3.3/4096 V; multiply by (R_t + R_b)/R_b.
2. (a) Variance of a uniform distribution of width w is w²/12. (b) Replace the divider by V_out behind R_t ∥ R_b; the meter makes a second divider. (c) Both the divider and the ADC are proportional to V_ref. (d) Differentiate R/(R_lo + R) − R/(R_hi + R) with respect to R.
3. Ratio k = (R_t + R_b)/R_b ≈ V_max/3.0. Worst case: top high and bottom low (and the reverse). RSS: `propagate` on k with u = 1 %/√3 of each resistor. Thevenin: R_t ∥ R_b.
4. `analyse_thermistor(make_synthetic_bench_data(s))`; for the budget write T as a function of (raw, R_f) and use `propagate` with one input's uncertainty at a time, then both.
5. R_b′ = V_r R_t/(V_in − V_r) is the parallel combination; solve R_b′ = R_b R_m/(R_b + R_m) for R_m. The derivative is large because R_m is the small difference between R_b and R_b′.
6. For each symptom ask: which component's error would produce exactly this pattern, and which single measurement isolates it (a resistor unpowered, a voltage across a part, a second instrument, a changed condition)?
