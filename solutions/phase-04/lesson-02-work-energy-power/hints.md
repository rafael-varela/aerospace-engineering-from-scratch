# P04.02  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. h_E = h + V²/2g; compare mgh with ½mV². For (b) the unknown is the top altitude h_t: h_t = h_E − (1.2 V_s(h_t))²/2g. For (c) dh/dt = P_s/(1 + (V/g) dV/dh).
2. (a) Lift is perpendicular to the velocity. (b) dD/dq = 0, or minimise C_D/C_L. (c) Minimise D·V: write it as aV³ + b/V. (d) Energy at constant h_E gives δV in terms of δh; lift at constant C_L gives δL/L = 2δV/V + δρ/ρ; vertical equation m ḧ = δL.
3. Rate of climb = excess power / weight. Stored energy = W Δh = P × time. For (c) `min_power_speed`, `stall_speed`, `power_required_available`; the sink rate is P_r,min/W; dt = dh/sink(h).
4. `np.linspace(hE0, hE1, 161)`; `rutowski_path` returns V*, h*, P_s*. For (c) `best_climb(ac, h)[1]` on a grid of altitudes and `np.trapezoid(1/ps, h)`. Weight enters P_s twice: in D (through C_L) and in the division by W.
5. `specific_excess_power(ac, 5000, 150, n=n)`; bank = arccos(1/n). `brentq(lambda n: specific_excess_power(...), 1, 6)`. For (c) maximise the root over V. For (d) Δh_E = P_s Δt and at constant h, ΔV ≈ gΔh_E/V.
6. Ask in each case which terms in h_E, P_s or (L/D) contain the mass, and what the energy method assumes about the exchange between h and V.
