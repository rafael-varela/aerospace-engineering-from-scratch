# P07.12  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Resistances in series: 1/h, L/k, L/k, 1/h. The largest one sets the hot-face temperature.
2. (b) Sum every row of K T = f; what is Σ_i Σ_j G_ij (T_i − T_j) when G is symmetric? (d) Substitute T_i = ξⁿ e^{ikiΔx} and require |ξ| ≤ 1 for every k.
3. ρc must be an array of node values; give the interface node the average of the two materials. The semi-infinite result needs the heat not to have reached the aluminium.
4. Use `solve_card(kt, n, T_wall, P_chip=...)`; `brentq` on the frame thickness.
5. Fo = αΔt/Δx². Crank–Nicolson's amplification factor tends to −1 for the shortest wavelengths when Fo is large.
6. For a semi-infinite convective wall the hot-face history depends on h√t/√(kρc).
