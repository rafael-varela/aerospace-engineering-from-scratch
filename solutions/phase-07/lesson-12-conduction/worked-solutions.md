# P07.12  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15; numbers are from `solution.py` and the probe scripts that call it, given to the figures shown. Material data are Lienhard Tables A.1/A.2 at 20 °C; h, loads and limits are illustrative.

## 1. A regeneratively cooled wall (s = 15: 3.5 mm copper, coolant h = 25 000 W/(m²·K))

**Assumptions.** Steady 1D, constant k, perfect contact between layers.

(a) Resistances (m²·K/W): gas film 1/3000 = 3.333 × 10⁻⁴ (**71.2 %**), Inconel 0.001/11.6 = 8.621 × 10⁻⁵ (**18.4 %**), copper 0.0035/398 = 8.79 × 10⁻⁶ (**1.9 %**), coolant film 1/25 000 = 4.0 × 10⁻⁵ (**8.5 %**). q″ = 2700/4.684 × 10⁻⁴ = **5.765 MW/m²**; hot face **1078.3 K**, Inconel/copper interface **581.3 K**, cold face **530.6 K**.

(b) Copper only: q″ = **7.066 MW/m²**, hot face **644.8 K**. The Inconel layer is the one to remove: it takes 497 K across 1 mm. The copper's resistance is under 2 % of the total, so its thickness hardly matters (which is why liners can be made thick enough for strength). The gas film then dominates: the hot face is set mostly by h_gas.

**Independent check.** The same flux crosses each resistance: (3000 − 1078.3) × 3000 = 5.765 MW/m². ✓

**Common mistakes.** Adding conductances instead of resistances; forgetting the films.

## 2. Derivations

(a) Interior node, V = ΔxΔy: conductances k t Δy/Δx to the east and west neighbours, k t Δx/Δy north and south. Net inflow k t[Δy(T_E − 2T_P + T_W)/Δx + Δx(T_N − 2T_P + T_S)/Δy] = k t V [(T_E − 2T_P + T_W)/Δx² + (T_N − 2T_P + T_S)/Δy²]: k t V times the five-point Laplacian.

(b) Write row i as (K T)_i − f_i = r_i, with r_i = 0 on free rows. K = K_c + D where K_c (conduction) is symmetric with zero row sums and D is diagonal (face and edge-convection conductances). Σ_i (K_c T)_i = Σ_i Σ_j G_ij(T_i − T_j) = 0 because each pair appears twice with opposite signs. So Σ_i r_i = Σ_i D_ii T_i − Σ_i f_i = [face losses + edge convection out] − [generation + flux in + face/edge terms at the fluid temperature], and Σ_i r_i = Σ_fixed r_i. Rearranged: −Σ_fixed r_i (heat out through fixed nodes) + face losses = generation + edge inflows. It uses no property of the solution except that the free rows are satisfied.

(c) Strip of width W, length dx: conduction in, k t W θ′(x), minus out, k t W θ′(x + dx), equals face loss h_f W dx θ: k t θ″ = h_f θ, so θ″ = m²θ with m² = h_f/(k t). General solution A cosh m(L − x) + B sinh m(L − x); θ′(L) = 0 ⇒ B = 0; θ(0) = θ_b ⇒ θ = θ_b cosh m(L − x)/cosh mL; Q = −k t W θ′(0) = k t W m θ_b tanh mL.

(d) Explicit: T_iⁿ⁺¹ = T_iⁿ + Fo(T_{i+1}ⁿ − 2T_iⁿ + T_{i−1}ⁿ). With T_iⁿ = ξⁿ e^{ıκiΔx}: ξ = 1 − 4 Fo sin²(κΔx/2). |ξ| ≤ 1 for all κ requires 1 − 4 Fo ≥ −1, i.e. **Fo ≤ ½**; the shortest wave (κΔx = π) fails first, which is the grid-scale oscillation seen at Fo = 0.55.

## 3. An insulating layer over aluminium

**Assumptions.** 1D, constant properties, perfect bond, adiabatic back; interface node with averaged ρc; Δx = 0.1 mm, Δt = 0.05 s, Crank–Nicolson.

(a)

| t (s) | hot face (K) | aluminium back (K) | √(αt) in silica (mm) |
|---|---|---|---|
| 30 | 710.7 | 337.5 | 5.0 |
| 60 | 883.6 | 453.3 | 7.1 |
| 120 | 1174.1 | 729.7 | 10.1 |

(b) With e = √(kρc) = 1504 J/(m²·K·s^½) for fused silica, T_s = T_i + 2q″√(t/π)/e = **710.9, 881.2, 1121.9 K**. Agreement within 0.3 K at 30 s, 2.5 K at 60 s, and 52 K too low at 120 s: by then the heat has reached the aluminium (√(αt) ≈ 10 mm = the layer thickness), and the aluminium, though a good conductor, has little heat capacity behind 2 mm and an adiabatic back, so it heats up and stops absorbing; the silica's hot face then rises faster than for an infinitely deep layer.

(c) Stored energy **3.000000 × 10⁶, 6.000000 × 10⁶, 1.2000000 × 10⁷ J/m²** at 30, 60, 120 s, equal to q″ t and to `Q_in` to round-off.

**Common mistakes.** Giving ρc per interval instead of per node; forgetting the interface node's mixed capacity; comparing at the aluminium instead of the hot face.

## 4. Your own avionics card (s = 15: 11 W processor, wall 67.5 °C)

(a) Hot spot **131.61, 131.92, 131.83 °C** on n = 41, 81, 161.

(b) Frame **0.330 mm** of 6061 for 100 °C.

(c) The generation uses the exact overlap area of the chip with each control volume, so the total is exactly 16 W on every grid and the global balance closes to round-off (the lesson's card: imbalance ≤ 1.8 × 10⁻¹¹ W). The temperature field is the response to a smoothed source; the footprint edges cause a local, not a global, error, so the peak moves by only 0.3 K across a fourfold refinement.

**Independent check.** Linearity: compared with the lesson (10 W, 70 °C wall, 129.20 °C), the rise here (64.4 K) is larger by about the extra watt's share; the frame needed is almost the same because the target and the rise both moved.

## 5. Time accuracy

(a) Maximum error against the Δt = 10⁻⁴ s Crank–Nicolson reference:

| Δt (s) | implicit (K) | Crank–Nicolson (K) |
|---|---|---|
| 0.04 | 1.221 | 2.283 |
| 0.02 | 0.611 | 0.879 |
| 0.01 | 0.306 | 0.160 |
| 0.005 | 0.153 | 0.0014 |

Orders: implicit **1.00, 1.00, 1.00**; Crank–Nicolson **1.38, 2.46, 6.85**.

(b) Fo = αΔt/Δx² = 1.1575 × 10⁻⁴ × 0.04/(1.25 × 10⁻⁴)² = **296**. At such Fo the shortest grid modes, excited by the sudden start, are damped by implicit Euler (amplification → 0) but not by Crank–Nicolson (amplification → −1): they flip sign each step and decay slowly. CN's error is then dominated by these modes, not by its O(Δt²) truncation error, until Δt is small enough that they die out within the run; that is why its observed "order" jumps around before settling. Implicit Euler has no such modes and shows clean first order.

(c) For the firing time either works at Δt = 0.01 s (errors below 0.3 K against a 500 K rise); Crank–Nicolson is more accurate once the start-up ringing has died out, which it has after a few tenths of a second here. For a sudden step with large Fo and an output needed at early times, prefer implicit Euler (or a few small implicit steps before switching to CN).

## 6. How much should you trust the heat-sink time?

(a) h = 1500: **10.26 s**; 3000: **4.34 s**; 4500: **2.35 s**. Halving h multiplies the time by 2.36, and +50 % divides it by 1.85: weaker than the 1/h² of a semi-infinite wall (whose hot face depends on h√t/e), because the finite wall's back face shortens the long runs more.

(b) Lower k at temperature lowers the effusivity, so the hot face heats faster: 398 W/(m·K) makes the prediction **long** (unconservative). With an effective k around 375 W/(m·K) over the heated region, e² falls by about 6 %, and for a nearly semi-infinite wall the time scales with e², so roughly 6 % shorter (about 0.25 s); Extension 1 computes it properly.

(c) The lesson's sweep: 25 mm gives 4.34 s, 40 mm 5.39 s, 60 mm 5.65 s. A 60 mm wall buys 30 %, not 100 %: within 5 s heat penetrates only about 24 mm, so the extra copper is never used. To double the time, lower h (a lower chamber pressure or a film-cooling ring) or raise the allowable temperature.

(d) The gas-side h (the time is nearly inversely proportional to it), the hot-face temperature limit (and the property it protects), then copper's k(T).

**Common mistakes.** Assuming more thickness always means more time; using a room-temperature k without saying so.

## Check-yourself answers

1. q″ = −k ∇T. k in W/(m·K); α = k/(ρc) in m²/s; k t in W/K (sheet conductance).
2. In series the same heat flows through each layer and the temperature drops add, so resistances (drop per unit flux) add. Side by side each layer sees the same temperature gradient and the heat flows add, so conductances (k t) add.
3. Each row is an exact energy balance of its control volume with symmetric pairwise exchanges, so summing rows cancels all internal conduction exactly; only round-off remains.
4. Where the heated edge meets the cold edges the boundary temperature jumps, the true gradient is singular, the series converges slowly and the numerical error is largest; away from the corners both are accurate.
5. Fo = αΔt/Δx² ≤ ½. The implicit scheme's amplification factor 1/(1 + 4Fo sin²(κΔx/2)) is below 1 for every Fo.
6. Bi = hL/k compares the wall's internal resistance with the film's: small Bi means nearly uniform temperature (lumped). √(αt)/L compares the heat's penetration with the thickness: small means semi-infinite.
7. The hot face's temperature rise under a given h falls with the effusivity √(kρc): copper (36 993) spreads heat away from the surface quickly; Inconel (6605) keeps it at the surface, so its surface reaches the limit in a fraction of a second even though its melting point is higher.
8. Order 2 (log 4/log 2). An order of 1 would suggest a first-order error somewhere: a misplaced source or boundary (for example a point load or edge snapped to the wrong node), not the central differences.
