# P05.09  Worked solutions

Numbers are for seed **s = 15** unless stated; rerun the code with your own seed. Convention throughout: forces and distributed loads positive up, couples positive counterclockwise, V = sum of upward forces to the left of the cut, M positive sagging; g = 9.80665 m/s². Every number below was reproduced with `solution.py` (the problem numbers by a probe script that calls its functions). Aircraft data are illustrative, not the data of any type.

## 1. Simply supported equipment beam

**Assumptions.** Rigid supports, a pin and a roller; loads in one plane; self-weight included in the 1.5 kN/m.

**Units.** N, m, N·m.

**Derivation and result.** L = 4.5 m, point load at a = 1.5 m.
(a) ΣM about x = 0: R_B·4.5 = 1500·4.5·2.25 + 4000·1.5, so R_B = **4708.3 N**; ΣF: R_A = 6750 + 4000 − 4708.3 = **6041.7 N**. The sum 10 750 N equals the total load.
(b) V = 6041.7 − 1500x on 0 ≤ x < 1.5 (3791.7 N just left of the load); V = 2041.7 − 1500x on 1.5 < x ≤ 4.5 (−208.3 N just right of it). V never passes through zero smoothly in either segment: it changes sign by jumping across the point load, so **M_max = 6041.7·1.5 − 1500·1.5²/2 = 7375.0 N·m at x = 1.5 m**. (For L ≥ 16/3 m, s ≥ 24, V just right of the load is still positive and the maximum moves into the right segment, at x = (R_A − 4000)/1500.)
(c) `solve_reactions` gives 6041.7 and 4708.3 N; `extreme_moments` gives (1.500 m, 7375.0 N·m).

**Independent check.** Moments of the right free body about x = 1.5 m: 4708.3·3.0 − 1500·3.0²/2 = 14 125 − 6750 = 7375 N·m, the same.

**Common mistakes.** Setting V = 0 in the left segment and getting x = 4.03 m, outside the segment; forgetting that a point load produces a jump in V and a kink in M.

## 2. Derivations

(a) On a slice dx with V, M on the left face and −(V + dV), M + dM (opposite sense) on the right face, vertical equilibrium V + w dx − V − dV = 0 gives **dV/dx = w**; moments about the right face, M + V dx + w dx²/2 − M − dM = 0, give **dM/dx = V** to first order. Across a point force P on [a − ε, a + ε]: ΔV = P + ∫w dx → **P**, ΔM = ∫V dx → **0** as ε → 0. Across a counterclockwise couple: the moment equation gains −M₀, so **ΔM = −M₀** and ΔV = 0.

(b) With k = (w_b − w_a)/(b − a): w_a⟨x − a⟩⁰ + k⟨x − a⟩¹ − w_b⟨x − b⟩⁰ − k⟨x − b⟩¹. For x > b it is w_a + k(x − a) − w_b − k(x − b) = w_a + k(b − a) − w_b = 0. Integrating once (V = 0 at x = 0⁻): w_a⟨x − a⟩¹ + (k/2)⟨x − a⟩² − w_b⟨x − b⟩¹ − (k/2)⟨x − b⟩²; again: (w_a/2)⟨x − a⟩² + (k/6)⟨x − a⟩³ − (w_b/2)⟨x − b⟩² − (k/6)⟨x − b⟩³, the table's row.

(c) Load −w₀x/L (downward of magnitude w₀ at x = L). Resultant w₀L/2 at 2L/3, so R_A = w₀L/6, R_B = w₀L/3. V = w₀L/6 − w₀x²/(2L) = 0 at **x = L/√3**; M = w₀Lx/6 − w₀x³/(6L), and at L/√3, M = (w₀L²/(6√3))(1 − 1/3) = **w₀L²/(9√3)**. Check: L = 3 m, w₀ = 1000 N/m gives 577.35 N·m at 1.7321 m, and `extreme_moments` returns the same.

(d) With x = s sin θ, ∫₀ˢ √(1 − (x/s)²) dx = s∫₀^(π/2) cos²θ dθ = πs/4, so the lift is **πw₀s/4**. With u = 1 − (x/s)², x dx = −(s²/2) du, ∫₀ˢ x√(1 − (x/s)²) dx = (s²/2)∫₀¹ √u du = s²/3, so **M_root = w₀s²/3**. Centre: (s²/3)/(πs/4) = **4s/(3π)**.

**Common mistakes.** Taking moments about the left face (it works, but the w dx²/2 term must still be dropped consistently); a sign error in the couple's jump; putting the triangle's resultant at L/3.

## 3. Wing on trestles

**Assumptions.** Wing as a straight beam, trestles as point supports, structure weight uniform, n = 1, empty tanks.

**Result.** Loads: structure 80g = 784.5 N as 117.1 N/m; engine 150g = 1471.0 N at 1.9 m; total 2255.5 N. With x₂ = 4.75 m:

(a) **R₁ = 1410.9 N** (x = 1.0 m), **R₂ = 844.7 N** (x = 4.75 m), sum 2255.5 N. The inboard trestle carries more because the engine sits 0.9 m from it and 2.85 m from the other.

(b)

| station | V just left (N) | V just right (N) | M (N·m) |
|---|---|---|---|
| trestle 1, 1.0 m | −117.1 | 1293.8 | −58.5 |
| engine, 1.9 m | 1188.4 | −282.6 | 1058.4 |
| trestle 2, 4.75 m | −616.3 | 228.3 | −222.6 |

Largest **sagging 1058.4 N·m at the engine**, largest **hogging −222.6 N·m over the outboard trestle** (the 1.95 m outboard overhang is a cantilever).

(c) Sweeping x₂ from 2.5 to 6.7 m in 1 mm steps: the optimum is **x₂ ≈ 3.41 m**, where the sagging and hogging peaks balance at **634.7 N·m** (reactions 940.1 N and 1315.4 N).

(d) The pull-up root moment, 75 158 N·m, is **71 times** the largest moment on the trestles. Strength in bending is not the concern; the concerns are local: the bearing load of each trestle pad on the skin (the trestle must sit under a rib or spar, not mid-panel), and stability: the centre of gravity (Σ(load·x)/2255.5 = 2.40 m from the root here, between the trestles at 1.0 and 4.75 m) must stay well between them while the engine is removed or work is done.

**Independent check.** Moments of the right part about x = 1.9 m: R₂·2.85 − 117.1·4.8²/2 = 2407.4 − 1349.0 = 1058.4 N·m.

**Common mistakes.** Treating the overhangs as unloaded; refusing the problem because the supports are not at the ends (any two vertical supports make a determinate beam).

## 4. Your own pull-up

n = 3.75, engine at 1.8 m.

(a) **w₀ = 6639.1 N/m**; exact root shear **34 936 N**, root moment **99 343 N·m**; centre of lift **2.8436 m**, independent of n: n scales the lift but not its shape, and the centre is 4s/3π.

(b)

| panels | root M (N·m) | relative error |
|---|---|---|
| 20 | 98 433.4 | 9.16 × 10⁻³ |
| 80 | 99 225.1 | 1.19 × 10⁻³ |
| 320 | 99 328.4 | 1.52 × 10⁻⁴ |

Observed orders 1.47 and 1.49: **order 1.5**, from the square-root tip. The relative errors are the lesson's, because only the shape of the load matters.

(c) Root shear **−23 718 N** (the root reaction pushing down), root moment **74 721 N·m**, **75.2 %** of the lift-alone value. `sections_check` on the exact ellipse gives 74 731.8 N·m: a difference of **10.8 N·m (1.4 × 10⁻⁴)**, the 400-panel error.

(d) V just inboard of the engine −14 137.7 N, just outboard −19 654.0 N; M at the engine **41 074.9 N·m** (continuous). The jump −5516.2 N equals −3.75·150·9.80665 = −5516.2 N.

**Common mistakes.** Measuring the error of the root shear and expecting it to match the moment's; reading the negative root shear as an error (it is the convention: everything to the left of a cut near the root is the root reaction).

## 5. Burn order

Outboard tank on 3.0–5.25 m (centroid 4.125 m), inboard on 0.5–3.0 m (centroid 1.75 m).

(a) Full tanks, 1900 kg: **71 292 N·m**.
(b) 1820 kg: inboard tanks emptied first **69 662 N·m**; outboard first **73 203 N·m**. The difference, 3540 N·m, is the 40 kg of fuel per wing times 3.8 g (1490.6 N) times the difference of the tank centroids (2.375 m).
(c) No wing fuel, 1740 kg: **71 573 N·m**, about 280 N·m more than full. Burning 80 kg per side removes 2981 N of lift acting at 2.844 m (−8477 N·m) but also removes 1490.6 N of relief at 1.75 m and 1490.6 N at 4.125 m (+8757 N·m). The fuel's mean lever arm (2.94 m) is longer than the lift's, so losing the fuel costs more relief than it saves lift.
(d) With empty tanks the root moment rises by 52.98 N·m per kg of aircraft mass. The full-fuel moment is reached at **1734.7 kg**. A maximum zero-fuel weight limits the mass carried in the fuselage (payload) when the wings are dry, because fuselage mass adds lift without adding relief; above it, the dry wing's root moment exceeds the value the spar was sized for.
(e) **Burn the inboard tanks first**: it keeps the relief outboard until the last and lowers the root moment by 3.5 kN·m at the midpoint of the burn. The cost is roll inertia and an outboard fuel system that must keep feeding at the end of the flight, plus a higher moment at the tank boundary (check the whole diagram, not just the root).

**Common mistakes.** Keeping the aircraft mass constant when fuel is burned; forgetting that the lift falls too.

## 6. Judgment

(a) The two-fitting load minus the point load is self-equilibrated, so V and M are unchanged outside 1.6–2.2 m. Inside, M changes by at most (P/2)(d/2) = 5590·0.6/4 = **838 N·m** at the engine station, lowering the sagging moment there: 2.1 % of the 39.7 kN·m. The point load is a fine idealisation for the spar; not for the fittings.

(b) No: the 10-panel root moment is 98 155 N·m against the exact 100 668 N·m, **2.5 % low**. Linear panels through points on a concave curve lie under it, so they lose lift; a coarse model is conservative only when its error has the safe sign, and here it does not.

(c) n = 1.5: root moment = −1.5·9.80665·(80·3.35 + 75·1.75 + 150·1.9) = **−10 065 N·m**, hogging, with a root shear of 4487 N. In hogging the **lower spar cap is in compression** (in the pull-up it is the upper one), so both caps must be checked in compression, the lower one at a much smaller load.

(d) The lift-alone root moment is the lift times the centre-of-lift distance, so it rises by 5 %: **about 5.0 kN·m**, more than the 3.9 kN·m gained by moving the engine 0.7 m outboard. Knowing the lift distribution better matters more than fine-tuning the engine position.

## Check yourself: answers

1. dV/dx = w and dM/dx = V with w positive up. With w positive down (Baker & Haynes), dV/dx = −w and dM/dx = V.
2. Reactions 4 kN each; M_max = wL²/8 = 2·4²/8 = 4 kN·m at mid-span, x = 2 m.
3. Across P, V jumps by P and M is continuous (with a kink). Across a counterclockwise M₀, V is continuous and M jumps by −M₀.
4. 0 at x = 1 and 1 at x = 3; its integral is ⟨x − 2⟩²/2.
5. Just beyond the right end the left free body is the whole beam, which is in equilibrium: the sum of all forces (V) and of all moments (M) vanishes.
6. V(0) = wL (the root reaction, up, for a downward load), M(0) = −wL²/2 (hogging).
7. In flight the engine's inertia load acts opposite to the lift, so its moment about the root subtracts from the lift's. On the ground there is no lift; the engine's weight is the only moment and bends the wing down (hogging).
8. The square-root behaviour of w at the tip: its slope is infinite there, so linear interpolation is less accurate near the tip than for a smooth load, and the error falls as n^(−1.5) rather than n^(−2).
