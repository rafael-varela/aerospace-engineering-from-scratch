# P05.15  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. L = 0.4 + 0.02s. 2. `tube_properties(0.0254, 0.00124)`; λ_c from (6). 3. For K = 1 and K = 0.6992: KL/r; compare with λ_c; P_cr = A × `critical_stress`; also `euler_load`. 4. Euler error = P_Euler/P_curve − 1 where Johnson applies.

## 2.
(a) Base fixed: w(0) = 0 → B + D = 0; w'(0) = 0 → Ak + C = 0. Tip free: w''(L) = 0 and E I w'''(L) + P w'(L) = 0; the latter reduces to C = 0 because the trigonometric terms cancel (k² = P/EI). Then A = 0 and B cos kL = 0. (b) Fixed at 0 (B + D = 0, Ak + C = 0), pinned at L (w = 0, w'' = 0): w''(L) = 0 gives A sin kL = −B cos kL; w(L) = 0 then gives tan kL = kL. Solve with `brentq` on (4.4, 4.6). Mode with A = 1: B = −β, C = −k, D = β. (c) F_cy − Cλ_c² = π²E/λ_c² and −2Cλ_c = −2π²E/λ_c³; eliminate C.

## 3.
1. P(d) = A(d) × `critical_stress(E, F_cy, K L / r(d))` with t = 1 mm. 2. `brentq(lambda d: P(d) - 40e3, 0.01, 0.2)`. 3. Compare KL/r with λ_c. 4. d/t. 5. Repeat with K = 1. 6. Mass per metre = ρA.

## 4.
(a) Loop over the end conditions; P₅₀, P₁₀₀; (4P₁₀₀ − P₅₀)/3; relative errors of P₁₀₀ and of the extrapolation against (4). (b) Derive as in 2(b); evaluate on the FD grid; normalise both by their maximum absolute value with the same sign; maximum difference; position of the maximum. (c) n + 5 unknowns; n + 1 equation rows and 4 boundary rows. Find the rows of B that are zero (boundary rows without P) and the columns of B that are zero (the outermost ghosts, which only the fourth difference reaches). Count `~np.isfinite` of `scipy.linalg.eigvals(A, B)` after rebuilding A and B (copy the assembly part of `fd_buckling`).

## 5.
(a) E I (a − w₀)(−π²/L²) sin(πx/L) = −P a sin(πx/L), so (a − w₀)P_E = Pa. (b) σ(P) = (P/A)(1 + a(P) c/r²) with a(P) = w₀/(1 − P/P_E); `brentq` on σ(P) − F_cy for P in (0, P_E). (c) Repeat for the other w₀; compare with the straight-strut Johnson load.

## 6.
Points: λ_c for 4130 is 96.4, so L/r = 60 is in the Johnson range and Euler overstates the strength by tens of per cent whatever K is; "stiffer than pinned" is a statement about K, not about the formula; welding distorts tubes and offsets the load line, adding eccentricity; recommend the column curve with K = 1, an imperfection allowance or a test-based allowable, and ask for the joint drawings or a strut test.
