# P03.11  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact, which are numerical, and what the numerical error is.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: the numbers that size a run)

An aluminium wing skin (k = 167 W/(m·K), ρ = 2700 kg/m³, c_p = 896 J/(kg·K)) is L = 3 + s/10 mm thick and is to be marched for 10 s of a transient.

(a) Give α, the diffusion time L²/α, and the number of diffusion times in the 10 s run.
(b) On n = 50 cells, give the FTCS stability limit Δt_max and the number of FTCS steps the run needs. Repeat for n = 100 and state the scaling law.
(c) At Δt = 1 ms with BTCS or Crank–Nicolson on n = 50 cells, give r and the number of steps. Which of the two would you choose, and what would you check before trusting it?

Answer form: six numbers with units, one scaling law, one sentence.

## 2. (Derivation: truncation error and amplification factors)

(a) Starting from the Taylor expansions of u(x ± Δx, t) and u(x, t + Δt), derive the truncation error of FTCS, equation (4) of the lesson, and show that for a solution of the heat equation the leading terms cancel at r = 1/6. Then derive the truncation error of BTCS by expanding about t^{m+1} and explain why no value of r makes it cancel.
(b) Derive the von Neumann amplification factors of FTCS, BTCS and Crank–Nicolson. Prove that FTCS is stable iff r ≤ ½, that BTCS satisfies 0 < g ≤ 1 for all r, and that Crank–Nicolson satisfies |g| ≤ 1 for all r.
(c) Show that for the exact factor e^{−z}, with z = 4r sin²(θ/2), the three schemes give 1 − z (below), 1/(1 + z) (above) and (1 − z/2)/(1 + z/2) (error of order z³). Which way does each scheme bias the decay of a mode, and what does that do to the time at which a transient reaches a threshold?

Answer form: three derivations with every expansion written out; one paragraph for (c).

## 3. (Unfamiliar scenario: an insulated face by a ghost node)

A liquid-hydrogen tank wall is covered by 25 mm of closed-cell foam (k = 0.03 W/(m·K), ρ = 40 kg/m³, c_p = 1500 J/(kg·K)). At t = 0 the tank is filled: the inner face of the foam drops to 20 K and stays there, the outer face is (for this problem) perfectly insulated, and the foam starts uniform at 300 K.

(a) Give α and the diffusion time of the foam. Is one hour long compared with it?
(b) The insulated face T_x(L) = 0 is imposed with a **ghost node** u_{n+1} = u_{n−1}. Derive the FTCS update for the boundary node n and show that the scheme's stability limit is unchanged at r ≤ ½ (hint: the boundary row of the update matrix has entries 2r and 1 − 2r; use the Gershgorin circle theorem, or the discrete cosine eigenvectors of the Neumann problem).
(c) Modify `step_ftcs` (in a copy, keeping the tested one intact) to handle this boundary and march to t = 3600 s on n = 25, 50 and 100 cells at r = 0.4. Report the mid-plane and outer-face temperatures and show they have converged.
(d) Check against the mixed-boundary series of `P02.11`: T = T_s + (T₀ − T_s) Σ_m 4/((2m + 1)π) sin(λ_m x) e^{−αλ_m² t}, λ_m = (2m + 1)π/(2L).

Answer form: two numbers for (a), a derivation for (b), a 3-row table for (c), one comparison for (d).

## 4. (Numerical, uses your code: the bulkhead with a hotter duct)

Repeat the lesson's bulkhead case with T_hot = 600 + 10s K (all else as in `titanium_bulkhead`), n = 100 cells.

(a) Give the steady mid-plane temperature and confirm the 450 K limit is reached. Find the exact time from `time_to_limit_exact`.
(b) Run `time_to_limit` with FTCS at r = 0.45, Crank–Nicolson at Δt = 0.1 s and BTCS at Δt = 0.1 s. Tabulate time, error against (a), and steps.
(c) Give the cold-face heat flux at the moment the limit is reached, and the steady flux. What fraction of the steady coolant load has arrived by then?
(d) Predict, before running, how the BTCS error changes at Δt = 0.05 s; then run it.

Answer form: a 3-row table, three numbers with units, one prediction with its check.

## 5. (Numerical, uses your code: buying accuracy with order)

For mode 1 on the bulkhead's L and α, with t_f = L²/(απ²), run `temporal_convergence` for BTCS and Crank–Nicolson on n = 200 cells with Δt = t_f/5, t_f/10, t_f/20, t_f/40.

(a) Tabulate the errors and the observed orders.
(b) Using the observed first-order law, estimate the Δt at which BTCS would match Crank–Nicolson's error at Δt = t_f/5, and the ratio of step counts. Check by running it.
(c) Run `spatial_convergence` for FTCS at r = 0.4 and r = 1/6 on n = 10, 20, 40, 80. By what factor is the r = 1/6 error smaller at n = 80, and what order do you observe? Explain from equation (4).
(d) The r = 1/6 result depends on three things being true. Name them, and say which one fails first in the bulkhead problem with temperature-dependent conductivity.

Answer form: two tables, one estimate with its check, one short paragraph.

## 6. (Judgment: a colleague's time step)

A colleague proposes to march the bulkhead step problem with Crank–Nicolson at Δt = 1 s on n = 100 cells, arguing that "it is unconditionally stable and second order, so the big step is free".

(a) Compute r and g(π). How many steps would it take for the grid-scale mode to decay to 1 % of its initial amplitude? Compare with the number of steps in the first 5 s.
(b) Run it to t = 5 s and report the maximum error against the series and the maximum |second difference| of the profile. Run BTCS at the same step and report the same two numbers.
(c) Propose and test a remedy that keeps second-order accuracy for the rest of the run (hint: start with two BTCS steps of Δt/2, then switch).
(d) Independently of the numerics: Ti-6Al-4V's conductivity is about 7 W/(m·K) at 300 K and about 10.5 W/(m·K) at 700 K. Using the mean value over the range as a crude correction, estimate how much earlier the mid-plane reaches 450 K, and compare that modelling error with every numerical error in this problem set. Which should the colleague work on?

Answer form: three numbers for (a), a 2-row table for (b), a 1-row table for (c), one paragraph with a number for (d).
