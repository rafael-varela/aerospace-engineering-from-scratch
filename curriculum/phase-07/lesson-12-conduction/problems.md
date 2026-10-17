# P07.12  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer, say which inputs are illustrative, and give the observed order for every convergence claim.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

Data: `MATERIALS` (Lienhard Tables A.1 and A.2, 20 °C). Heat-transfer coefficients, loads and temperature limits are typical, illustrative values.

## 1. (Warm-up: a regeneratively cooled wall)

A chamber wall has 1 mm of Inconel X-750 on a copper liner of (2 + s/10) mm; gas side h = 3000 W/(m²·K) at 3000 K, coolant side h = (10 000 + 1000 s) W/(m²·K) at 300 K.

(a) Compute the four resistances per unit area, their shares, the heat flux and the three surface/interface temperatures.
(b) Repeat without the Inconel layer. Which layer would you remove first to cool the hot face, and why does the copper thickness hardly matter?

Answer form: a table and four numbers; two numbers and a sentence.

## 2. (Derivation)

(a) Write the energy balance of an interior control volume of the finite-volume model and show it equals k t V times the five-point Laplacian.
(b) Prove that, for any grid, the heat out through the fixed-temperature nodes plus the face losses equals the generation plus the edge inflows, using only the symmetry and zero row sums of the conduction matrix.
(c) Derive the fin equation θ″ = m²θ for a sheet with face losses, m² = h_f/(k t), and its solution with an adiabatic tip.
(d) By a von Neumann (Fourier-mode) argument, show the explicit 1D scheme is stable only for Fo ≤ ½.

Answer form: four derivations.

## 3. (Unfamiliar scenario: an insulating layer over aluminium)

A 10 mm layer of fused silica (Lienhard Table A.2) is bonded to 2 mm of 6061 aluminium. The silica face receives 100 kW/m² (illustrative, a short heating pulse) and the aluminium's back is adiabatic; everything starts at 300 K. (Fused silica stands in for a heat-shield material here; real tiles are porous and far more insulating.)

(a) Build the layered 1D system (k per interval, ρc per node; put the interface on a node) and march it with Crank–Nicolson for 30, 60 and 120 s. Report the hot-face and aluminium temperatures.
(b) Compare the hot-face temperature with the semi-infinite constant-flux result T_s = T_i + 2q″√(t/π)/√(kρc). When does it stop being accurate, and why?
(c) Check the energy budget against q″ t.

Answer form: a table; three comparisons and a sentence; a check.

## 4. (Numerical, uses your code: your own avionics card)

The lesson's card with a processor of (8 + s/5) W and a chassis wall at (60 + s/2) °C.

(a) Find the hot spot with copper planes only on n = 41, 81, 161 grids.
(b) Find the 6061 frame thickness that brings it to 100 °C.
(c) Show the energy balance closes and explain why the hot spot hardly moves with the grid even though the chip footprint has sharp edges.

Answer form: three temperatures; a thickness; a residual and a sentence.

## 5. (Numerical, uses your code: time accuracy)

A 50 mm copper slab (401 nodes) with convection h = 3000 W/(m²·K) at 3000 K from 300 K, to t = 1 s.

(a) With a reference solution from Crank–Nicolson at Δt = 10⁻⁴ s, find the maximum error for Δt = 0.04, 0.02, 0.01, 0.005 s with implicit Euler and with Crank–Nicolson, and the observed orders.
(b) Explain why implicit Euler shows a clean order 1 and Crank–Nicolson does not show a clean order 2 at these steps. What is the Fourier number at Δt = 0.04 s?
(c) Which scheme would you use to find the heat-sink firing time, and with what Δt?

Answer form: a table with orders; an explanation and a number; a choice with a reason.

## 6. (Judgment: how much should you trust the heat-sink time?)

(a) Recompute the 25 mm copper firing time for h = 1500 and 4500 W/(m²·K). How does the time scale with h?
(b) Copper's k falls from 401 W/(m·K) at 0 °C to 366 at 600 °C (Lienhard Table A.1). Does using 398 W/(m·K) make the predicted time long or short? Estimate the size of the effect with the effusivity argument.
(c) The test engineer proposes a 60 mm wall "to double the run time". What do you tell them?
(d) List the three inputs you would most want measured before a test, in order.

Answer form: two times and a scaling; a direction and an estimate; a sentence with a number; a list.
