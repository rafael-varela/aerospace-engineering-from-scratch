# P08.19  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. V = πD²/4 × (0.18 + 0.03). 2. `predict_drain_time(0.18, 0.03, 0.61, 0.03, d, 1e-4, D, 5e-4, u_h=5e-4)`; read the budget. 3. `jet_reynolds([0.18, 0.03], d)`.

## 2.
(a) Bernoulli surface → vena contracta; continuity A V_s = C_d a V_j; V_j² (1 − r²) = 2gh with r = C_d a/A; then A dh/dt = −C_d a V_j; integrate. Evaluate r²/2. (b) Contraction (area) times velocity loss (speed); a rounded entry removes the contraction. (c) Geometry of a tube along a tilted board; V ∝ (L sin θ)^½; partial derivatives. (d) Laplace pressure σ(1/R₁ + 1/R₂) for a cylinder and a sphere of radius d/2; divide by ρg.

## 3.
(a) For each run: `fit_discharge_coefficient(t, h, 0.5e-3, 0.08, 4e-3, 1e-4, 0.09, 5e-4)`. (b) `en_number` with the mean and one run's u (the geometry is shared); time difference of marks 200 and 30 against `predict_drain_time`. (c) 0.62 (4.05/4.0)² (90.0/90.4)². (d) L = (right − left) − (zero right − zero left); `airspeed_with_uncertainty(L, 0.6e-3, angle, 0.5, rho, 0.01 rho)`; compare dryer high upright and tilted heads with `en_number`; flag rows with Δh < u_Δh.

## 4.
(a) p_sc + ½ρV_sc² = p_ts + ½ρV_ts² with V_sc = V_ts/CR; Δh = L sin 15°; Δp; V_ts; uncertainty as in Problem 2(c); compare with √(2Δp/ρ). (b) ISA ρ at 3000 m; Δp = ½ρV²; Δh = Δp/((ρ_w − ρ)g).

## 5.
(a) Loop over angles with `airspeed_with_uncertainty`. (b) n = (single-reading u_V/V ÷ 0.1)²; think about quantisation and systematic parallax. (c) Find the angle and the angle uncertainty that give 5 %; how to measure an angle to 0.1° with a ruler.

## 6.
(a) Re range of each hole; where sharp-orifice C_d is known to be nearly constant; a second bottle with a 1 mm (or 2 mm) hole. (b) Probe Re against the threshold of about 100. (c) The fan's Re says its flow is turbulent; what matters is the jet's unsteadiness and its decay with distance; the procedure fixed the distance and averaged.
