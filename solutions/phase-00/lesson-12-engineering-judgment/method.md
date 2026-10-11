# P00.12  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. (a) Multiply the three dimensions, multiply by 1.2 kg/m³. 2. (b) Multiply seats by about 100 kg. 3. (c) $\tfrac12 m V^2$ in joules; divide by $43 \times 10^6$ J/kg to get kilograms of fuel equivalent. 4. For each, round to the nearest power of ten (use `order_of_magnitude`) and separately to one significant figure. 5. State the uncertainty factor of each input and combine them.

## 2.
1. $d(T/W)/dC_L = -C_{D0}/C_L^2 + k = 0$. 2. Solve for $C_L^*$; substitute back for $(T/W)_{\min}$; invert for $(L/D)_{\max}$. 3. Write $C_L^* = 2W/(\rho V^{*2} S)$ and solve for $V^*$. 4. Read off: $V^* \propto \sqrt{W/S}$, $\propto 1/\sqrt{\rho}$, independent of $C_L$'s own value except through $C_{D0}, e, AR$. 5. Evaluate with ISA density at 11 km.

## 3.
1. Units: $V/c$ → (m/s)/(kg/(N·s)) = m·N/kg = m·(kg·m/s²)/kg = m²/s². Not a length. 2. $V/(cg)$ → divide by m/s²: m. A length. 3. Formula (c) has $W_i/W_f$ instead of its logarithm: dimensionally fine but gives nonzero range with no fuel burned. 4. Evaluate (b). 5. Compare with 3000–15 000 km for airliners.

## 4.
1. Compute mass and area from the seed. 2. Loop over speeds, call `check_airliner`, store `report.ok`, $T/W$, $C_L$, Mach. 3. Passing interval = contiguous speeds with `ok`. 4. Minimum of the $T/W$ array; analytic $V^*$ from problem 2 with the seed's wing loading. 5. Plot with `matplotlib` (`Agg`), mark the analytic minimum, save.

## 5.
1. Run the checker; it will flag mass, altitude, Mach, wing loading, probably $T/W$. 2. For each, compute the value by hand to confirm it is a true statement. 3. Build the UAV table: mass 1–50 kg, altitude 0–3 km, Mach 0.03–0.15, wing loading 50–400 N/m², $C_L$ 0.4–1.0, $T/W$ 0.08–0.3 (justify each from drone practice: low speed means high $C_L$ and low $L/D$). 4. Apply it to the drone by hand or by temporarily substituting the table.

## 6.
1. Run both datasets, list the flags. 2. Confirm each flag describes the aircraft correctly. 3. Choose an option. 4. For the cost, construct one error that your choice lets through (for (i)) or one burden it adds (for (ii), (iii)), and one error it still catches.

## 7.
1. Identify the quantity: the velocity change per firing, $\Delta v = I/m$ where $I$ is impulse. 2. Reference value: rated thrust × commanded duration / spacecraft mass, from the propulsion specification. 3. Comparison: the $\Delta v$ the navigation filter needed to fit the tracking data. 4. Note the factor between them and that it is constant, which fingerprints a unit conversion. 5. Explain why dimensions cannot distinguish lbf·s from N·s.
