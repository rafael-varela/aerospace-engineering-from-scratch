# P04.15  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
(a) r_total = 10(1 + s/10) kΩ; `divider_design(9, 3, r_total)`; round top up, bottom down; `divider_output(9, R_t, R_b)`; I = 9/(R_t + R_b); P = 9²/(R_t + R_b). (b) `lsb(3.3, 12)` × (R_t + R_b)/R_b. (c) `led_resistor(3.3, 2.2, i)`; round up; `led_current`.

## 2.
(a) Error e uniform on [−LSB/2, LSB/2]; mean 0; E[e²] = (1/LSB)∫e² de = LSB²/12. (b) V_th = V_in R_b/(R_t + R_b), R_th = R_tR_b/(R_t + R_b); reading = V_th R_m/(R_th + R_m); error = −R_th/(R_th + R_m); R_th = R/2 for an equal divider. (c) V_out = V_ref R_f/(R_th + R_f), raw/FS = V_out/V_ref; solve for R_th. (d) Set the derivative R_lo/(R_lo + R)² − R_hi/(R_hi + R)² to zero; simplify to R² = R_loR_hi.

## 3.
(a) k = V_max/3.0; R_b = 100 kΩ/k, R_t = 100 kΩ − R_b; nearest E12 each; full scale = 3.3k; LSB = k × 3.3/4096. (b) Evaluate k with (1.01R_t, 0.99R_b) and (0.99R_t, 1.01R_b); RSS via `propagate`. (c) `thevenin_resistance(R_t, R_b)`. (d) The bus is not the ADC reference: V_ref errors enter directly.

## 4.
(a) `analyse_thermistor(make_synthetic_bench_data(s))`. (b) u_raw = √((std/√16)² + (16/√12)²); f(raw, R_f) = `temperature_from_u16(raw, R_f, sh)`; `propagate` with (u_raw, 0), (0, u_R), (u_raw, u_R). (c) `adc_calibration_fit(pot_u16, pot_dmm_V)`; compare residual with 3.3/4096 V. (d) Label.

## 5.
(a) R_b′ = V_rR_t/(V_in − V_r); R_m = R_b′R_b/(R_b − R_b′). (b) `dmm_standard_uncertainty` for V_r and V_in; `propagate` R_m(V_r, V_in); then with R_t, R_b as extra inputs. (c) Make R_th comparable to R_m (larger resistors make the reading move more), and measure R_t, R_b with the meter's ohms range first, so their tolerance is replaced by the meter's resistance spec.

## 6.
Map each symptom to a component: resistor value (unpowered ohms reading), LED V_f or shunt (voltage across the resistor without the ammeter), ADC nonlinearity (calibrate the middle only), self-heating (vary the divider current or the stirring), scaling vs resolution (look at the low bits), ratiometric vs absolute (correlate drift with the reference).
