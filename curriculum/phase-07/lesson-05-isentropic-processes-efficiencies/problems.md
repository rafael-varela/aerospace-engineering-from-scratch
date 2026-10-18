# P07.05  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact, which are closed forms and which are numerical. Write the steady-flow energy equation and the entropy balance for the device before any other equation.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

Gas data: the NASA dry-air model of the lesson (`AIR` in your module), or constant cp with γ = 1.4 and R = 287.05 J/(kg·K) where a problem says so. Machine efficiencies are typical, illustrative values, not data for a particular engine.

## 1. (Warm-up: one axial stage)

A single axial compressor stage takes air at 288.15 K and raises its pressure by π = 1.2 + s/100 with an isentropic efficiency of 0.90.

(a) Find the isentropic and actual exit temperatures with constant cp and with NASA air.
(b) Find the work per kilogram, the entropy generated per kilogram and T₀ s_gen with T₀ = 288.15 K.
(c) Find the stage's polytropic efficiency. Why is it so close to 0.90?

Answer form: four temperatures; three numbers per gas; a number and a sentence.

## 2. (Derivation: isentropes and polytropic efficiency)

(a) Starting from T ds = dh − v dp, derive s°(T₂ₛ) − s°(T₁) = R ln(p₂/p₁) for an ideal gas with cp(T), and reduce it to T₂ₛ/T₁ = π^((γ−1)/γ) for constant cp.
(b) From the definition dh = v dp/η_p, derive s°(T₂) − s°(T₁) = (R/η_p) ln π for a compressor, and the constant-cp closed form η_c = (π^m − 1)/(π^(m/η_p) − 1). Derive the turbine's version from dh = η_p v dp.
(c) Show that along a polytropic compression ∫T ds = (1 − η_p) w and s_gen = R ln π (1/η_p − 1), for any cp(T).
(d) Prove that η_c < η_p for every π > 1 and η_p < 1 (constant cp). Hint: with x = π^m > 1 and k = 1/η_p > 1, compare x^k − 1 with k(x − 1).

Answer form: four derivations.

## 3. (Unfamiliar scenario: an electric cabin-air compressor)

Some airliners pressurise the cabin with electrically driven compressors instead of air bled from the engines. At cruise, outside air at 216.65 K and 22 632 Pa (ISA, 11 km) is compressed to 75 kPa with an isentropic efficiency of 0.80 (typical, illustrative), at 1.2 kg/s. Ram compression in the inlet is ignored.

(a) Find the exit temperature, the specific work and the electrical power (motor and drive losses ignored), with NASA air.
(b) Find the entropy generated and the exergy destroyed per kilogram, using the outside air as the dead state (T₀ = 216.65 K).
(c) What polytropic efficiency does the compressor have? The cabin needs air near 290 K: what must happen to the compressed air before it enters the cabin, and what does that say about the work spent heating it?

Answer form: a temperature, a work and a power; two numbers; a number and two sentences.

## 4. (Numerical, uses your code: a multi-stage compressor)

A compressor has n = 8 + (s mod 7) stages of equal pressure ratio and the same isentropic efficiency, 0.88 each, and an overall pressure ratio of 20 + s, from 288.15 K (NASA air).

(a) Compute the exit temperature, the overall isentropic efficiency and the overall polytropic efficiency with `multistage_compressor`.
(b) Compute the exit temperature of a compressor with polytropic efficiency 0.88 and the same overall π. How far is your machine from it?
(c) How many stages of 0.88 are needed to come within 1 K of that limit?
(d) Compare the work of the first and the last stage. Explain the difference in one sentence.

Answer form: three numbers; a temperature and a difference; an integer; two works and a sentence.

## 5. (Numerical: a high-pressure turbine)

A high-pressure turbine expands gas from 1600 K through an expansion ratio of 4 with a polytropic efficiency of 0.90.

(a) Treating the gas as NASA dry air, find the exit temperature, the isentropic efficiency and the work per kilogram.
(b) Repeat with constant cp at γ = 1.4 and at γ = 1.33 (same R). Which is closer to the variable-cp answer, and why? What is γ of NASA air at 1600 K?
(c) Why is the turbine's isentropic efficiency higher than 0.90 while a compressor's is lower?

Answer form: three numbers; two pairs of numbers, a choice and a value; a sentence.

## 6. (Judgment: what limits the answer?)

The lesson's core compressor (π = 40, η_p = 0.90, 288.15 K) is being used to size a combustor.

(a) Compute the delivery temperature, work and isentropic efficiency for η_p = 0.88, 0.89, 0.91 and 0.92 (NASA air). Compare the spread with the difference between constant and variable cp.
(b) For which purposes is the constant-cp model adequate here, and for which is it not? Use numbers from the lesson's sweep.
(c) A test-cell engineer derives η_s from measured inlet and delivery temperatures. Name two physical effects that would make that "temperature-rise" efficiency differ from the efficiency defined by the shaft work, and say in which direction each moves it.

Answer form: a table of four rows and a comparison; two lists; two effects with directions.
