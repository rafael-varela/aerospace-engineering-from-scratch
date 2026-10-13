# P00.13  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. Write E = m g h with g = 9.81 m/s² and h = 1.5 m; evaluate for each of the four masses, keeping units.
2. Write v = √(2 g h); evaluate once (it is the same for all four).
3. For h = 2.0 m, check each object against three rows of the envelope: mass ≤ 2 kg, height ≤ 2 m, m g h ≤ 40 J. An object is inside only if all three hold.

## 2.
1. (a) Substitute the mass limit and the height limit into m g h together, with g = 9.80665 m/s². Compare with 40 J. Argue that any smaller mass or height gives a smaller energy, so the energy row never binds for a dropped object.
2. (b) Set ½ m v² = 40 J and solve for v = √(2 E / m). Evaluate for m = 0.5 kg and m = 0.05 kg.
3. Compare each v with 5 m/s. Argue which limit binds for light fast objects and which for heavy slow ones; conclude why both are needed.

## 3.
1. (a) E = ½ k x²; then ½ m v² = E gives v = √(2 E / m) = x √(k / m).
2. Compare v with the 5 m/s launched-object limit. If it exceeds, compute the stretch x that gives exactly 5 m/s: x = v √(m / k). That is your "substitute" control.
3. (b) Hazards: go down the energy table (elastic, kinetic / projectile) and add eye and pinch. Controls: eliminate, substitute (the reduced x), engineer (rail, aim, tether), administer (downrange exclusion zone, "clear" call), protect (glasses). Name the hazard beside each control.
4. (c) Each stop criterion must be a sentence with a visible or audible test. Check each by asking whether two bystanders could disagree.

## 4.
1. (a) Run the checker and copy its lines: which rows are MISSING (absent or placeholder), which limits were parsed, the RESULT.
2. (b) Review each cell for validity: numbers with units in Energy limits, controls that name a hazard, binary stop criteria, a complete hazard list, a location that says what matters.
3. (c) In `check_safety_box`, after parsing limits: if the Energy limits row exists and is not a placeholder but `parse_energy_limits` returns an empty dictionary, append the agreed string to `missing`. Write a fixture whose Energy limits cell has words only and assert the string is in `report.missing`. Re-run the whole file.
4. (d) Rewrite every cell; re-run; confirm PASS.

## 5.
1. (a) E = Q · V with Q in mA·h, V in volts, times 3.6 J per mA·h·V. Divide by 40 J.
2. (b) Separate stored energy from the rate and path of release. Identify the engineered control inside the product. Contrast with gravity, which has an open path and an immediate release.
3. (c) Energy limits must carry numbers with units (voltage, current, mass of anything that moves, tip speed if a fan). Stop criteria: at least three observable conditions including battery heat or swelling. Disposal: where the power bank goes when it fails.

## 6.
1. (a) r = diameter / 2; ω = 2π n / 60; v = ω r. Compare with 5 m/s. Check the voltage and current against the envelope as well.
2. State which single quantity would need a higher limit and what controls justify it (guard both ends, disconnect before reaching in, fixed mounting).
3. (b) Write the nine rows; run `check_safety_box(text)` and record the flag; run again with `dict(HOME_LAB_ENVELOPE, speed_m_s=...)` and record PASS; write the justification.

## 7.
1. Read rules 1 to 11 of the NAR code; for each, write the row(s) of the safety box it answers.
2. Pick one rule for its *form* (a numeric table, a waiting time, an interlock) and say what the equivalent would be for a home experiment.
3. Two reasons it is not this course's code: scope (what it governs) and the activities this course never requires.
