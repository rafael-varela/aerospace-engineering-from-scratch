# P00.13  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every number.

## 1. (Warm-up: the energy inventory)

A tennis ball (m = 0.058 kg), a baseball (0.145 kg), a 0.5 kg bag of rice, and a 2.0 kg bag of flour are each held at h = 1.5 m above the floor.

Find, for each: the stored energy m g h in joules and the impact speed √(2 g h) in m/s (ignore drag). Then state which of the four would still be inside the envelope at h = 2.0 m. Answer form: a four-row table with units, then one sentence.

## 2. (Derivation: the shape of the envelope)

The envelope caps the falling mass at 2 kg, the height at 2 m, and the stored energy at 40 J.

(a) Show that the mass and height limits alone already keep m g h below 40 J, so the energy row is redundant for a *dropped* object.
(b) Now consider a thrown or launched object, which has kinetic energy that does not come from a drop. Derive the speed v at which a 0.5 kg object reaches 40 J, and the speed at which a 0.05 kg object does. Compare both with the 5 m/s speed limit and explain in two sentences why the course needs a speed limit *and* an energy limit rather than either alone.

## 3. (Unfamiliar scenario: a rubber-band catapult)

A later lesson launches a foam glider with a rubber band that behaves like a spring with k = 60 N/m, stretched by x = 0.25 m.

(a) Compute the elastic energy ½ k x² and the launch speed of a 0.020 kg glider if all of it becomes kinetic energy.
(b) Write the Hazards row and the Controls row of a safety box for this launch, with the controls in hierarchy order and each control naming the hazard it addresses.
(c) Write three stop criteria. For each, say how a bystander would know it was met.

## 4. (Numerical: uses your code)

Save the following text as `bad-box.md` and run `python starter/main.py bad-box.md`.

```
| | |
|---|---|
| **Hazards** | falling mass. |
| **Energy limits** | a small ball from head height, 9 V battery. |
| **Controls** | be careful. |
| **Supervision** | … |
| **Location** | kitchen. |
| **Stop criteria** | if it looks dangerous. |
| **Disposal** | TBD |
| **Not permitted in this project** | ladders. |
```

(a) List every finding the checker reports.
(b) List every problem with this box that the checker does *not* report, and say for each why a text checker cannot see it.
(c) Add one honest improvement to the checker: when the Energy limits cell contains no parsable number at all, report `"Energy limits (no numeric limit)"` in `missing`. Write a test for it in `tests/test_main.py`, make it pass, and confirm all the existing tests still pass.
(d) Rewrite the box so that it passes *and* fixes every item from (b).

## 5. (Judgment: a battery that stores a thousand times the envelope)

Look up (or read off the label of) your own phone's battery: a capacity Q in mA·h and a nominal voltage V (typically 3.8 V). Its stored energy is E ≈ Q · V, converted to joules (1 mA·h × 1 V = 3.6 J).

(a) Compute E for your phone. Compare it with the 40 J stored-energy limit.
(b) The course permits a phone on a tripod in the drop experiment. Explain, in terms of *energy* and *path*, why an object storing roughly a thousand times the envelope is permitted while a 2.5 kg mass at 2 m is not.
(c) Write the Energy limits row, the Stop criteria row, and the Disposal row for a project that powers a small fan from a USB power bank (5 V, up to 2 A).

## 6. (Design and defence: a guarded fan)

A Phase 10 project uses a 120 mm computer fan (12 V, 0.3 A, up to 2000 rev/min) inside a cardboard duct to measure the drag on a model.

(a) Compute the blade tip speed in m/s and compare it with the 5 m/s speed limit. Is the fan permitted under the default envelope? Under what single change to the envelope *and* which added controls could a project justify it?
(b) Write the complete nine-row safety box for this project. Run the checker on it with the default envelope, note which quantity it flags, then run it with an `envelope` dictionary that raises only that one quantity and explain in writing why the raise is justified by the controls you wrote.

## 7. (Reading a code: the NAR safety code as a model)

Read `@nar-safety-code` rules 1 to 11. Produce a two-column table mapping each rule to the row(s) of our safety box it most resembles, and name one rule whose *form* (not its content) you would want to copy into a home-experiment code and why. Then state, in one sentence each, two reasons why the NAR code itself is *not* the safety code for this course.
