# P00.13  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. g = 9.81 m/s² unless the problem says otherwise.

## 1. The energy inventory

**Assumptions.** Point masses released from rest; air drag ignored (a limiting-case idealisation; `P00.15` shows it changes the fall of a tennis ball by a few percent, far less than the margins here).

**Derivation.** E = m g h with h = 1.5 m, so E = m × 14.7 J/kg. The impact speed is v = √(2 g h) = √(2 × 9.81 × 1.5) m/s = 5.42 m/s for every object, because the mass cancels.

**Result.**

| Object | m (kg) | E at 1.5 m (J) | v at impact (m/s) |
|---|---|---|---|
| tennis ball | 0.058 | 0.85 | 5.42 |
| baseball | 0.145 | 2.13 | 5.42 |
| rice | 0.5 | 7.4 | 5.42 |
| flour | 2.0 | 29.4 | 5.42 |

At h = 2.0 m: every mass is ≤ 2 kg, the height is exactly the 2 m limit (inclusive), and the largest energy is 2.0 × 9.81 × 2.0 = 39.2 J ≤ 40 J. All four are inside the envelope, the flour bag with 0.8 J to spare.

**Check.** Dimensions: kg·(m/s²)·m = J. Flour at 1.5 m is 2.0/0.5 = 4 times rice, and 29.4/7.4 = 4.0.

**Common mistakes.** Using grams without converting (gives millijoules by accident). Listing the 5.42 m/s impact speed against the 5 m/s envelope row: that row is for *powered or launched* objects; a dropped object is governed by the height and energy rows.

## 2. The shape of the envelope

**Assumptions.** g = 9.80665 m/s² (the value the checker uses).

**(a)** With m = 2 kg and h = 2 m, m g h = 2 × 9.80665 × 2 = 39.2 J < 40 J. Since m g h increases with m and with h, any permitted mass and height gives at most 39.2 J. For a dropped object the energy row can never bind before the mass or height row does. (It exists because other projects store energy in springs or batteries, where there is no mass-times-height to cap.)

**(b)** ½ m v² = 40 J gives v = √(2 × 40 / m):
- m = 0.5 kg: v = √160 = 12.6 m/s.
- m = 0.05 kg: v = √1600 = 40 m/s.

Both are far above 5 m/s. For a light object the energy limit alone would allow a 40 m/s pellet, which is a serious eye hazard; the speed limit is what stops it. For a heavy object the speed limit alone would allow 2 kg at 5 m/s, which is 25 J and fine, but a 4 kg mass at 5 m/s is 50 J and is stopped by the energy (and mass) rows. Energy measures how much damage is possible; speed measures how a light object can concentrate it. Each catches a case the other misses.

**Check.** Rearranging v² = 2E/m: halving m by ten should multiply v by √10 ≈ 3.16, and 40/12.6 = 3.17.

**Common mistakes.** Forgetting the ½. Concluding that the energy limit is "useless" because it never binds for drops.

## 3. The rubber-band catapult

**Assumptions.** The band is an ideal spring (k = 60 N/m) and all elastic energy becomes the glider's kinetic energy (no losses, no energy left in the band). Real bands lose a third or more, so this is the *worst case*, which is the right case for a safety calculation.

**(a)** E = ½ k x² = 0.5 × 60 × 0.25² = 1.88 J. v = √(2 E / m) = √(2 × 1.875 / 0.020) = √187.5 = 13.7 m/s. This exceeds the 5 m/s launched-object limit by nearly a factor of three, so the launch as described is not permitted at home. The substitute control: the stretch that gives exactly 5 m/s is x = v √(m / k) = 5 × √(0.020 / 60) = 0.091 m. A stop block that limits the pull to 90 mm makes the design compliant.

**(b)**
- **Hazards:** elastic energy in the stretched band (snap to hand or eye on release or breakage); projectile (glider at up to 5 m/s with a hard nose; bystander and eye); pinch (fingers in the release hook); trip (launch rail on the floor).
- **Controls, in hierarchy order:** substitute: stretch limited to 90 mm by a fixed stop (projectile speed); engineer: launch rail bolted to a board and aimed at a wall or a hung sheet, foam nose on the glider, band inspected and replaced at the first nick (band snap); administer: a 5 m downrange exclusion zone marked on the floor, one person launches and calls "clear" first (projectile to bystander); protect: safety glasses for everyone in the room (eye).

**(c)** Stop criteria, each with a bystander test:
1. "Anyone is downrange of the rail inside the marked zone." Test: look at the tape on the floor.
2. "The band shows a nick, a white stress line, or has been stretched past the stop." Test: look at the band before each shot; the stop is a physical block.
3. "The glider has struck anything other than the sheet or the wall." Test: where the glider is lying.
4. "The launcher is tired or has been asked to hurry." Test: the launcher says so; the rule is that saying so ends the session.

**Check.** Energy at the reduced stretch: ½ × 60 × 0.091² = 0.25 J = ½ × 0.020 × 5² ✓.

**Common mistakes.** Computing with x in centimetres. Writing "be careful with the band" as a control. Writing "if it looks unsafe" as a stop criterion.

## 4. Numerical, with the checker

**(a)** The checker reports:
```
MISSING   Supervision          (cell is the placeholder "…")
MISSING   Regulations          (row absent)
MISSING   Disposal             (cell is "TBD")
LIMIT     voltage_v   9  (max 12)  ok
RESULT    FAIL
```
Exactly three rows are missing. The one number present (9 V) is inside the envelope.

**(b)** Not reported, and why a text checker cannot see it:
- "a small ball from head height": no number, so no envelope check happened at all. The checker only compares numbers it can parse.
- "be careful": a control that controls nothing. Validity of a sentence is a judgment.
- "if it looks dangerous": not binary. Same reason.
- Hazards lists only "falling mass"; the bounce (eye, equipment) and the trip hazard of the measuring tape are missing. The checker cannot know what hazards the experiment has.
- "kitchen": says nothing about the floor, the clearance, or who else uses the room.

**(c)** In `check_safety_box`, after `limits = parse_energy_limits(limits_cell)`:
```python
if "Energy limits" not in missing and not limits:
    missing.append("Energy limits (no numeric limit)")
```
Test:
```python
def test_energy_limits_without_numbers_is_a_finding():
    text = COMPLETE.replace("mass 0.5 kg, drop height 1.2 m, 9 V battery, 0.3 A, 3 m/s.", "a small ball from head height.")
    assert "Energy limits (no numeric limit)" in impl.check_safety_box(text).missing
```
All existing tests still pass because every complete fixture carries numbers.

**(d)** A passing rewrite (one of many):
```
| **Hazards** | falling mass (foot); bounce (eye, phone on tripod); trip (tape measure on floor). |
| **Energy limits** | mass 0.20 kg, drop height 1.8 m (3.5 J), phone battery 3.8 V, no other electricity. |
| **Controls** | mat under the drop point (falling mass); glasses (bounce); tape measure taped flat and removed after marking (trip); phone 1.5 m from impact (bounce). |
| **Supervision** | none needed for an adult; children only with an adult present. |
| **Location** | indoors on a hard floor with a mat, 1 m clear of furniture, door closed. |
| **Regulations** | none apply. |
| **Stop criteria** | anyone inside the 1 m circle; ball bounces toward a person; any damage to the floor; experimenter tired or rushed. |
| **Disposal** | none (no consumables). |
| **Not permitted in this project** | ladders, chairs, windows; metal or glass objects; drops onto people or pets. |
```

**Common mistakes.** Fixing only what the checker reported and declaring the box done.

## 5. The phone battery

**Assumptions.** E ≈ Q V with nominal voltage; 1 mA·h × 1 V = 3.6 J. Example phone: Q = 4000 mA·h, V = 3.8 V.

**(a)** E = 4000 × 3.8 × 3.6 J = 54.7 kJ, about 1400 times the 40 J limit. (Your own phone will differ; the ratio is of order 10³ for any modern phone.)

**(b)** The envelope caps the energy that has an *open path*: a raised mass has gravity pulling on it continuously and the only control between it and a foot is a hand. The phone's chemical energy has no path to a person unless the cell is punctured, shorted, crushed or over-charged, and the phone's protection circuit (an engineered control) limits the current that can flow out to a few amperes at 5 V. The energy is large; the rate at which it can be released through any permitted use is tiny, and the stored-energy row of the envelope is about energy that *can move*. That is also why the course adds handling rules (no soldering to cells, no unattended charging) and a stop criterion (heat, swelling): those are the abuses that open a path.

**(c)**
- **Energy limits:** 5 V, 2 A maximum from a protected USB power bank; fan 80 mm diameter, 1500 rev/min, tip speed 6.3 m/s inside a guard; no object heavier than 0.1 kg moves.
- **Stop criteria:** power bank warm to the touch, swollen or smelling; any rattle or rubbing from the fan; the guard displaced or a finger-sized gap appears; anything drawn into the inlet.
- **Disposal:** a failed power bank goes to a battery collection point in a non-flammable container; it is never put in household waste.

**Check.** 15.2 W·h × 3600 s/h = 54.7 kJ ✓.

**Common mistakes.** Concluding that the phone is "more dangerous" than the 2.5 kg mass and must be banned. Confusing mA·h (charge) with energy.

## 6. The guarded fan

**(a)** r = 0.060 m, ω = 2π × 2000 / 60 = 209 rad/s, v = ω r = 12.6 m/s. This is 2.5 times the 5 m/s limit, so under the default envelope the fan is not permitted. Voltage (12 V) and current (0.3 A) are inside. The single quantity that must be raised is the speed, to about 15 m/s, and only for a *guarded* rotor. Controls that justify it: finger-proof mesh (openings under 8 mm) on both ends of the duct, fixed with tape or screws (engineer); fan screwed to the duct so it cannot walk (engineer); power disconnected at the connector, not just switched off, before any hand enters the duct (administer); long hair and sleeves tied back (administer); glasses (protect).

**(b)** A complete box, then the two checker runs:

```
| **Hazards** | rotating blades (cut, entanglement); objects drawn into the inlet and ejected (eye); 12 V supply (wire heating in a short); duct tipping (pinch). |
| **Energy limits** | 12 V, 0.3 A; fan 120 mm at 2000 rpm, tip speed 12.6 m/s behind guards; nothing heavier than 0.05 kg moves. |
| **Controls** | mesh guards both ends (blades, ejection); fan screwed to duct, duct clamped to bench (tipping); fused 12 V supply, disconnected before reaching in (short, blades); hair and sleeves tied back (entanglement); glasses (eye). |
| **Supervision** | none needed for an adult; children only with an adult present. |
| **Location** | indoors on a stable bench, inlet 0.5 m from any wall, nothing loose within 1 m of the inlet. |
| **Regulations** | none apply. |
| **Stop criteria** | any rattle or rub; guard displaced or any gap over 8 mm; anything enters the duct; supply warm; experimenter tired. |
| **Disposal** | none routine; a failed fan goes to electronic-waste collection. |
| **Not permitted in this project** | running without both guards; hand-held operation; mains supply; speeds above 2000 rpm. |
```
With the default envelope: `OVER ENVELOPE speed_m_s 12.6 (max 5)`, `FAIL`. With `dict(HOME_LAB_ENVELOPE, speed_m_s=15.0)`: `PASS`. The written justification: the 5 m/s limit protects against contact with a moving part; with finger-proof guards on both ends and a disconnect-before-reach rule, contact is prevented by an engineered control rather than by speed, and the raised limit is still well below the speed of a desk fan sold without a safety mark.

**Check.** 2000 rev/min × 2π × 0.06 m / 60 = 12.57 m/s ✓.

**Common mistakes.** Using the diameter for r. Raising the whole envelope rather than one key. Treating the raised limit as the new default.

## 7. Reading the NAR code

| NAR rule | Safety-box row(s) it answers |
|---|---|
| 1 Materials | Not permitted (what may be built from) |
| 2 Motors | Not permitted; Energy limits (certified, bounded energy) |
| 3 Ignition System | Controls (electrical launch with interlock) |
| 4 Misfires | Stop criteria (a wait before approaching) |
| 5 Launch Safety | Controls; Supervision (countdown, distance, clearance) |
| 6 Launcher | Controls (rod, deflector, tip protection) |
| 7 Size | Energy limits (mass and impulse cap) |
| 8 Flight Safety | Location; Not permitted (no targets, no cloud, clear sky) |
| 9 Launch Site | Location (a distance table scaled to motor power) |
| 10 Recovery System | Controls (how the energy comes back down) |
| 11 Recovery Safety | Not permitted; Stop criteria (no recovery from power lines or dangerous places) |

Form worth copying: the launch-site table in rule 9, which scales a required distance to a quantitative measure of energy. The home-experiment equivalent is our envelope table, and a per-project version could scale an exclusion radius to m g h or to launch speed. Another is rule 4's fixed waiting time after a misfire: a stop criterion with a built-in, non-negotiable recovery step.

Why it is not this course's code: (1) its scope is model rocketry with certified motors, outdoor launches and recovery, none of which this course requires, so most rules have no object here; (2) it says nothing about the hazards the course does have (batteries, rotors, springs, drops, hot tools), so it would leave the real home-lab risks unaddressed.

## Check-yourself answers

1. A hazard is stored or moving energy with a path to a person (or anything else that can cause harm); a control is a measure that removes the energy, blocks the path, or reduces the harm.
2. E = 0.5 × 9.81 × 1.5 = 7.4 J; inside the envelope on all three rows (0.5 kg ≤ 2 kg, 1.5 m ≤ 2 m, 7.4 J ≤ 40 J).
3. Voltage sets the current through a body (shock), current sets the heating of a wire in a short circuit (fire). A low-voltage source that can deliver high current is a fire hazard with no shock hazard, so one limit does not cover the other.
4. v = 2π n r / 60 with n in rev/min and r in m, giving m/s.
5. (b) and (c): each has a visible test and is met or not. (a) needs a judgment and two people could disagree.
6. Lighter ball (substitute) → landing mat (engineer) → exclusion circle (administer) → glasses (protect). There is no "eliminate" option that still answers the question.
7. A false pass: the checker reports PASS for a box that a reviewer would reject. It is the checker's main limitation.
8. Any two of: warm or hot to the touch; swelling or bulging; a dent or puncture; a sweet or chemical smell; visible leakage.
