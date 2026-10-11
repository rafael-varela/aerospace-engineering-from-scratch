# P00.13  Safety & Standards for Home Experiments

> Every hazard is stored energy looking for a path. Name the energy, cap it, block the path, and decide in advance when you stop.

| | |
|---|---|
| **Learning outcomes** | 1. Identify the hazards in a proposed home experiment by listing every energy source and the path it could take to a person. 2. Compute the stored or kinetic energy of a setup and decide whether it sits inside the course's home-lab envelope. 3. Write controls in the order eliminate → substitute → guard → procedure → protective equipment, and write stop criteria that are observable and binary. 4. Complete the nine-row safety box used by every 🔬 lesson in this course and check it with a script. 5. Decide whether a proposed experiment is permitted at home, permitted with changes, or not permitted. |
| **Prerequisites** | `P00.01`, `P00.07`, `P00.10`; skill check: can compute m·g·h and ½·m·v² with SI units; can write a Python function that returns a list. |
| **Engineering question** | For the falling-object timing experiment of `P00.15` (a 0.2 kg ball dropped from up to 1.8 m, filmed with a phone), what can hurt someone, how much energy is involved, which controls remove the risk, and when must the session stop? |
| **Assumptions and limits** | Hazard analysis here is qualitative plus simple energy bookkeeping. It assumes an adult learner, an ordinary room, and the course's envelope (below). It does not replace a regulatory risk assessment, and it cannot judge things a text checker cannot see: a wobbly table, a tired experimenter, a curious child. |
| **Analytical work** | Energy inventory of a setup: gravitational, kinetic, elastic, electrical, thermal, chemical. Conversion between forms. Tip speed of a rotor. Why the envelope numbers are what they are. |
| **Computational deliverable** | `starter/main.py`: a checker that parses a `safety.md` table, reports missing or placeholder rows, extracts the numbers in "Energy limits" with their units, and flags anything above the envelope. Tests prove it on fixture strings: complete boxes pass, each missing row is named, units convert correctly, the boundary is inclusive, and a derived stored energy is checked. |
| **Practical activity** | **design-build**: write the complete safety box for the `P00.15` drop experiment before that lesson, run the checker on it, and fix every finding. No physical apparatus in this lesson. |
| **Equipment** | None (a computer). The safety box you write will be used with Level 1 equipment in `P00.15`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Your completed `safety.md` for the drop experiment, the checker's PASS report, your energy inventory table, and the problem set. |
| **Assessment** | `assessments/phase-00/grading-rubric.md`; competencies 1 (reasoning about energy), 3 (the checker), 4 (engineering interpretation: the permitted / not-permitted decision). A safety box with any missing row is not graded. |
| **References** | `@nar-safety-code` rules 1 "Materials" through 11 "Recovery Safety" (required-free): read it as a model of a *good* safety code, not as the rules for this course, which uses no motors. |
| **Extensions** | ★ Extend the parser to battery energy in W·h and to rotor tip speed from "rpm" and a diameter. ★ Write a safety box for a tabletop fan-driven wind tunnel and defend it. |
| **Time** | unmeasured |

## 1. The engineering problem

In `P00.15` you will drop a ball from several heights and time its fall with a phone camera. It sounds harmless, and done properly it is. But "done properly" has content. A 0.2 kg ball dropped from 1.8 m arrives at 5.9 m/s. A steel ball of the same size would be heavier and would dent a floor or a toe. A ball that bounces into a phone on a tripod breaks the phone. A tape measure stretched across the floor is a trip hazard for whoever walks in. None of these is a disaster, which is exactly why they get skipped. The habit this course needs is the habit of writing them down *before* every session, including the harmless ones, so that it is already automatic when the energies are larger: a tensioned rubber-band catapult in Phase 04, a fan in a cardboard wind tunnel in Phase 10, a lithium battery driving a motor in Phase 14.

The course handles this with one standard document: the **safety box**, a nine-row table that sits at the top of every project's `safety.md` and is copied into every report. The rows are fixed: Hazards, Energy limits, Controls, Supervision, Location, Regulations, Stop criteria, Disposal, and Not permitted in this project. A report whose safety box has an empty or placeholder row is not graded. That is a hard rule, so this lesson teaches you how to fill each row and gives you a script that refuses to pass an incomplete box.

The success criteria are concrete. You can list every energy source in a setup and put a number on each. You can show each number is below the **home-lab envelope**: the limits below which the course allows an experiment at home with no supervision. You can write a stop criterion that a stranger could apply without judgment. And you can run `python main.py my-safety.md` and get `PASS`.

A good model to study first is the National Association of Rocketry's Model Rocket Safety Code (`@nar-safety-code`). It governs a hobby with real stored energy, it has been refined for decades, and it is short. Notice its structure as you read: it names *what may be used* (materials, certified motors), *how* (an electrical launch system, a launch rod, a minimum angle), *where* (a launch-site size table scaled to the energy of the motor), *what to do when something goes wrong* (a wait after a misfire), and *how the energy comes back down* (a recovery system, no recovery from power lines). Every one of those ideas has a row in our safety box. We borrow the structure, not the rules: this course never requires combustion, a motor, or an outdoor flight.

## 2. Understand and derive

### 2.1 Hazard = energy + path

A **hazard** is something that can cause harm. Almost every physical hazard is stored or moving energy plus a path by which it can reach a person. The list of energy forms is short, and you should run down it every time:

| Energy form | Where it hides at home | Formula (SI) |
|---|---|---|
| Gravitational | a raised mass | E = m g h |
| Kinetic | anything moving, including a rotor tip | E = ½ m v² |
| Elastic | stretched rubber, a compressed spring, a bent ruler | E = ½ k x² |
| Electrical | batteries, chargers, mains | battery E ≈ V · Q; capacitor E = ½ C V² |
| Thermal | soldering irons, hot glue, resistors | heat flow, surface temperature |
| Chemical | glues, solvents, battery chemistry | not computed; identified and avoided |
| Pressure | inflated or pumped containers | E ≈ p V for a gas at modest overpressure |

Two more hazards do not fit the energy picture but belong in the Hazards row every time: **eye** (anything that can fly, snap, splash, or spray) and **pinch / cut** (clamps, sharp edges, knives).

**Worked example (the drop test).** A ball of mass m = 0.20 kg is released from h = 1.8 m. Its gravitational potential energy relative to the floor is

E = m g h = 0.20 kg × 9.81 m/s² × 1.8 m = 3.5 J.

Check the units: kg · m/s² · m = kg·m²/s² = J. At the floor, all of that is kinetic (ignoring air drag, which `P00.15` shows is a small correction for this ball), so ½ m v² = 3.5 J gives v = √(2 g h) = 5.9 m/s. Is 3.5 J dangerous? A dropped 350 g smartphone from 1 m carries about the same. It hurts a bare foot, it does not injure a shod one, and it can break a screen. Hazard: falling mass, foot. Hazard: projectile (bounce), eye and equipment.

### 2.2 Why the envelope has the numbers it has

The course's **home-lab envelope** is the set of limits inside which a 🔬 lesson may be done alone at home:

| Quantity | Limit | Reason |
|---|---|---|
| Mass that can fall | 2 kg | a 2 kg mass is the heaviest thing one hand can place and release in a controlled way |
| Drop height | 2 m | reach of an adult standing on the floor; no ladders, no windows |
| Stored energy m g h | 40 J | 2 kg from 2 m is 39 J; the same as a 1 kg object at 9 m/s, which bruises but does not fracture |
| Voltage | 12 V DC | below the threshold at which dry skin conducts enough current to be felt; ordinary batteries and USB supplies qualify, mains never does |
| Current | 2 A | enough for small motors and sensors; above it, thin wires heat and short circuits burn |
| Speed of a powered or launched moving mass | 5 m/s | a brisk walk; a 50 g part at 5 m/s carries 0.6 J. Applies to rotors, springs, catapults and anything driven; a *dropped* object's fall speed is already bounded by the height and energy rows and is not listed in the Energy limits cell |
| Surface temperature | 60 °C | hot to touch, no burn with a brief contact; hot-glue guns and soldering irons are out of this envelope and get their own project rules |
| Gauge pressure | 0 kPa | no pressure vessels at home, ever |

These are deliberately conservative. They are not the point where injury begins; they are the point below which a careful adult alone at home has margin. Later projects that need more (a guarded fan, a soldering iron) carry their own safety box with a justified, larger limit for that single quantity and extra controls. The checker accepts a different `envelope` dictionary for that case, but the default is this table.

**Why voltage and current are separate limits.** Harm from electricity is set by current through the body, and current is voltage divided by body resistance. Dry skin has a resistance of tens of thousands of ohms, so 12 V drives well under a milliampere: not felt. Mains at 230 V (or 120 V) drives tens of milliamperes, which can stop breathing. The separate current limit protects the *equipment*: a 2 A source into a short circuit heats a thin wire to the point of melting insulation. A battery that can deliver 20 A into a short (a lithium pouch cell can) is a fire hazard even at 3.7 V. That is why batteries get their own paragraph below.

**Why speed is its own limit.** Kinetic energy alone misleads for rotors. A 100 mm propeller at 3000 rev/min has a tip speed of

v = ω r = (2π × 3000 / 60 s⁻¹) × 0.050 m = 15.7 m/s,

far above 5 m/s, yet its kinetic energy is small because the blade is light. What hurts is the edge speed and the cut, not the joules. So an unguarded rotor is never permitted at home, and any rotor at all appears only inside a mesh enclosure in a project with its own safety box.

### 2.3 Controls: the hierarchy

A **control** is a measure that reduces risk. List them in a fixed order, because the order is also a ranking of how much they can be trusted:

1. **Eliminate** the hazard. Do not use a steel ball when a tennis ball answers the same question.
2. **Substitute** a lower energy. Drop from 1.2 m instead of 1.8 m if the camera still resolves the fall.
3. **Engineer** a barrier. A landing mat; a clamp instead of a hand; a mesh guard; a fuse; a tether on anything that could fly.
4. **Administer** a procedure. An exclusion circle on the floor; a "clear!" call before release; one person releases and the same person says when it is safe.
5. **Protect** the person. Safety glasses, closed shoes. Last, because it only works if worn, and only for the person wearing it.

Safety glasses are required for any experiment with a projectile, a snapping or flying part, a spring under load, a drill, or a battery being connected. For the drop test that means the ball could bounce toward your face: glasses on.

### 2.4 Stop criteria

A **stop criterion** is a condition that, when it is met, ends the session immediately and without discussion. It must be *observable* (you can see, hear, smell, or read it) and *binary* (it is met or it is not). "Be careful if the ball bounces oddly" is not a stop criterion. "Anyone, including a pet, enters the 1 m circle" is. So is "any reading outside the envelope", "the battery is warm to the touch", "a crack appears in the fixture", and "I am tired or rushed". Write the last one. It is the most common real cause of an accident.

When a stop criterion is met: stop, make the setup safe (put the mass down, disconnect the battery, release the spring), and only then think about what happened.

### 2.5 Two topics that need their own words

**Batteries.** Alkaline AA cells store roughly 10 kJ each and release it slowly; they are safe to handle, must not be shorted (a shorted AA gets hot enough to burn skin), and go to a collection point when dead. Lithium-ion cells (phone, power bank, hobby packs) store similar energy per cell but can release it in seconds if shorted, punctured, over-charged, or crushed, and the result is fire. Rules in this course: use lithium cells only inside their original protected product (a phone, a power bank with a USB socket) or a pack with a protection circuit; never solder to a cell; never charge unattended; a swollen, hot, or dented cell is a stop criterion and goes to a battery collection point in a metal or ceramic container, not in household waste.

**Moving parts.** Anything that rotates or reciprocates gets a guard the fingers cannot pass through before it is powered. Long hair and loose sleeves are tied back. The power is disconnected, not just switched off, before a hand goes near a mechanism.

### 2.6 Reading the NAR code as a code

Go through `@nar-safety-code` rules 1 to 11 and, for each, write the row of our safety box it corresponds to. You will find: materials and motors → Not permitted (only certified, only lightweight), ignition system and misfires → Controls and Stop criteria (an electrical system with a safety interlock, a wait after a misfire), launch safety → Controls and Supervision (a countdown, a minimum distance), launcher → Controls (a blast deflector, a launch rod with its tip protected), size and flight safety → Energy limits (a mass and impulse cap, no flight at targets or into clouds), launch site → Location (a distance table scaled to motor power), recovery → Controls and Disposal. The point of the exercise is to see that the nine rows are not bureaucracy. They are the questions every mature safety code has ended up answering.

## 3. Build the computational model

The deliverable is a **completeness checker**: a script that turns the rule "no report is graded with an incomplete safety box" into something a computer applies the same way every time. It is a text tool, not a physics model, so its tests are fixture strings.

### Step 1: Parse the table

A safety box is a Markdown table with two columns. Each data row is `| **Label** | text |`. We want a dictionary from label to cell text.

```python
def parse_safety_box(text: str) -> dict[str, str]:
    rows: dict[str, str] = {}
    for line in text.splitlines():
        line = line.strip()
        if not (line.startswith("|") and line.endswith("|")):
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        if len(cells) < 2:
            continue
        label = cells[0].strip("*").strip()
        if not label or set(label) <= set("-: "):   # blank header or |---|---| separator
            continue
        rows.setdefault(label, cells[1].strip())
    return rows
```

`strip("*")` removes the bold markers; the `set(label) <= set("-: ")` test discards separator rows. `setdefault` keeps the first occurrence if a label is repeated.

### Step 2: Find the missing rows

A row is missing if its label is absent, or if its cell still holds a placeholder such as `…` or `TBD`. Compare case-insensitively so `energy limits` and `Energy limits` are the same row.

```python
PLACEHOLDERS = frozenset({"", "…", "...", "tbd", "todo", "n/a", "na", "-", "?"})

def missing_rows(box, required=REQUIRED_ROWS):
    lowered = {k.lower(): v for k, v in box.items()}
    out = []
    for label in required:
        value = lowered.get(label.lower())
        if value is None or value.strip().strip(".").lower() in PLACEHOLDERS:
            out.append(label)
    return out
```

### Step 3: Read the numbers in "Energy limits"

The cell is free text: "mass 0.5 kg, drop height 1.2 m, 9 V battery, 500 mA". We want `{"mass_kg": 0.5, "drop_height_m": 1.2, "voltage_v": 9.0, "current_a": 0.5}`. A regular expression finds `number + unit` pairs; a table maps each unit to a key and a conversion factor. Order matters: `m/s` must be tried before `m`, and `mA` before `m` and `A`, or "500 mA" is read as 500 m. Keep the maximum when a quantity appears twice, because the envelope is a ceiling.

```python
_UNIT_TABLE = (
    (r"m/s", "speed_m_s", 1.0), (r"km/h", "speed_m_s", 1/3.6),
    (r"°\s*C", "temperature_c", 1.0), (r"kPa", "pressure_kpa", 1.0),
    (r"mA", "current_a", 1e-3), (r"kg", "mass_kg", 1.0), (r"g", "mass_kg", 1e-3),
    (r"mm", "drop_height_m", 1e-3), (r"cm", "drop_height_m", 1e-2), (r"m", "drop_height_m", 1.0),
    (r"V", "voltage_v", 1.0), (r"A", "current_a", 1.0),
)
```

Build one regex with a named group per unit, then for each match find which group fired. Add a negative look-ahead `(?![A-Za-z/])` after the unit so `m` does not match the start of `mass`.

### Step 4: Compare with the envelope

```python
def check_envelope(limits, envelope=HOME_LAB_ENVELOPE):
    exceeded = [k for k, v in limits.items() if k in envelope and v > envelope[k]]
    if "mass_kg" in limits and "drop_height_m" in limits:
        if limits["mass_kg"] * G0 * limits["drop_height_m"] > envelope["energy_j"]:
            exceeded.append("energy_j")
    return exceeded
```

The comparison is strict (`>`), so a value equal to the limit passes: the envelope is inclusive. The stored-energy line is the one piece of physics in the tool: a mass and a height together imply an energy, and the energy has its own limit.

### Step 5: Report

`check_safety_box` ties the steps together and returns a `SafetyReport` with `missing`, `limits`, `exceeded`, an `ok` property, and a `format()` method that prints one finding per line and ends with `PASS` or `FAIL`.

### Verify

The tests in `tests/test_main.py` use fixture strings and check behaviour, not numbers from a file:

- a complete box passes with nothing missing and nothing exceeded;
- deleting one row names exactly that row; a `…` cell counts as missing; a template-like box lists all nine rows in order;
- each unit converts to the right key and factor (`250 g` → 0.25 kg, `36 km/h` → 10 m/s, `1,5 kg` → 1.5 kg);
- `3 m/s` does not create a height, `500 mA` does not create a height or a mass;
- values exactly at the envelope pass and 1 % above fail;
- 2 kg at 2.0 m (39.2 J) passes but 2 kg at 2.1 m trips both the height and the derived energy;
- a custom envelope is respected.

Run `pytest curriculum/phase-00/lesson-13-safety-for-home-experiments` and make every test pass.

### Compare

There is no established library for this job, so the comparison is against a human review. Give your safety box for the drop test to someone else with the template open and ask them to mark every row they would reject. Then run the checker. The checker will catch every structural problem the reviewer catches (missing rows, numbers over the envelope) and none of the substantive ones (a control that does not actually control the hazard, a stop criterion that needs judgment). That gap is the subject of part 5.

## 4. Experiment at home, or investigate real data

**Mode: design-build (a document).** There is no apparatus in this lesson. The product is the safety box you will use in `P00.15`.

A safety box is easier to write with the floor plan in front of you. Sketch yours; the drop test's looks like this from above:

```
   ┌──────────────────────────────────────────────┐
   │  wall with height marks                      │
   │     ●  release point (ball drops here)       │
   │   ╭────╮                                     │
   │  (  mat )  ← 1 m exclusion circle: nobody    │
   │   ╰────╯     inside it while a ball is held  │
   │                                              │
   │             ▣ phone on stand, 1.5 m away     │
   │                                              │
   │  ✕ experimenter stands here, outside circle  │
   │  door closed ──────────────────────────────  │
   └──────────────────────────────────────────────┘
```

Every control in the box should be visible in this picture (the mat, the circle, the distance to the phone, the closed door) or be something a person wears or does (glasses, the "clear" call). If a control cannot be drawn or acted, it is probably a wish, not a control.

1. Copy the table from `projects/_template/safety.md` into a new file `my-drop-test-safety.md`.
2. Fill the **Energy inventory** first, on paper: the ball's mass, the greatest height you will use, m g h, the impact speed √(2 g h), the phone's battery voltage. Write each number with its unit.
3. Fill the nine rows. For Hazards, go down the table in §2.1 and write one line per energy form that is present, plus eye and pinch. For Controls, write them in hierarchy order. For Stop criteria, write at least three, each observable and binary, and include "I am tired or rushed".
4. Run `python curriculum/phase-00/lesson-13-safety-for-home-experiments/starter/main.py my-drop-test-safety.md`. Fix every `MISSING` and every `OVER ENVELOPE` line until the report says `PASS`.
5. Keep the file. `projects/p00-first-experiment/safety.md` (written for the course) is the reference; compare yours against it when you reach `P00.15` and note what you missed.

## 5. Validate and improve

Run the checker on the raw template `projects/_template/safety.md`. It should fail with every row reported missing, because every cell is an instruction rather than content. Now run it on your drop-test box. It passes. Does that make the experiment safe?

No. Compare what the two checks can see. The script can see *structure* (all rows present, no placeholders) and *magnitude* (numbers inside the envelope). It cannot see *validity*. A Controls cell that says "be careful" passes the checker and controls nothing. A Stop criteria cell that says "if it looks dangerous" passes the checker and is not a criterion. An Energy limits cell that says "a small ball from head height" passes the checker (no numbers, so nothing exceeds) while hiding the fact that you have not computed anything.

Improve the box, not the checker, to close that gap: every control names the hazard it controls ("mat: falling mass; glasses: bounce to eye"), every stop criterion could be applied by a stranger, and Energy limits always carries numbers with units so the envelope check has something to bite on. Then improve the checker in one honest way: make it report when the Energy limits cell contains *no* parsable number at all, since a limit without a number is not a limit. (That is Problem 4.)

Where does this approach fail? The envelope is tuned for an adult alone in an ordinary room. It is not tuned for a child, for a cramped space, or for a shared flat where someone can walk in. Supervision and Location rows exist to carry exactly that information, and no number in the envelope substitutes for them.

## 6. Defend your work

Submit: your energy inventory table, your completed `my-drop-test-safety.md`, the checker's `PASS` output, your passing test run, and the problem set. A short note (five sentences) on what the checker could not catch in your box and how you closed that gap by writing.

### Adversarial self-checks

- What if the mass doubles? The energy doubles; does any control change, and does the envelope still hold? What if the height doubles instead?
- What if the ball is replaced by a steel bearing of the same diameter (about eight times the mass)? Which row changes first?
- Which assumption fails first in the real world: the adult alone, the clear floor, or the stated height?
- What observation would show that your stop criteria were not binary: that two people would disagree on whether one was met?
- Does the checker's verdict survive rewording the Energy limits cell ("1.8 metres" instead of "1.8 m")? If not, is that a bug in the parser or a rule for how the cell must be written?
- Which decision would change if your estimate of the drop height were off by a factor of two?

## Check yourself

1. Define *hazard* and *control* in one sentence each.
2. A 0.5 kg bag of rice is lifted to 1.5 m. What energy does it store, in joules, and is that inside the envelope?
3. Why does the envelope limit voltage *and* current separately?
4. Give the formula for the tip speed of a rotor of radius r turning at n revolutions per minute.
5. Which of these is a valid stop criterion, and why: (a) "if the setup seems unstable", (b) "if the clamp moves by any visible amount", (c) "if someone enters the room"?
6. Put these controls in hierarchy order: safety glasses, using a lighter ball, a landing mat, an exclusion circle.
7. The checker passes a box whose Controls cell reads "take care". What kind of error is that: a false pass, a false fail, or neither?
8. Name two things about a lithium-ion cell that are stop criteria.
