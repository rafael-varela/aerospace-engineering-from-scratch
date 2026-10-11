# Aerospace Engineering from Scratch — Consolidated Curriculum Recommendations

**Status:** Proposal for revising `ROADMAP.md` (Draft 1, 2026-10-10)  
**Purpose:** Turn an excellent topic roadmap into a rigorous, self-directed aerospace engineering education centered on affordable, safe, home-accessible experimentation.  
**Scope:** Recommendations only. The existing 27-phase, 448-lesson roadmap is the baseline, not replaced here.

---

## 1. Executive assessment

The original curriculum has unusually strong breadth across mathematics, mechanics, fluids, thermodynamics, aerodynamics, propulsion, structures, aircraft performance, flight dynamics, GNC, astrodynamics, spacecraft systems, avionics, experiments, systems engineering, and design. Its **computational-first philosophy**—derive a method, implement it independently, then compare with professional tools—is worth retaining.

Its main problems are educational rather than topical:

1. **Coverage is not mastery.** Completing a numerical implementation does not establish physical understanding or engineering judgment.
2. **Instruction-hour estimates are too optimistic for beginners.** The plan has 448 lessons and roughly 1,113 estimated lesson hours before its separate 200+ capstone allocation. Study, exercises, experiment setup, rework, and reporting require separate estimation.
3. **Physical experimentation arrives too late.** Hands-on work needs to appear within the foundations and recur throughout, not mainly in the experimental-methods phase and final capstones.
4. **Independent learning needs stronger assessment.** Clear rubrics, validation expectations, and opportunities for external review are necessary.
5. **Several supporting foundations are thin.** Precalculus, general chemistry, technical CAD/drawings, early electronics labs, and early engineering safety/standards need reinforcement.
6. **Specialized paths have prerequisite gaps.** A direct check of their listed lesson dependencies finds 19 missing prerequisite links across five specialist/discipline paths.
7. **“Degree-equivalent” needs careful wording.** Comparable topic coverage is not the same as independently demonstrated competence, institutional accreditation, or a formal credential.

**Recommendation:** Preserve most of the existing technical topics. Strengthen the *learning architecture*, make safe home experiments central, add integrated projects, repair dependencies, and calibrate workload using pilot learners.

## 2. Product and educational thesis

> **Learn aerospace engineering from first principles. Design it on your computer. Build it at home. Test it against reality.**

The curriculum should enable a motivated self-learner—with a computer, basic tools, a modest budget, and no university lab—to progressively develop theoretical knowledge, computational capability, physical intuition, experimental skill, and defensible engineering judgment.

The governing principle is:

> **Every major concept should connect to something a learner can calculate, simulate, measure, build, or observe.**

A home experiment is **not** a craft tutorial. It must require the learner to make a prediction, carry out a safe experiment, record data, analyze uncertainty, compare prediction and reality, investigate discrepancies, and defend a conclusion.

**Do not weaken academic depth in pursuit of accessibility.** Lower the equipment and institutional barriers, not the standards of reasoning, documentation, and evidence.

## 3. Preserve the roadmap's existing strengths

Keep these aspects of Draft 1:

- The common core (Phases 00–09), domain-specific phases (10–23), design phases (24–25), and capstones (26).
- The combination of aeronautics and astronautics, with narrower learning paths for specialists.
- First-principles mathematics, mechanics, and physical modeling.
- Building algorithms independently before relying on professional software, where doing so provides substantial educational value.
- Examples such as implementing a panel method before using XFOIL, an orbit propagator before GMAT, and a basic FEM solver before industrial FEM tools.
- Engineering ethics, reliability, safety, verification, technical communication, and design reviews.
- Phase exams, mini-projects, explicit cross-phase prerequisites, and optional advanced `★` material.

**Important nuance:** “From scratch” should mean understanding assumptions, governing equations, numerical methods, and validation—not rebuilding every production tool or algorithm regardless of educational value.

## 4. Rebuild the lesson model around engineering evidence

### 4.1 New default lesson sequence

Replace the current five-part arc (`The Problem → The Concept → Build It → Use It → Check Yourself`) with a six-part teaching model:

1. **The engineering problem.** Present a real physical question, design trade, or observable phenomenon, along with constraints and success criteria.
2. **Understand and derive.** Explain the physics; derive governing equations; identify units, assumptions, boundary conditions, and applicability; work representative analytic problems.
3. **Build the computational model.** Implement an educational model where useful; test numerical correctness; characterize error and failure modes; compare to established tools.
4. **Experiment at home / investigate real data.** Construct and measure something safe and accessible when appropriate. Otherwise use an authentic dataset, observational activity, or well-defined simulation study.
5. **Validate and improve.** Compare analytical, computational, and experimental results; quantify uncertainty; explain discrepancies; iterate.
6. **Defend your work.** Submit calculations, source code, raw data where applicable, plots, engineering interpretation, limitations, and a concise technical report.

Not every lesson needs hardware. **Every major technical phase needs applied evidence**, and every physical experiment needs a theoretical basis.

### 4.2 Standard lesson metadata

Each developed lesson should specify:

| Field | Required content |
|---|---|
| Learning outcomes | Observable skills, not vague topics |
| Prerequisites | Prior lesson IDs plus math/physical skill checks |
| Engineering question | The prediction or decision to be made |
| Assumptions and limits | Model conditions and known failure cases |
| Analytical work | Derivation and calculation exercises |
| Computational deliverable | Code, tests, numerical accuracy checks where meaningful |
| Practical activity | Home experiment, dataset exercise, observation, or simulation |
| Equipment | Required vs optional tools, substitutes, approximate budget tier |
| Safety | Hazard identification, risk controls, and stop conditions |
| Evidence | Data, plots, uncertainty, photos/videos only when helpful, report |
| Assessment | Rubric, worked examples or test cases, pass/fail expectations |
| References | Textbook chapters, lectures, papers, standards, and documentation |
| Extensions | Optional advanced or specialist challenges |

### 4.3 Meaningful distinctions between modes

Use consistent labels so students know what experience they have completed:

- **Physical experiment:** Student measures a real system using sensors, instruments, or observable motion.
- **Computational experiment:** Student studies a mathematical model, numerical method, or simulation, including sensitivity and validation.
- **Data-based investigation:** Student analyzes authentic measurements collected by others.
- **Design/build activity:** Student fabricates, assembles, or configures a physical artifact.

These modes can complement but should not be presented as perfectly interchangeable. A simulation does not expose a learner to all the calibration, tolerance, noise, wiring, and physical test issues of real hardware.

## 5. Define real mastery, not lesson completion

### 5.1 Four independently assessed competencies

Every major technical lesson should assess:

1. **Mathematical and physical reasoning:** Derive and explain governing relationships and assumptions.
2. **Analytical problem-solving:** Solve unfamiliar engineering problems, check units, and make order-of-magnitude estimates.
3. **Computational competence:** Implement and test a model; assess convergence, numerical error, and sensitivity.
4. **Engineering interpretation:** Decide what results mean, when the model is valid, what could fail, and how to validate it.

A passing code test does **not** automatically pass the lesson. A learner who produces plausible results without understanding the physics should not be certified as proficient.

### 5.2 Recommended assessment layers

- **Quick checks:** Low-stakes questions testing terminology, units, limiting cases, and intuition.
- **Problem sets:** Novel calculation and derivation tasks with worked reference solutions or structured hints.
- **Implementation tests:** Correctness and numerical-regression tests, including edge cases and convergence checks.
- **Experiment evaluations:** Calibration, raw measurements, uncertainty, reproducibility, and quality of conclusions.
- **Phase exam:** Timed or otherwise independent assessment that mixes conceptual and technical questions.
- **Phase mini-project:** Integration of multiple lessons into a documented engineering investigation.
- **Major design reviews:** Written and oral/video defense with a checklist resembling preliminary and critical design review expectations.
- **External challenge:** Optional peers, mentors, clubs, or independent reviewers assess the learner's reasoning.

### 5.3 Rubric for a representative technical project

Suggested grading dimensions (weights are **proposed**, not taken from Draft 1):

| Dimension | Illustrative weight |
|---|---:|
| Correct physics, derivation, and assumptions | 25% |
| Engineering calculation and numerical correctness | 20% |
| Experiment design, data quality, and uncertainty | 20% |
| Verification, validation, and investigation of discrepancies | 20% |
| Communication, reproducibility, and safety documentation | 15% |

For simulation-only projects, substitute appropriate numerical or data-based validation for the physical-experiment dimension; identify this explicitly. Use minimum competency thresholds in critical categories, rather than allowing presentation quality to compensate for unsafe methods or incorrect physics.

### 5.4 Require adversarial self-checks

Each major lesson should ask questions such as:

- What if the load doubles? What if velocity approaches zero? What are the limiting cases?
- Which assumptions would fail first in the real world?
- What measurement could falsify the model?
- Does the result survive a change in mesh size, time step, solver, or input parameters?
- What decision would change if measurement uncertainty were twice as large?
- What independent evidence would make this result credible?

## 6. Build a home aerospace laboratory progressively

**Accessibility is a hard requirement:** a student should never need a university wind tunnel, machine shop, costly instrumentation, specialized software license, or large workshop to complete the *core* path.

### Level 1 — Basic engineering lab

**Typical equipment:** ruler, measuring tape, digital scale, weights, common hand tools, safety glasses, small clamps, low-voltage breadboard equipment, inexpensive multimeter, basic microcontroller and sensors. A computer and phone camera can substitute for dedicated motion-capture equipment.

**Activities:** mechanics, beam deflection, pendulum dynamics, motion measurement, basic circuits, unit calibration, data logging, and feedback control.

### Level 2 — Prototyping lab

**Optional equipment:** soldering tools with appropriate ventilation, small servos/motors, low-power electronic supplies, IMU, load cell, simple fixtures, and optional 3D printer (or pre-cut/cardboard/laser-cut alternatives).

**Activities:** materials testing, instrumented structures, actuators, embedded systems, basic controls, simple aerodynamic models, repeatable fabrication, and measurement calibration.

### Level 3 — Aerospace experimentation

**Optional equipment:** guarded, low-speed bench wind tunnel, instrumented glider models, flight-data logger, desktop attitude-control platform, CubeSat bench prototype.

**Activities:** lift/drag comparisons, structural validation, flight-performance measurement, sensor fusion, closed-loop control, subsystems testing, and integrated engineering trade studies.

### Equipment and budget requirements

- **Buy incrementally:** Hardware is introduced at the point of use, not as an upfront purchase list.
- **Give three choices per project:** minimum-cost method, recommended method, and enhanced method.
- **Offer substitutions:** household objects and phone sensors where scientifically appropriate; community makerspaces or clubs for optional advanced work.
- **No mandatory 3D printer:** Patterns, hand-built parts, cardboard/foam, off-the-shelf components, or local fabrication can be alternatives.
- **Publish and maintain bills of materials (BOMs):** Quantity, specification, safe substitutions, estimated price ranges, sources, and tools required. Validate prices periodically rather than hardcoding them as permanent facts.
- **Support limited space and low resources:** Small-footprint apparatuses, low-noise activities, reuse of materials, and accessible storage.
- **Support different accessibility needs:** Alternate measurement methods and equivalent learning outcomes when reasonable.
- **Require a calibration plan:** A cheap instrument with known uncertainty is more useful educationally than an uncalibrated expensive device.

## 7. Make safety a design constraint for every hands-on project

Every lab guide must include an explicit safety box **before** its procedure: hazards, controls, supervision requirements, safe operating envelope, and criteria for stopping the experiment.

Core home-lab activities should favor low energies, low voltages, non-toxic materials, small masses, and guarded moving parts. For gliders and RC aircraft, require safe outdoor locations, respect for bystanders, preflight checks, and compliance with applicable local regulations. For electronics, include battery handling, electrical protection, mechanical guarding, and thermal precautions.

**Do not require at-home combustion, rocket motor firing, high-pressure vessels, dangerous propellants, uncontrolled high-speed rotors, or other high-consequence experiments.**

For rocket propulsion and high-energy aerospace subjects:

- At home: use derivations, nozzle calculations, airflow demonstrations, instrumented low-energy analogs, authentic datasets, simulations, and design reviews.
- At qualified facilities/clubs: allow supervised, legally compliant activities with trained personnel, appropriate infrastructure, and site-specific safety protocols.
- Always offer a rigorous non-combustion learning route for core credit.

A project is not acceptable just because the materials are available; foreseeable risks must be manageable in the intended setting.

## 8. Progressive home-project portfolio

These are **proposed additions and integrations**, not existing requirements in Draft 1. Each should have clear specifications, analysis questions, equipment tiers, safety requirements, a measurement protocol, and evaluation rubric.

| Project | Phases to integrate | Build / investigate | Required engineering evidence |
|---|---|---|---|
| **01. Cantilever beam** | 01, 04, 05 | Load a small beam and measure deflection/strain | Predicted vs measured displacement, load curve, uncertainty, error explanation |
| **02. Dynamic-system identification** | 02, 06, 09 | Build a pendulum or spring–mass apparatus; log motion | Equations of motion, estimated damping, time-series fit, sensitivity |
| **03. Heat-transfer investigation** | 04, 07, 22 | Measure heating/cooling of a safe small object | Lumped vs distributed model assumptions, temperature curves, fitted parameters |
| **04. Simple fluid experiment** | 03, 08, 22 | Observe flow/pressure with a safe low-energy setup | Reynolds/similarity reasoning, calibrated measurement, comparison with model |
| **05. Low-speed wind tunnel** | 08, 10, 22 | Build a guarded small tunnel or use equivalent accessible setup | Flow characterization, blockage discussion, repeatability, lift/drag comparison |
| **06. Airfoil geometry and test** | 10, 11 | Generate/construct different airfoil sections and test safely | Theory vs computation vs experiment; stall/viscosity limitations |
| **07. Instrumented glider** | 10, 13, 16, 17, 22 | Design, build, and evaluate a small unpowered glider | CG, stability, estimated glide performance, structures, flight-test data |
| **08. Materials and structure comparison** | 05, 12, 13, 22 | Test low-load coupons, simple joints, or built-up beams | Stiffness, failure-mode reasoning, manufacturing variability, safety margins |
| **09. Feedback-control demonstrator** | 06, 09, 18, 21 | Control a low-energy motor/servo/rotating system | Controller model, measured transient response, robustness and limits |
| **10. Satellite attitude demonstrator** | 06, 09, 18, 20, 21 | Desktop platform with sensor and guarded low-energy reaction wheel | Attitude estimation, torque budget, disturbance model, tracking error |
| **11. CubeSat bench prototype** | 19, 20, 21, 25 | Non-flight model of power, telemetry, sensors, software | Power budget, fault cases, environmental assumptions, interface verification |
| **12. Autonomous flight in simulation** | 09, 17, 18, 21 | Closed-loop aircraft/UAV simulation with optional safe outdoor hardware extension | Stability, navigation errors, mode transitions, failure injection |
| **13. Propulsion modeling and safe analogs** | 07, 08, 11, 14, 15, 22 | Cycle/nozzle models, low-energy flow demonstrators, published test data | Thrust/efficiency predictions, thermodynamic assumptions, validation |
| **14. Mission systems trade study** | 19, 20, 23, 25 | Design a constrained spacecraft mission using existing tools and data | Requirements traceability, mass/power/communications budgets, risks |

**Design rule:** The point is not to collect hardware projects. Every project must lead to a testable prediction and a defensible engineering conclusion.

### Example: minimum project package for the beam experiment

- **Question:** How accurately does linear beam theory predict small deflections under a known end load?
- **Theory:** Free-body diagram, material/section assumptions, moment of inertia, beam-deflection derivation.
- **Model:** Unit-tested Python calculator; sensitivity to beam length, thickness, stiffness, and load.
- **Build:** Safely clamped beam, known small masses or force, ruler or phone-based measurement.
- **Procedure:** Measure the unloaded baseline; take repeated measurements across several low loads; avoid failure or permanent deformation.
- **Analysis:** Estimate measurement uncertainty; graph predicted and observed deflections; test proportionality and limiting cases.
- **Deliverables:** Assumptions sheet, calculations, code, raw data, plot, error analysis, report, and suggested redesign.
- **Check:** Learner must explain why discrepancies arise; merely matching one data point does not pass.

## 9. Add two evolving flagship projects (spiral curriculum)

Instead of isolating each discipline, maintain two integrated systems that students revisit as their knowledge increases.

### Track A — Progressive aircraft engineering

1. **Foundations:** Choose a simple glider configuration; establish geometry, mass, coordinate frames, basic loads, and requirements.
2. **Aerodynamics:** Estimate lift and drag; compare analytical approximations, numerical models, and a simple physical test.
3. **Structures and materials:** Size the wing and spar; measure beam deflection; account for tolerances and mass growth.
4. **Performance and flight dynamics:** Model glide range, trim, center-of-gravity sensitivity, and static/dynamic stability.
5. **Electronics and control:** Add a safe lightweight data logger; optionally build a simulation-tested controller before any real flight.
6. **Verification:** Conduct controlled, legally compliant low-risk tests; reconcile results and redesign.
7. **Final defense:** Report requirements, models, assumptions, test evidence, decisions, uncertainties, and limitations.

The starting physical object can be inexpensive foam/cardboard. Sophistication comes from the engineering process, not the cost of the aircraft.

### Track B — Progressive spacecraft systems engineering

1. **Foundations:** Define a notional mission; calculate orbit, geometry, basic mass, and environmental assumptions.
2. **Astrodynamics:** Propagate orbits; investigate perturbations and ground contacts; validate against professional software/data.
3. **Subsystem budgeting:** Size power, thermal, attitude-control, communication, and onboard computing needs.
4. **Physical demonstrator:** Build a desktop sensor/power/telemetry prototype; add attitude-control demonstration if safely practical.
5. **Integration:** Handle interfaces, fault cases, time synchronization, operational modes, and telemetry.
6. **Verification:** Exercise subsystem tests and reconcile predicted and measured performance.
7. **Final defense:** Produce a mission-design package and a prototype verification report.

**Choice architecture:** Students can select one flagship after the common foundation, do both for broader coverage, or complete equivalent supervised/community alternatives. This should not force everyone to purchase two large hardware kits.

## 10. Shift practical training earlier without disrupting conceptual prerequisites

The original Phase 22 is valuable, but it should be the **advanced integration and formal experimental methods phase**, not the learner's introduction to experimentation.

Recommended changes:

- **P00–P04:** Measurement units, calibration, uncertainty, safe experiment setup, first beam/pendulum measurements, reproducible datasets.
- **P05–P09:** Instrumented mechanics, thermal/flow measurement, circuits, system identification, simple feedback control.
- **P10–P17:** Low-speed aero testing, glider design, structures, performance and stability measurements.
- **P18–P21:** Sensors, data fusion, embedded systems, simulated autonomous flight, desktop attitude control and spacecraft subsystems.
- **P22–P23:** Formal experiment design, statistical methods, test plans, quality control, hazard analyses, verification plans, and technical reviews.
- **P24–P26:** Fully integrated design and build/test capstones with documented validation.

Begin requirements, margins, reliability thinking, risk assessment, and engineering design reviews in early mini-projects; formalize them in Phase 23.

## 11. Reinforce missing or underrepresented foundations

| Topic | Recommended treatment | Why |
|---|---|---|
| Precalculus / trigonometry | Optional placement assessment plus bridging lessons before P01 | A true beginner should not silently be assumed to know algebra, functions, trigonometry, and geometry |
| General chemistry | Add foundational modules or prerequisites before materials/combustion | Bonding, reactions, thermochemistry, and material behavior need systematic foundations |
| CAD and technical drawings | Dedicated early-to-intermediate sequence integrated with 05, 12, 13, 24 | Parts must be dimensioned, toleranced, assembled, fabricated, and inspected |
| GD&T and tolerance analysis | Introduce practical basics, then deepen in manufacturing/design | Design-for-manufacture is not just computational sizing |
| Introductory electronics lab | Start in foundational physics, build toward P21 | Students need wiring, measurement, calibration, power, and sensor intuition |
| Safety procedures and standards | Begin in P00; recur in every physical project | Safety and compliance cannot be deferred until late phases |
| Systems-engineering discipline | Early requirements and verification exercises; formal P23 | Real projects require traceability and trade-offs from the outset |

Do not automatically add an entire new phase for every gap: use placement modules and integrations when they preserve a manageable path.

## 12. Correct workload estimates and learning-time claims

### Current baseline

- 27 phases, 448 listed lessons.
- About **1,113 lesson hours in Phases 00–25** plus **200+ capstone hours** in Phase 26, based on the roadmap table.
- The roadmap states an estimated **2,500–3,000 total hours** after additional problem sets/exams. This estimate has not been validated with actual learners.
- A simple average for the 440 non-capstone lessons is about **2.5 lesson hours each**, which is not credible as a complete mastery-time measure for many major subjects.

### Recommendation

Report separate estimates for:

1. Reading/video/instruction.
2. Worked derivations and problem sets.
3. Programming, debugging, and numerical validation.
4. Equipment preparation and physical experimentation.
5. Data analysis and written reporting.
6. Exams, remediation, and design-review feedback.
7. Major design/build capstones.

Provide ranges for learners with different starting knowledge. Clearly distinguish an experienced programmer's time for Python-based assignments from a novice's time to acquire physical and mathematical intuition.

**Do not replace the current estimates with a different invented precise number.** Pilot representative lessons (especially calculus, numerical solvers, aerodynamics, CFD, structural analysis, labs) with real beginners and intermediates, measure completion and failure rates, and revise the forecast based on evidence.

The standard should be demonstrated competence, not finishing within a prescribed number of hours.

## 13. Fix prerequisite and path integrity

The listed cross-phase lesson prerequisites were checked against the five non-full learning paths. **Nineteen prerequisite links are not included in the corresponding path.** The following list is the direct-link audit of Draft 1 (it is not a full audit of implicit same-phase or conceptual prerequisites):

| Path | Lesson in path | Required lesson missing from path |
|---|---|---|
| Aeronautical Engineer | P21.08 | P20.15 |
| Aeronautical Engineer | P21.14 | P20.03 |
| Aeronautical Engineer | P22.10 | P15.15 |
| Astronautical Engineer | P13.17 | P10.10 |
| Astronautical Engineer | P13.18 | P10.19 |
| Astronautical Engineer | P18.20 | P17.16 |
| Astronautical Engineer | P21.09 | P17.02 |
| Astronautical Engineer | P21.13 | P17.12 |
| Astronautical Engineer | P22.05 | P10.01 |
| Astronautical Engineer | P22.12 | P16.14 |
| Propulsion Specialist | P14.14 | P10.17 |
| Propulsion Specialist | P22.12 | P16.14 |
| GNC & Autonomy Specialist | P17.03 | P10.12 |
| GNC & Autonomy Specialist | P17.06 | P10.11 |
| GNC & Autonomy Specialist | P21.08 | P20.15 |
| GNC & Autonomy Specialist | P21.14 | P20.03 |
| Structures & Materials Specialist | P13.18 | P10.19 |
| Structures & Materials Specialist | P22.10 | P15.15 |
| Structures & Materials Specialist | P22.12 | P16.14 |

Also resolve these inconsistencies:

- **P24.13** requires **P03.16** (Multidisciplinary Design Optimization), but P03.16 is marked optional (`★`) even in degree-equivalent paths. Decide whether to make the foundation mandatory, change the later requirement, or introduce a non-advanced substitute.
- **P26** requires at least one *design* capstone and one *build/test* capstone. The **Aeronautical Engineer** path lists `01 or 05, plus 07`, which could lead to `05 + 07` (two build/test projects, no design capstone). Make the allowed combinations satisfy the rule.
- Check the same design/build requirement for every path and whether the selected project is actually supported by included lessons.
- Audit **same-phase** sequencing too: the `Needs` column intentionally lists only external phase dependencies, so missing explicit links do not mean a lesson is truly independent.
- Check all `★` topics against downstream mandatory requirements.

### Definition of done for learning paths

A learning path must include the transitive closure of all **mandatory** prerequisites, or explicitly identify a legitimate equivalent bridging module. Validate this automatically from machine-readable lesson metadata, and run the check whenever the roadmap changes.

## 14. Make academic standards and the benchmark claims precise

Replace broad “degree-equivalent” assertions with three separate categories:

1. **Curriculum scope comparable:** Major subject matter maps to a relevant program benchmark.
2. **Demonstrated technical competence comparable:** A student passes independent, demanding assessments and completes validated projects comparable in complexity to undergraduate engineering work. This is a goal to demonstrate, not a status automatically earned by completing pages.
3. **Credential equivalence:** A self-directed roadmap does not grant an accredited degree, professional authorization, or professional engineering license.

Check the exact **current ABET accreditation criteria and program-specific aerospace criteria** when preparing public claims; do not rely indefinitely on a historical topic checklist. Maintain a versioned mapping to ABET student outcomes and verified university programs. A syllabus comparison is evidence of topic coverage, not accreditation.

Because self-study lacks inherent external oversight, strengthen:

- Teamwork through optional cohorts, group projects, reviews, and club collaboration.
- Communication through formal reports, public design logs, and oral/video defenses.
- Experimentation through measured work and authentic data, not just simulations.
- Professional responsibility through repeated hazard, ethics, regulations, and failure-analysis exercises.
- Independent inquiry through paper reading, reproduction of published results, and an open research capstone.

## 15. Reference and instructional-content standards

For every phase, build a *verified reference map* that connects lessons to actual chapters/sections rather than listing entire textbooks generically. **All required readings must be free to access; paid textbooks may be recommended as optional enrichment only.**

Potential sources already contemplated in Draft 1 include MIT OpenCourseWare, OpenStax, NASA technical reports and handbooks, and Leishman's *Introduction to Aerospace Flight Vehicles*. Add appropriate university texts, professional standards, peer-reviewed papers, open datasets, and technical documentation where useful. Recommend respected paid textbooks alongside the free resources when they offer genuinely valuable depth or another pedagogical approach, but never make purchasing them a requirement.

For each reference, document:

- Exact title, author/organization, edition/version, and section.
- Why it is included, what it teaches, and what the learner is expected to extract.
- Access label (**required and free** or **optional and paid**), approximate price or access options if applicable, and licensing status (do not assume permission to reproduce content).
- Whether it supports theory, computation, laboratory practice, or professional engineering.
- Last verification date for links, standards, and software instructions.

Solutions should provide enough worked reasoning for self-study without making rote answer matching the default. Consider hints, partial derivations, separate instructor keys, randomized parameters, and independent validation cases.

## 16. Delivery model for self-learners

For a usable website/repository or course, provide:

- **Diagnostic entry test:** Determines whether the student needs math, programming, or science bridge material.
- **Prerequisite-aware navigation:** Unlocks or recommends lessons based on competencies, not just phase numbers.
- **Visible deliverables:** Code, measurements, design artifacts, and reports tied to each lesson.
- **Reproducible starter repositories:** Unit-tested examples, standardized datasets, plotting templates, and report scaffolds.
- **Lab manuals:** Illustrated assembly plans, material substitutions, safety instructions, calibrations, data sheets, and troubleshooting.
- **Equipment planning:** A cumulative inventory showing what can be reused in future experiments.
- **Mastery tracking:** Separate achievement indicators for theory, implementation, experimentation, and engineering interpretation.
- **Public project portfolio option:** Make it easy to publish design reviews, technical reports, and experiments without exposing personal information by default.
- **Community support:** Optional peers, engineering clubs, peer review, and mentor feedback, while allowing a learner to proceed independently.
- **Offline/low-bandwidth access:** Text-first lessons, downloadable worksheets, and reusable code/data packages.

The interface should make it obvious **what the student will know**, **what they will build**, **what evidence they must produce**, and **what equipment they need** *before* starting a lesson.

## 17. Implementation roadmap — by capabilities, not calendar dates

### Priority 0 — Make the existing plan internally correct

- [ ] Correct all 19 direct prerequisite gaps or explicitly document bridging equivalents.
- [ ] Fix the mandatory/optional dependency at P24.13 → P03.16.
- [ ] Make capstone combinations consistent with the required design + build/test pairing.
- [ ] Audit additional same-phase and transitive prerequisites.
- [ ] Rewrite degree-equivalence language; verify the latest accreditation and university benchmarks.
- [ ] Separate lesson time, assessed study time, laboratory time, and capstone time.

**Done when:** Every published path is prerequisite-complete, capstone-compliant, and accurately labeled.

### Priority 1 — Create the educational contract

- [ ] Define measurable competencies and the six-part lesson template.
- [ ] Create project and exam rubrics with minimum thresholds for physics and safety.
- [ ] Publish standard formats for engineering reports, data, code tests, and verification/validation.
- [ ] Build a diagnostic/bridging path for absolute beginners.
- [ ] Select and verify teaching references for the first phases.

**Done when:** Another author can write a lesson that consistently meets the same standard.

### Priority 2 — Prove the home-lab learning model

- [ ] Publish a Level 1 lab starter guide, substitutions, safety checklist, and equipment inventory.
- [ ] Develop and pilot the beam, pendulum, and early data-collection projects.
- [ ] Demonstrate an actual theory → model → experiment → validation → report workflow.
- [ ] Gather learner feedback on cost, confusion, measurement quality, and workload.
- [ ] Revise instructions and time estimates using the pilot results.

**Done when:** A beginner can safely reproduce early physical investigations with documented, assessable results and affordable resources.

### Priority 3 — Integrate engineering across phases

- [ ] Map home projects to Phases 00–25 at moments where prerequisites are satisfied.
- [ ] Introduce the progressive aircraft and spacecraft flagship projects.
- [ ] Connect test results and reusable code from earlier phases to later design work.
- [ ] Add a structured CAD/drawings/tolerances sequence, electronics practice, and chemistry bridges.
- [ ] Expand validation against professional tools and authoritative measurement datasets.

**Done when:** Students continually apply knowledge to evolving engineering systems and can explain the accuracy and limits of their predictions.

### Priority 4 — Validate outcomes at capstone level

- [ ] Require documented engineering requirements and verification plans.
- [ ] Require at least one substantial design capstone and one genuinely assessed build/test experience for paths that claim them.
- [ ] Offer safe, rigorous low-cost variants and supervised-club extensions.
- [ ] Require technical reports, uncertainty, failure analysis, and final design defenses.
- [ ] Pilot external evaluation and compare demonstrated skills to defined university-level learning outcomes.

**Done when:** Completion demonstrates independently reviewable engineering work, rather than merely a finished lesson checklist.

## 18. Acceptance criteria for the next version of ROADMAP.md

A revised curriculum is ready to publish when all of the following are true:

1. **Topic integrity:** The existing major aerospace foundations and disciplines remain covered; changes are intentional and mapped.
2. **Mastery:** Each lesson has observable outcomes and meaningful evaluation, beyond code execution.
3. **Accessibility:** Each required practical competency has a feasible home, data-based, or appropriately supervised path; limitations are disclosed.
4. **Experimentation:** Physical experiments begin in the foundation sequence and recur throughout relevant phases.
5. **Scientific method:** Major projects involve prediction, observation/data, uncertainty, comparison, and interpretation.
6. **Safety:** Every hands-on project includes risk controls and does not presume unsafe at-home aerospace testing.
7. **Reproducibility:** Calculations, code, data, test plans, BOMs, and results can be independently reviewed.
8. **Progression:** Aircraft and/or spacecraft flagship systems mature as the student learns new disciplines.
9. **Correct dependencies:** All paths have valid prerequisite closure and consistent capstone selection rules.
10. **Realistic effort:** Workload estimates distinguish activity types and carry caveats until validated.
11. **Verified sources:** Lessons connect to cited, versioned, reliable teaching materials and technical references.
12. **Honest claims:** The distinction between subject coverage, demonstrated competence, and formal credentials is explicit.

## 19. Guardrails — what not to do

- **Do not turn the course into a hardware shopping list.** The models, measurements, and engineering decisions are the educational value.
- **Do not require costly machinery or proprietary software.** Keep the core accessible with open tools and alternatives.
- **Do not substitute impressive-looking code for physics.** Learners must be able to justify models and their validity.
- **Do not claim simulation-only work is identical to experimentation.** Acknowledge the competencies each develops and misses.
- **Do not make the curriculum fully linear just because the table is linear.** Revisit ideas in progressively richer projects.
- **Do not add topics solely to inflate the lesson count.** Improve teaching depth and integrate overlooked foundations.
- **Do not demand unsafe home propulsion or high-energy demonstrations.** Supervised external options must remain optional.
- **Do not promise an accredited degree or a professional credential.** Demonstrate work and learning outcomes instead.
- **Do not invent precise completion times.** Measure them.

---

## 20. Decisions on the three open questions in ROADMAP.md

These are **settled curriculum requirements for Step 2**, not unresolved discussion items. They resolve the questions under “Open questions for step 2” in Draft 1 of `ROADMAP.md` while preserving its existing 27-phase structure.

### 20.1 Hardware — budgeted BOMs, alternatives, and honest experimental credit

**Decision: YES.** Every hands-on lesson and capstone, especially P26.05–P26.07, must publish a reusable, region-agnostic bill of materials (BOM) with three procurement approaches and an accessible computational/data-based fallback.

Each project package must provide:

1. **Three equipment tiers:** **minimum-cost**, **recommended**, and **enhanced**. Assign itemized local-currency/region-adjustable prices; record the price-check date. Target broad accessibility rather than rigid universal dollar ceilings. A student must not have to own a 3D printer, precision machine tools, or a university laboratory.
2. **An itemized BOM:** item, quantity, functional specification, cost range, safety rating/limits where relevant, vendor-neutral search term, suggested substitutes, optional/required flag, required tools, reusable inventory ID, and any consumables.
3. **A reuse map:** show which prior lessons supplied the same sensors, microcontroller, multimeter, fixtures, or materials; calculate *incremental* spend as well as total spend.
4. **Complete build documentation:** schematics, mechanical drawings or cutting patterns, wiring, firmware, test setup, calibration, safe handling, troubleshooting, and disposal.
5. **A computational fallback:** required models and validation datasets, reproducible scripts, parameter studies, sensitivity analyses, and a comparable technical report. If possible, also offer a **safe physical surrogate** (a small home experiment exploring part of the same physics).
6. **An evidence and safety contract:** prediction before measurement; raw data, calibration, uncertainty, comparison, limitations, and stop criteria. Identify any required local permissions, approved flying locations, club supervision, or adult supervision.

For the three existing capstones:

| Capstone | Affordable hands-on entry point | Recommended/enhanced route | Accessible fallback | Critical boundary |
|---|---|---|---|---|
| **P26.05 — Fixed-wing UAV/autopilot** | Build a small foam glider or control-surface rig; instrument a bench test using a microcontroller/IMU and servo | Add a compatible flight controller, telemetry, and supervised outdoor fixed-wing testing where legal | 6-DOF flight dynamics + software-in-the-loop, sensor logs, control-system hardware-in-the-loop where accessible | Simulation or bench-only work does **not** claim that an independently developed autopilot was successfully flight-tested |
| **P26.06 — Hybrid/solid motor** | Perform an instrumented, guarded **low-voltage electric-fan thrust experiment** to practice thrust measurement and data reduction | Actual rocket-motor firing only through an appropriately qualified and equipped supervised club/facility following local laws and club procedures | Motor performance and trajectory simulation, public thrust datasets, uncertainty propagation, and comparison with the safe electric-thrust experiment | **No home propellant fabrication, motor ignition, or improvised pressurized firing equipment as a curriculum requirement**; safe surrogates do **not** count as a rocket-motor fire |
| **P26.07 — Wind tunnel** | Make a small, guarded, low-speed tabletop airflow rig; calibrate airspeed and use a simple force balance | Improve flow conditioning, mounting, calibrated instrumentation, data acquisition, and repeatability | Panel/VLM/CFD experiments with published wind-tunnel datasets, convergence and model-error reports | Document fan/electrical hazards, rotating parts, debris, model-mount integrity, and the limits of low-cost measurements |

**Completion rule:** A simulation/data-only fallback earns explicit **modeling and simulation competencies**, not automatically **physical experimentation competencies**. To claim a full design + build-and-test completion, the student must perform **at least one safe, accessible physical build-and-measure project** meeting a common experiment rubric. This need not be motor firing or autonomous outdoor flight: a sufficiently rigorous home wind tunnel, instrumented structure, sensor/attitude-control rig, or approved equivalent can meet the physical requirement. Other capstone alternatives should be offered if a particular capstone is inaccessible or unsafe.

**Procurement rule:** Prices are inherently regional and change over time. BOM budgets must be verified during authoring rather than presented as immutable amounts. Equipment is acquired progressively, not through a mandatory expensive starter kit.

### 20.2 Exams — transparent self-study solutions, robust checks, and independent mastery

**Decision: use all three assessment formats; keep solutions in the same public repository but in a clearly separated `solutions/` tree.** Do **not** rely on a private “instructor-only” directory, which undermines accessibility for independent learners. A folder named `instructor/` in a public repository is not access control. Conversely, do not place worked answers directly below the exercise prompts.

For each lesson, provide:

- **Auto-graded checks:** conceptual multiple-choice or short-answer questions, dimensional reasoning, boundary cases, and parameterized numerical checks, each with feedback explaining the result.
- **Substantial problem sets:** derivations, unfamiliar engineering scenarios, independent calculations, numerical analysis, and judgment about model limitations.
- **Engineering programming tasks:** public unit tests plus unseen or alternate validation cases for deeper evaluation; emphasize residuals, limits, conservation laws, and convergence rather than exact output matching only.
- **Three-step assistance:** targeted hint → conceptual/method outline → complete worked solution with assumptions, units, derivation, alternatives, checks, and common mistakes.
- **Phase exams:** a substantial integrated test, with a suggested timed/closed-reference *self-assessment mode*, mixing analytic problems, conceptual explanations, computational checks, and interpretation.
- **Phase mini-project:** a separate rubric-graded artifact; pass the phase based on both exam evidence and the project, rather than solely quiz completion.

**Recommended repository layout:**

```text
curriculum/
  phase-XX/
    lesson-YY/
      lesson.md
      problems.md
      starter/
      tests/
solutions/
  phase-XX/
    lesson-YY/
      hints.md
      worked-solutions.md
assessments/
  phase-XX/
    exam.md
    grading-rubric.md
    parameterized-cases/
projects/
  <project-id>/
    bom.yaml
    safety.md
    procedure.md
    data/
    report-template.md
references/
  sources.yaml
  phase-XX.md
```

**User experience:** hide worked solutions by default in the website interface, with an explicit “Reveal hint / Reveal method / Reveal solution” control. Solutions must remain available without payment, registration, or a human instructor, including in an offline copy. Encourage learners to attempt first and record self-reported completion before revealing answers; do not enforce artificial time locks that punish self-learners.

**Integrity:** A public, self-graded assessment cannot provide the security or independent proof of competence of a proctored examination. Use new variants, parameter randomization, alternate numerical cases, lab evidence, critical-category minimum scores, and optional peer/mentor design defenses to make superficial answer copying less useful. Clearly distinguish **practice completed**, **self-assessed mastery**, and **externally reviewed mastery**; none alone confers an accredited credential. Do not rehost restricted third-party instructor answer keys or test banks; author original assessment content and comply with source licenses.

### 20.3 Textbook anchors — required free references and optional paid textbooks

**Decision: YES.** Every phase gets a small curated set of academically credible, **free-to-access required references**, with **paid textbooks offered only as optional recommendations**. Every referenced resource maps to **specific sections or chapters**, not to broad book homepages. This curation is a required deliverable of Step 2 and must not be treated as completed simply because the reference providers below have been identified.

**Source selection policy:**

1. Prefer **one free primary learning anchor** per phase and **one or two complementary free references** when they supply missing depth, problems, applications, or a different explanation. Additional specialist sources are allowed if the phase cannot honestly be covered by a tiny list. Curate paid titles separately as optional enrichment.
2. Identify **exact chapter/section/lecture numbers and edition or access date** for every lesson; distinguish **required/free**, **supplementary/free**, and **optional/paid** reading. Never label a paid resource essential.
3. Guarantee a **complete no-cost learning route** using legally accessible original sources (open textbooks, freely accessible course content, public technical reports, agency handbooks) plus original instructional material. Where a reliable free source does not provide sufficient depth, flag the gap and author original instruction; **never use a paid book to patch a required-learning gap**. Do not falsely claim that a reference teaches material it merely mentions.
4. Verify the **specific linked content**, depth, preconditions, language/accessibility, and availability of example problems. Note when lecture materials are incomplete or course-specific.
5. Record **license and reuse restrictions separately from reading cost**. Free-to-read does not automatically mean free to copy, redistribute, sell, or modify. Link to original material and write original explanations rather than copying protected chapters or restricted instructor resources.
6. Provide a machine-readable **reference manifest** (`references/sources.yaml`) with title, creator, version/date, canonical URL, license or usage terms, **access tier (required-free / supplementary-free / optional-paid)**, mapped lesson IDs/sections, verification date, and known caveats. For paid titles, add ISBN where relevant, a short explanation of the book's added value, and non-purchase options (library, used copy, legitimate rental) when known. Prices should be treated as estimates, not guarantees. Run link checking in CI and conduct periodic human content-quality reviews. A working HTTP response alone is not content verification.
7. Provide a **free-resource fallback** for broken links when available, and mark coverage gaps transparently.
8. **Paid books are never prerequisites.** No lesson, exercise, exam, project, or capstone may require purchasing a text or accessing paywalled material. The course must include enough original explanation, examples, problems, and reference material for a motivated learner to succeed entirely for free.
9. Present paid books in a clearly marked **“Optional textbooks for deeper study”** section after the free learning resources. Recommend them only where they add specific pedagogical value; avoid long generic shopping lists and undisclosed affiliate incentives.

**Verified starting catalog (access checked 2026-10-10; each phase/lesson mapping still requires Step 2 curation):**

| Source | Suitable initial use | Official access | Reuse caveat |
|---|---|---|---|
| MIT OpenCourseWare | Mathematics, physics, dynamics, controls, engineering design, CAD, and specialist lectures | https://ocw.mit.edu/ | Check the license and completeness for each course and embedded asset |
| OpenStax | Foundational calculus, calculus-based physics, chemistry, and statistics as applicable | https://openstax.org/subjects | Verify the specific textbook edition and rights for any reused material |
| J. Gordon Leishman, *Introduction to Aerospace Flight Vehicles* (Embry-Riddle) | Broad aerospace survey, fundamentals, aircraft systems, structures, aerodynamics, and design | https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/ | Site currently marks the book **CC BY-NC-ND**; link to it, do not assume commercial reuse, adaptation, or derivative republication is permitted |
| NASA Technical Reports Server | Aerospace technical reports, worked analyses, reference data, experiment methods | https://ntrs.nasa.gov/ | Inspect rights, technical relevance, and age for each individual report |
| NASA Systems Engineering Handbook | Requirements, design and verification processes, reviews, and systems practice | https://www.nasa.gov/reference/systems-engineering-handbook/ | Treat as engineering guidance, not a substitute for current applicable regulatory/technical standards |

**Pilot anchor examples to author in Step 2:**

- **P00.11 (Reading a Technical Paper):** teach how to read a NASA NTRS report; curate one freely available report and annotate its assumptions, methods, limitations, and evidence.
- **P04 (Physics Foundations):** use selected OpenStax *University Physics* sections appropriate to the individual lessons, adding aerospace-specific examples and the mathematics they require.
- **P23.02 (Requirements Engineering):** pair the NASA Systems Engineering Handbook's **Appendix C, “How to Write a Good Requirement”**, with original aerospace requirement-writing exercises.

**Acceptance test:** Before a phase is marked ready, a reviewer can pick any lesson and find a **free** accessible primary explanation, a relevant section or worked example, original assessed exercises, and a clear explanation of any reference coverage gap. The reviewer must be able to complete all required work **without buying or borrowing a paid textbook**. Paid recommendations, if included, are visibly optional and mapped to concrete learning benefits.

### 20.4 Required change to ROADMAP.md

Replace the three bullets under `## Open questions for step 2` with `## Decisions for step 2` summarizing these requirements:

- **Hardware:** Publish three budget-tier BOMs, reuse/substitution guidance, safety and build documents, a simulation or open-data fallback, and explicit distinctions between simulated and physically measured outcomes. Motor firing only in qualified supervised settings; accessible alternative physical projects preserve the overall build-and-test route.
- **Exams:** Auto-graded checks + original rigorous problem sets + phase exams + mini-projects, with hint/method/full-solution scaffolding. Keep full solutions in the public repo's separate `solutions/` directory, hidden by default in the interface but freely accessible. Never claim public unsupervised exams establish credential equivalence.
- **Textbook anchors:** Curate, verify, and map **free-to-access required** primary/supplementary references to exact lesson sections; provide **optional paid textbook recommendations** for additional depth, never as prerequisites or substitutes for core free instruction. Maintain a reference manifest with access tier, source, edition/date, license, link, and last-verification date. Seed free sources with MIT OCW, OpenStax, Leishman, NASA NTRS, and the NASA Systems Engineering Handbook; fill required coverage gaps with original instruction.

**Definition of done for Step 2:** Each phase has a verified reference map (complete free required reading plus clearly marked optional paid books where beneficial) and original assessment plan; every physical project has a safety-reviewed, cost-tiered build package and an honest accessible alternative; solution availability and achievement types are explicit to independent learners.

---

## Final direction

The key change is from **“448 aerospace topics to finish”** to **“an engineering education whose graduates can model, build, measure, validate, and defend technical decisions.”**

Preserve the academic ambition and numerical-first rigor of Draft 1. Make affordable, safe, progressive home experiments part of the curriculum's identity. Demand clear evidence of physical understanding and trustworthy results. Build enduring aircraft and spacecraft projects that turn isolated lessons into connected engineering competence.

**North star:** A learner confronted with an unfamiliar aerospace problem should be able to establish requirements, identify governing physics, make justified assumptions, solve and simulate, design a practical test, evaluate risk and uncertainty, improve a solution, and defend it under technical scrutiny.
