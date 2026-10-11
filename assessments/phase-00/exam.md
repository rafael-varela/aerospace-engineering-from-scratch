# Phase 00 Exam: Setup, Engineering Computing & First Measurements

Integrated assessment of the whole phase. Suggested self-assessment mode: closed reference, 2 hours chosen and recorded in advance, no code except in Part C. There is no enforced lock; the point is to find out what you can do unaided. Solutions and pass thresholds are in `solutions/phase-00/exam/`.

Report your result per competency: (1) mathematical and physical reasoning, (2) analytical problem solving, (3) computational competence, (4) engineering interpretation. Questions are tagged with the competencies they assess.

## Part A: Concepts

**1.** (1, 4) A colleague says their new altimeter model "matches the standard atmosphere exactly." Explain, in five sentences or fewer, why that sentence cannot be literally true, and what they should say instead. Use the words *model*, *assumption*, and *uncertainty*.

**2.** (1) The Mars Climate Orbiter was lost because one piece of software produced thruster impulse in pound-force seconds and another consumed it as newton-seconds. By what factor was the impulse wrong, and in which direction (too much or too little effect modeled)? What single engineering practice from this phase would have caught it, and at which interface?

**3.** (4) Distinguish, with one aerospace example each: a *physical experiment*, a *computational experiment*, and a *data-based investigation*. Which competency does each one build that the others do not?

## Part B: Analysis

**4.** (1, 2) You read a digital scale ten times with the same object: 52, 53, 52, 52, 53, 52, 52, 53, 52, 52 g. The display resolution is 1 g. Compute the Type A standard uncertainty of the mean and the Type B standard uncertainty from resolution. Which dominates? State the mean with its combined standard uncertainty and with an expanded uncertainty at k = 2, to the correct number of significant figures.

**5.** (1, 2) The fall time of an object from height h is t = √(2h/g). You measure h = 1.500 ± 0.005 m and take g = 9.81 ± 0.03 m/s². Derive the expression for the standard uncertainty of t by the law of propagation, then evaluate it. Which input contributes more, and by what ratio?

**6.** (2) Using only the ISA sea-level values (T₀ = 288.15 K, p₀ = 101325 Pa, lapse rate 6.5 K/km) and the hydrostatic relation for a linear temperature profile, estimate the pressure at 5 km. Then state, without computing, whether the density at 5 km is more or less than half the sea-level value, and why.

## Part C: Computation

**7.** (3) Write a Python function `drag_force(rho, v, cd, area)` that returns the aerodynamic drag ½ρv²C_D A, and three pytest tests for it: one tolerance-based check against a hand-computed value, one limiting case, and one that asserts a `ValueError` for a negative area. Then explain in two sentences why a test that uses `==` on the floating-point result would be a bad test.

## Part D: Interpretation

**8.** (4) A learner's free-fall experiment gives normalized errors of 0.3, 0.4, 2.8, and 3.1 at heights of 0.5, 1.0, 1.5, and 2.0 m (measured times longer than predicted at the two higher drops). Rank three candidate causes, say which single follow-up measurement would best discriminate between the top two, and say what result from that follow-up would support each one.
