# Reference map: Phase 00: Setup, Engineering Computing & First Measurements

Verified 2026-10-11 by loading the live pages; section numbers and titles are copied from the publishers' own chapter outlines or section pages. Every source here is **required-free** unless marked otherwise; sources are identified by an `@`-prefixed identifier and recorded in [sources.yaml](sources.yaml) with license, version, and caveats. **GAP** marks a lesson that no listed free source covers adequately: the lesson must carry original instruction there, and any "candidate" named is unverified until it is added to the manifest.

Conventions: OpenStax section URLs follow `https://openstax.org/books/<book-slug>/pages/<n>-<m>-<title-slug>`; a URL written as `.../pages/...` continues the last full URL above it. MIT OCW video URLs follow `https://ocw.mit.edu/courses/<course-slug>/resources/lecture-<n>-<slug>/`. OpenStax and MIT OCW are **CC BY-NC-SA 4.0**, Leishman is **CC BY-NC-ND 4.0**: the curriculum links to them and writes its own explanations; nothing is copied or adapted.


## 00-01 Dev Environment for Engineers

Sources: `@numpy-userguide`, `@pytest-getstarted`

- NumPy user guide — "Installation" (numpy-userguide). https://numpy.org/doc/stable/user/index.html
  Extract: install via pip/conda; verify import.
- pytest — "Get Started" > "Install pytest" (pytest-getstarted). https://docs.pytest.org/en/stable/getting-started.html
  Extract: `pytest --version` check as the smoke test for the environment.
- GAP (editor/git/virtual-env workflow): no single authoritative free engineering-specific source; use
  the tool vendors' own docs.

## 00-03 NumPy as an Engineering Calculator

Sources: `@numpy-userguide`

- NumPy user guide — "NumPy quickstart"; "NumPy: the absolute basics for beginners"; "NumPy
  fundamentals" > "Broadcasting", "Data types", "Copies and views". https://numpy.org/doc/stable/user/index.html
  Extract: ndarray creation, vectorized arithmetic, broadcasting rules, float64 default.
- NumPy user guide — "NumPy for MATLAB users". Same URL. Extract: translation table for students with MATLAB background.

## 00-04 Plotting That Engineers Trust

Sources: `@matplotlib-userguide`

- Matplotlib "Using Matplotlib" — "Quick start guide"; "Axes and subplots"; "Text"; "Colors";
  "Customizing Matplotlib with style sheets and rcParams". https://matplotlib.org/stable/users/index.html
  Extract: explicit Figure/Axes API, labeled axes with units, log scales, saving vector output.

## 00-05 Units & Dimensional Consistency (pint, Mars Climate Orbiter)

Sources: `@nasa-mco-mib`, `@pint-tutorial`

- pint — "Tutorial": "Initializing a Registry", "Defining a Quantity", "Converting to different units",
  "String parsing", "String formatting". https://pint.readthedocs.io/en/stable/getting/tutorial.html
  Extract: UnitRegistry, Quantity arithmetic, DimensionalityError as a unit-check.
- NASA, "Mars Climate Orbiter Mishap Investigation Board Phase I Report" (1999-11-10).
  https://llis.nasa.gov/llis_lib/pdf/1009464main1_0641-mr.pdf
  Extract: root cause — "failure to use metric units in the coding of a ground software file, 'Small
  Forces'" (lbf-s vs N-s); contributing causes (interface spec not followed, anomalies not escalated).

## 00-06 SciPy Tour

Sources: `@scipy-userguide`

- SciPy User Guide (v1.18.0) — sections "Integration (scipy.integrate)", "Optimization
  (scipy.optimize)", "Interpolation (scipy.interpolate)", "Linear Algebra (scipy.linalg)", "Sparse
  Arrays (scipy.sparse)", "Statistics (scipy.stats)", "Special Functions (scipy.special)".
  https://docs.scipy.org/doc/scipy/tutorial/index.html
  Extract: one worked call per submodule; map to later Phase 03 lessons.

## 00-07 Testing Engineering Code (pytest)

Sources: `@numpy-userguide`, `@pytest-getstarted`

- pytest — "Get Started": "Create your first test", "Assert that a certain exception is raised",
  "Compare floating-point values with pytest.approx", "Group multiple tests in a class".
  https://docs.pytest.org/en/stable/getting-started.html
  Extract: approx() for tolerance-based assertions; raising on bad units/inputs.
- NumPy — `numpy.testing` (assert_allclose) via the NumPy reference (secondary). https://numpy.org/doc/stable/reference/routines.testing.html

## 00-09 Symbolic Math with SymPy

Sources: `@numpy-userguide`, `@sympy-tutorial`

- SymPy 1.14.0 — "Introductory Tutorial": "Gotchas" (Symbols, Equals signs), "SymPy Features" >
  "Basic Operations", "Simplification", "Calculus", "Solvers", "Matrices".
  https://docs.sympy.org/latest/tutorials/intro-tutorial/index.html
  Extract: symbols vs Python names, Eq(), diff/integrate/series, lambdify to NumPy.

## 00-10 Writing Technical Reports

Sources: `@psu-style-for-students`, `@psu-writing-engineer`

- Joe Schall, "Style for Students Online" (Penn State, Dutton Institute; CC BY-NC-SA) — Chapter 4
  "Equations, Figures, and Tables"; Chapter 5 "Using Sources"; Chapter 6 "Writing Documents for Classes".
  https://courses.ems.psu.edu/styleforstudents/  (free: "available for free use by students, faculty, and other interested parties")
  Extract: numbering/captioning figures and equations, citing sources, lab-report structure.
- Michael Alley et al., "Writing as an Engineer or Scientist" (Penn State Leonhard Center) — section
  "Reports" / "Tutorial: Reports". https://writing.engr.psu.edu/
  Extract: report structure, model documents. (Free web resource; no license statement on page.)

## 00-12 Engineering Judgment & Sanity Checks

Sources: `@ocw-18-335j`

- GAP (no single authoritative free text). Partial: NASA Lessons Learned Information System,
  https://llis.nasa.gov/ (public NASA site; the landing page itself states no terms), and the MCO
  report above as a case study. 18.335J Week 2 "Lecture 5: Condition Numbers" for the numeric side.

## 00-13 Safety & Standards for Home Experiments

Sources: `@nar-safety-code`

- National Association of Rocketry, "Model Rocket Safety Code" — rules 1 "Materials" through 11
  "Recovery Safety". https://nar.org/ModelRocketSafetyCode
  Extract: materials, motors, launch-site distance table, recovery; use as the model for a lab safety code.
- GAP for general electrical/chemical/workshop safety: no single authoritative free standard found
  (relevant NFPA/OSHA standards are either paywalled or regulatory text, not instructional).

## 00-14 Measurement, Calibration & Uncertainty

Sources: `@jcgm-100-2008`, `@nist-tn-1297`

- JCGM 100:2008, "Evaluation of measurement data — Guide to the expression of uncertainty in
  measurement" (GUM). Publications list: https://www.bipm.org/en/committees/jc/jcgm/publications ;
  PDF: https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf ; DOI https://doi.org/10.59161/JCGM100-2008E
  Extract: Type A/Type B evaluation, combined standard uncertainty (law of propagation), coverage factor k.
  Note: JCGM 100:2008/Amd.1:2026 and JCGM GUM-1:2023, GUM-6:2020, GUM-5:2026 also listed (free).
- B. N. Taylor and C. E. Kuyatt, NIST Technical Note 1297 (1994 ed.), "Guidelines for Evaluating and
  Expressing the Uncertainty of NIST Measurement Results". https://www.nist.gov/pml/nist-technical-note-1297
  Extract: the short, practical US-government restatement of GUM; reporting-format examples.



## 00-15 First Experiment: Predict, Measure, Compare

Sources: `@jcgm-100-2008` (supplementary-free here; required in 00-14)

- Original instruction. The propagation and reporting rules reuse the GUM extract from 00-14. Procedure,
  bill of materials and safety box live in `projects/p00-first-experiment/`.

## Lessons without a mapped free source

- **00-02 Git & Collaboration**: original instruction; the authoritative reference is the manual installed
  with Git (`git help <command>`), as the lesson states.
- **00-08 Reading Data**: `@numpy-userguide` for array I/O and `@nasa-mco-mib` as the worked telemetry
  case; Python `csv` and `h5py` documentation linked directly from the lesson.
- **00-11 Reading a Technical Paper**: see the NTRS pick in [phase-04.md](phase-04.md#00-11-reading-a-technical-paper-ntrs-pick), `@ntrs-naca-tr-824-summary-of-airfoil-data`.
