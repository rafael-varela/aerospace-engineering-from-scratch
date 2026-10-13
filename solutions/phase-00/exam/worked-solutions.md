# Phase 00 Exam: worked solutions and pass thresholds

Attempt the exam first, in the time you set yourself. Then grade each question against this key and record the result per competency (1 reasoning, 2 analysis, 3 computation, 4 interpretation). **Pass threshold:** proficient (3 of 4) or better on every competency, where "proficient" means the reasoning is right and the result is right to within rounding, even if the presentation is rough.

## Part A

**1.** A model is a set of equations plus assumptions (here: a fixed sea-level temperature and pressure, a fixed lapse rate, hydrostatic equilibrium, dry air as an ideal gas). The real atmosphere on any given day differs from every one of those, so an altimeter can at best match the *standard* atmosphere's equations, not the air outside. "Matches exactly" also hides numerical and measurement uncertainty in the altimeter's own pressure sensor. The honest sentence is: "The model reproduces the ISA pressure–altitude relation to within X Pa over 0–11 km, and the sensor adds ±Y Pa, so altitude is uncertain by ±Z m, before any weather correction." *Proficient:* identifies that a model embodies assumptions and that a claim without an uncertainty is not a claim.

**2.** 1 lbf·s = 4.448 N·s, so values in pound-force seconds read as newton-seconds were too small by a factor of about 4.45. The ground software therefore modeled *too little* thruster effect per firing; the real spacecraft received 4.45 times the modeled impulse, and over months of small firings its reconstructed trajectory drifted from the real one until the orbit-insertion altitude was far too low. The practice: carry units with the number (a unit-aware type, or at minimum the unit in the interface specification and a check at the boundary) at the interface between the two programs, which is exactly where the interface control document had already specified SI. *Proficient:* factor ≈ 4.4–4.5, direction correct, names the software interface as the place to check.

**3.** Physical experiment: timing a dropped object with phone video; it builds calibration, noise handling, and the humility of a result that disagrees with you. Computational experiment: sweeping the time step of an orbit propagator and watching energy drift; it builds judgment about numerical error and model sensitivity. Data-based investigation: reducing a published wind-tunnel dataset to a drag polar; it builds skill in reading others' measurements, their conditions, and their uncertainty statements. None substitutes for another: simulation never teaches you that your ruler was mis-zeroed. *Proficient:* three correct examples and a sensible distinction.

## Part B

**4.** Readings: seven at 52 g, three at 53 g. Mean = 52.3 g. Sample standard deviation s = √(Σ(xᵢ − x̄)²/(n−1)) = √((7·0.09 + 3·0.49)/9) = √(2.1/9) = 0.483 g. Type A standard uncertainty of the mean u_A = s/√n = 0.483/√10 = 0.153 g. Type B from a 1 g display, uniform over one step: u_B = 1/√12 = 0.289 g. Type B dominates (the display rounds away more than the readings scatter). Combined: u_c = √(0.153² + 0.289²) = 0.327 g. Result: m = 52.3 ± 0.3 g (k = 1); expanded, U = 2 × 0.327 = 0.65 g, so m = 52.3 ± 0.7 g (k = 2). Significant figures: the uncertainty to one (at most two) figures, the value to the same decimal place. *Proficient:* both uncertainty types computed, combination in quadrature, correct dominance, sensible rounding.

**5.** t = (2h/g)^{1/2}. ∂t/∂h = t/(2h); ∂t/∂g = −t/(2g). Law of propagation (independent inputs): u_t² = (t/2h)² u_h² + (t/2g)² u_g², or in relative form (u_t/t)² = ¼[(u_h/h)² + (u_g/g)²]. Numbers: t = √(3.000/9.81) = 0.5530 s. Height term: 0.5530/(2 × 1.5) × 0.005 = 0.00092 s. Gravity term: 0.5530/(2 × 9.81) × 0.03 = 0.00085 s. u_t = √(0.00092² + 0.00085²) = 0.00125 s, so t = 0.5530 ± 0.0013 s. The height contributes slightly more, by a ratio of about 1.09 (relative inputs 0.33% vs 0.31%). *Proficient:* correct partial derivatives, quadrature sum, ratio within ±10%.

**6.** For a linear lapse T = T₀ − L z, hydrostatics plus the ideal gas law integrate to p/p₀ = (T/T₀)^{g/(L R)} with g/(L R) = 9.80665/(0.0065 × 287.05) = 5.256. At 5 km, T = 288.15 − 32.5 = 255.65 K, T/T₀ = 0.8872, so p/p₀ = 0.8872^{5.256} = 0.533 and p ≈ 54.0 kPa. Density: ρ/ρ₀ = (p/p₀)/(T/T₀) = 0.533/0.887 = 0.60, so *more* than half, because the air at 5 km is colder and colder air is denser at a given pressure; density falls more slowly with altitude than pressure does. *Proficient:* pressure within ±1 kPa, correct direction and reason for density.

## Part C

**7.** A proficient answer looks like:

```python
def drag_force(rho: float, v: float, cd: float, area: float) -> float:
    """Aerodynamic drag 0.5*rho*v**2*cd*area, SI units."""
    if area < 0 or rho < 0 or cd < 0:
        raise ValueError("rho, cd and area must be non-negative")
    return 0.5 * rho * v**2 * cd * area

def test_hand_value():
    # sea level, 50 m/s, Cd 0.3, 1 m^2: 0.5*1.225*2500*0.3 = 459.375 N
    assert drag_force(1.225, 50.0, 0.3, 1.0) == pytest.approx(459.375, rel=1e-9)

def test_zero_speed_gives_zero_drag():
    assert drag_force(1.225, 0.0, 0.3, 1.0) == 0.0

def test_negative_area_rejected():
    with pytest.raises(ValueError):
        drag_force(1.225, 50.0, 0.3, -1.0)
```

Why not `==` on the float result: the product is evaluated in binary floating point, so a hand-computed decimal and the computed value can differ in the last bit depending on operation order; `pytest.approx` with a stated tolerance tests what you actually mean (agreement to within a precision), and makes the tolerance an explicit engineering choice. *Proficient:* the three tests present and correct, the explanation names rounding and a tolerance.

## Part D

**8.** Pattern: good agreement at 0.5 and 1.0 m, measured times *longer* than predicted at 1.5 and 2.0 m. Candidate causes, ranked:

1. **Release artifact at the higher drops.** Reaching up to 1.5–2.0 m, the learner releases less cleanly; the first visible motion in the video may be a pre-release wobble, so the release frame is read too early and the fall time comes out too long. It affects only the awkward heights, which matches.
2. **Height reference error at the higher drops.** If the tape measure was read at an angle or the wall mark was above the true release point, the measured h is too large for those two drops; that would make the *prediction* too long, not the measurement, so by itself it gives the opposite sign unless the mark was below the release point. Keep it on the list; it is testable.
3. **Air drag.** It lengthens the fall and grows with speed, so it does act in the right direction and more at larger h; but for a dense 100 g object the effect at 2 m is well under 1%, far smaller than the 2.8–3.1 normalized errors imply, unless the object is light.

Best discriminating follow-up between 1 and 3: repeat the 2.0 m drops with a clip release (hands out of frame) and the same object. If the normalized error falls below 1, the release artifact was the cause. If it persists, swap to a denser object of the same size at 2.0 m: if the error now falls, drag was the cause; if it still persists, re-measure the height from floor to the object's underside at release (cause 2). *Proficient:* at least two plausible causes with the right sign, one follow-up that actually separates them, and the predicted outcome for each.
