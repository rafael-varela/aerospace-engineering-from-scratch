# P04.12  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact, which are numerical and which are SYNTHETIC (Monte Carlo on simulated images). Angles in degrees or arcseconds as asked.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a fibre and a window)

An optical fibre carries a laser-communication signal to a detector. Its core has n₁ = 1.48 and its cladding n₂ = 1.44 + s/1000.

(a) Compute the critical angle at the core–cladding wall.
(b) Show that the largest angle from the axis at which a ray from air can enter the fibre's end face and still be trapped satisfies sin θ_a = √(n₁² − n₂²) (the numerical aperture, NA). Compute the NA and θ_a. Check your result by tracing such a ray with `refraction_angle`.
(c) The detector sits behind an uncoated sapphire window (n = 1.76). What fraction is reflected at each surface at normal incidence, and what fraction of the light passes both surfaces (ignore multiple reflections)?

Answer form: one angle; one number and one angle; two numbers.

## 2. (Derivation: four results of this lesson)

(a) Derive the field of view (3) and the pixel scale p/f from the geometry of a lens focused at infinity. Show that the arctangent form is smaller than np/f and that the two agree for small fields.
(b) A spot's N photons land at independent positions with standard deviation σ (pixels), and each is recorded only as the centre of its pixel. Derive the variance of the centroid, (σ² + 1/12)/N. State the assumption under which the pixel's contribution is uniform with variance 1/12.
(c) Derive the quarter-bridge output (6) for one gauge R(1 + x) and three resistors R, and its linear approximation. Show that a full bridge with two gauges at R(1 + x) and two at R(1 − x), placed correctly, gives exactly V_ex x.
(d) Show that the subsonic isentropic impact pressure q_c = p[(1 + (γ − 1)M²/2)^{γ/(γ−1)} − 1] reduces to ½ρV² at low Mach number, and find the first correction term in M².

Answer form: four derivations.

## 3. (Unfamiliar scenario: a CubeSat battery thermistor)

A CubeSat's battery must be kept between −20 °C and +60 °C, and its heater controller reads the temperature through the lesson's NTC thermistor (10 kΩ at 25 °C, B = 3950 K) at the bottom of a divider on 3.3 V, read by a 12-bit ADC with a 3.3 V reference. The tightest control is needed near +20 °C.

(a) Find R(20 °C). Choose the top resistor to maximise |dV/dT| at 20 °C; show analytically that the optimum is R_top = R(20 °C).
(b) Tabulate dV/dT in mV/K and the temperature resolution in K per count at 20 °C for R_top = 2, 5, 20 and 50 kΩ and for the optimum.
(c) With the optimum resistor, compute the resolution at −20, 0, 40 and 60 °C. Where in the range is the reading coarsest?
(d) A designer argues that self-heating is negligible. With the optimum resistor at 20 °C, compute the power dissipated in the thermistor, and say what dissipation constant (mW/K) would make the self-heating error equal to one count.

Answer form: one resistance, a five-row table, a four-row table, one power and one dissipation constant.

## 4. (Code: your own star tracker)

Design a star tracker with focal length f = 16 + s/2 mm at f/1.4, a 2048 × 2048 detector of 3.45 μm pixels, stars to magnitude 5.5 and a 50 ms exposure. Use `StarTracker` and the functions of `main.py`.

(a) Compute the aperture, the pixel scale in arcseconds, the field of view in degrees and steradians, and the expected number of stars in the field.
(b) Compute the Airy radius at 600 nm in μm and in pixels. Does this tracker need defocus?
(c) Compute the photo-electrons from a magnitude-5.5 star. With σ = 1 pixel, 3 e⁻ read noise and seed s, compute the SYNTHETIC rms centroid error (2000 images) for 3 × 3, 5 × 5 and 7 × 7 windows, and the per-star and cross-boresight attitude errors in arcseconds.
(d) Compare with the lesson's baseline (0.75″). Which change made the difference: pixel scale, photons or star count?

Answer form: four numbers and a star count; two numbers and yes/no; one number and a three-row table; a paragraph.

## 5. (Code and design: a wing-spar strain channel)

A wing spar is instrumented to measure bending strain up to ±1500 με with gauges of gauge factor 2, excited at 3.3 V and read by a 12-bit ADC with a 0–3.3 V range, biased so that zero strain reads mid-scale.

(a) Compute the bridge output at 1500 με for a quarter bridge (exact and linear, and the linearity error in per cent) and for a full bridge.
(b) Choose the amplifier gain for each so that 1500 με gives a 1.6 V swing from mid-scale.
(c) Compute the strain resolution in με per count for the full-bridge channel, and confirm with `bits_needed` how many bits are needed to resolve 1 με over the ±1500 με range.
(d) Which error sources in the real channel are larger than one count? Name three and say how each is reduced.

Answer form: three voltages and one percentage; two gains; one resolution and a bit count; a list.

## 6. (Judgment: review a star-tracker proposal)

A student team proposes, for a 3U CubeSat that must point a camera to 10″: "We will focus the star tracker as sharply as possible to maximise the peak signal, centroid with a generous 9 × 9 window to capture all the light, and use the brightest single star in the field to save processing time." Their detector has 10 e⁻ read noise.

Write a short review (one page at most). For each of the three decisions, say whether it is right, give the quantitative consequence using results of this lesson, and propose the change. Then say which assumption of the lesson's budget you would check first with real night-sky images, and what result would make you reject the design.

Answer form: three verdicts with numbers; one test with a rejection criterion.
