# Phase 10 shared datasets

Authentic, public-domain wind-tunnel data that several Phase 10 lessons read (`P10.06`, `P10.13`, `P10.14`, `P10.15`, `P10.21`, `P10.22`, `P10.23`). Nothing in this folder is synthetic, and nothing here was produced by this curriculum's own codes or by any home rig. Read the CSVs with any CSV reader that skips lines starting with `#`. Never edit a value: apply corrections in code, and say so.

Both files were transcribed on 2026-10-17 from the NASA Technical Reports Server (NTRS) scans by reading the rendered page images, not the scans' OCR text, which garbles decimal points and columns. The sections below are the transcription notes.

## NACA 0012 section data at M = 0.15 (Ladson, NASA TM-4074)

File: `ladson_naca0012_M015.csv` (161 data rows + header)

### Source

Ladson, Charles L.: *Effects of Independent Variation of Mach and Reynolds Numbers on the
Low-Speed Aerodynamic Characteristics of the NACA 0012 Airfoil Section.* NASA TM-4074,
NASA Langley Research Center, October 1988.

- NTRS ID: 19880019495
- Citation page: https://ntrs.nasa.gov/citations/19880019495
- PDF: https://ntrs.nasa.gov/api/citations/19880019495/downloads/19880019495.pdf
- Rights: US Government work, public domain (NASA).

### Tables transcribed

| CSV `table` | Report table | Condition | PDF pages | Printed pages | Reynolds-number blocks (rows) |
|---|---|---|---|---|---|
| I | Table I, "Force Coefficients at M = 0.15 for Free Transition" | free transition | 8–9 | 6–7 | 2.00e6 (14), 3.94e6 (12), 5.97e6 (18), 8.86e6 (16), 11.90e6 (16), 18.90e6 (17) |
| VII | Table VII, "Force Coefficients at M = 0.15 for Transition Fixed With No. 60-W Grit" | fixed transition | 18–19 | 16–17 | 2.00e6 (14), 4.00e6 (14), 6.00e6 (14), 8.95e6 (15), 12.10e6 (11) |

Columns: `table, transition, mach, reynolds, alpha_deg, cd, cl, cm, l_over_d`. Each value is
copied exactly as printed (same number of digits), with a leading zero added (".00730" becomes
"0.00730"). A blank `cd` or `l_over_d` means the report prints no value there (post-stall
points with no wake-survey drag).

### Facility, model and measurements (from the report text, pp. 1–3)

- Langley Low-Turbulence Pressure Tunnel (LTPT), a two-dimensional test.
- Model: machined aluminium. Chord c = 23.66 in. (60.10 cm). Span = 36.00 in. (91.44 cm),
  which is the full tunnel width. The model was polished with No. 400 carborundum paper, and
  its ordinates were within 0.0002c of the design values.
- `cl` and `cm` (about 0.25c) come from integrating the model surface pressures. The model
  has 81 orifices: 24 in each chordwise row on the upper and lower surfaces, plus 33 in spanwise rows.
- `cd` comes from integrating a wake survey. The fixed rake sits at midspan, about 1c behind
  the trailing edge, and has 91 total-pressure tubes and 5 static-pressure tubes.
- `l_over_d` is the tabulated ratio cl/cd.
- Reynolds number is based on chord.
- All data include the standard low-speed 2-D wall corrections of Allen & Vincenti (NACA
  Rep. 782). The report says these amount to about 2% of the measured coefficients.
- Repeatability (from repeated alpha = 0 points): cd within 0.0002, normal-force coefficient
  within 0.004, cm within 0.0002.
- Fixed transition: carborundum strips about 0.01c wide at 0.05c on both surfaces. **Table VII
  ("No. 60-W") uses the NACA-type *wraparound* strip.** This strip runs from 0.05c on the upper
  surface, around the leading edge, to 0.05c on the lower surface, with No. 60 grit. The
  report includes it for comparison with older NACA data. Its drag (roughly 0.0085–0.011
  near zero lift) is therefore much higher than with free transition.

### Transcription method and checks

1. The OCR text layer (pdftotext) is badly garbled: decimal points are read as "-" or "*",
   "10^6" as "lo6", and columns are interleaved. In the Table VII R = 6.00e6 block, for
   example, the OCR swaps the cm values of the first two rows. For this reason the OCR was
   used only to locate the tables.
2. Every value was read from the rendered page images (pdftoppm at 200 dpi, PDF pp. 8, 9,
   18, 19). The scan is a clean typewritten table, and no digit was ambiguous on the images.
3. Automated checks on the transcription:
   - |cl/cd − l/d| / |l/d| < 2% for all 158 rows that have cl, cd and l/d. For |l/d| < 1
     the check uses an absolute tolerance of 0.05. **No failures.** The largest relative
     residual is 0.49% (Table VII, R = 12.1e6, alpha = 0.01: cl/cd = −1.045 vs. printed
     −1.04), which is rounding of a small number. All other residuals are below 0.2%.
   - cd > 0 in every row. The sign of cl matches the sign of alpha for |alpha| > 0.5°.
     Alpha increases monotonically within every block. No failures.
   - As a cross-check, the digit string of every transcribed cm value was found in the OCR
     text.
4. cm has no redundant check, unlike cl and cd, which are checked through l/d. Its signs
   were read visually and are unambiguous on the images.

### Points to flag (printed this way in the report, not transcription errors)

- Table I, R = 3.94e6: cd at alpha = 10.02 (0.00934) is lower than at 8.17 (0.00985).
  Table I, R = 5.97e6: cd at 6.20 (0.00680) is lower than at 4.18 (0.00760). Both are
  consistent with the printed l/d, so they are as printed. They probably reflect
  laminar-bucket and transition behaviour in free-transition wake surveys.
- Table I, R = 3.94e6, alpha = 6.00: cl is printed with only three decimals ("0.636").
- Table VII, R = 8.95e6, alpha = 15.08: cl = 0.7151 with cm = **+0.0767**, followed at 16.27
  by cl = 0.9974 with cm = −0.1008. This is a post-stall point; the values are checked
  against l/d = 3.28, but the cm sign is unusual.
- Table VII, R = 12.10e6, alpha = 12.37: cl = 0.9704, cd = 0.02305, cm = **+0.1128**,
  l/d = 42.10. The l/d agrees with cl/cd, but this is probably an unsteady stall-onset
  point and should be used with caution.
- Table I includes an R = 18.90e6 block (alpha starts at −2.07). It was transcribed along
  with the other blocks.
- Post-stall rows without cd: Table I R = 2.00e6 alpha 16.25; R = 5.97e6 alpha 17.65 and
  18.65. These rows keep cl and cm and leave cd and l/d blank.



## Finite-wing lift-curve slopes for checking a vortex-lattice code

File: `finite_wing_cla.csv` (15 rows: 12 distinct planforms from 3 NACA reports).
All sources are NACA Technical Notes: US Government works, public domain.

Columns: `report, ntrs_id, page, planform, aspect_ratio, taper_ratio (tip/root),
sweep_quarter_chord_deg, section, reynolds, mach, cla_per_deg, cla_source_note`.

**Note on taper ratio.** Both TN-1422 and TN-2445/TN-1093 give "ratio of root chord to tip
chord". The CSV converts this to tip/root (TN-1422: root/tip = 2.5, so 0.40).

### 1. NACA TN-1669: 9 untapered wings (AR 1.34 / 2.61 / 5.16, sweep 0 / 45 / 60 deg)

Goodman, Alex; and Brewer, Jack D.: *Investigation at Low Speeds of the Effect of Aspect Ratio
and Sweep on Static and Yawing Stability Derivatives of Untapered Wings.* NACA TN-1669,
August 1948.

- NTRS 19930082295: https://ntrs.nasa.gov/citations/19930082295
- PDF: https://ntrs.nasa.gov/api/citations/19930082295/downloads/19930082295.pdf

Test and model details:
- Facility: Langley stability tunnel, 6- by 6-ft test section. Dynamic pressure 24.9 psf,
  M = 0.13 (p. 7 / PDF p.6).
- Wings: untapered, NACA 0012 section normal to the leading edge. Planform figure 2 is on
  PDF p.15 (printed p.14) and gives A, S, b, mean chord and mounting point for each wing:
  - 0°: wings 1, 4, 7
  - 45°: wings 2, 5, 8
  - 60°: wings 3, 6, 9
  
  Because the wings are untapered, the leading-edge and quarter-chord sweeps are equal.
  The tips are streamwise with rounded corners (from the drawing).
- Reynolds number based on mean aerodynamic chord and V: Table I, PDF p.13 (printed p.12),
  which is printed upside-down on the scan. It ranges from 0.76e6 to 1.58e6.
- Corrections: jet-boundary corrections to alpha and drag (Silverstein & White, NACA Rep. 517),
  an approximate alpha correction for model deflection, and a strut tare for rolling moment
  only. There are **no** blockage or turbulence corrections (p. 6).
- CL_alpha source: **digitized from Figure 5, PDF p.21 (printed p.20)**, "Variation with sweep
  angle of the experimental and theoretical values of lift-curve slope". The report does not
  state the C_L range over which the slope was taken.
- Digitizing method: the page was rendered at 300 dpi and the gridlines were found from pixel
  rows and columns. The grid spacing is 0.01 per deg = 119.6 px, and the sweep axis has
  0° at x = 720 px and 60° at x = 2180 px. Each value is the centre dot of a symbol on that
  calibrated scale. Estimated reading uncertainty: ±0.0005 per deg, which is about ±5 px
  including symbol size and scan skew.

Caveats:
- The Reynolds numbers are low (about 1e6).
- The text (p. 7) says that for the low-aspect-ratio wings (1, 2 and 4) "the lift-curve slope
  increases" with lift coefficient. CL_alpha at AR 1.34 therefore depends on the fitting
  range. This nonlinear lift is also why AR 1.34 at 45° (0.0315) comes out higher than at 0°
  (0.0280). A linear VLM will not reproduce this.

### 2. NACA TN-1671: 45-deg swept wings, AR 2.61, taper 1.0 / 0.5 / 0.25

Letko, William; and Cowan, John W.: *Effect of Taper Ratio on Low-Speed Static and Yawing
Stability Derivatives of 45° Sweptback Wings with Aspect Ratio of 2.61.* NACA TN-1671,
July 1948.

- NTRS 19930082298: https://ntrs.nasa.gov/citations/19930082298
- PDF: https://ntrs.nasa.gov/api/citations/19930082298/downloads/19930082298.pdf

Test and model details:
- Same tunnel and conditions as TN-1669: q = 24.9 psf, M = 0.13.
- Sections are NACA 0012 normal to the quarter-chord line. Reynolds number based on mean
  aerodynamic chord is about 1.10e6 for taper 1.00 and 1.23e6 for taper 0.25 (text,
  PDF p.6). The value for taper 0.5 is not stated.
- Corrections: approximate jet-boundary corrections based on unswept-wing theory, applied to
  alpha and the longitudinal force. There are no blockage or turbulence corrections (PDF p.7).
- CL_alpha source: **digitized from Figure 6 (upper panel), PDF p.17 (printed p.15)**, for
  C_L < 0.6. Grid: 0.01 = 138 px at 300 dpi. Uncertainty ±0.0003 per deg.
- The text (PDF p.7) says taper "had very little effect on the lift-curve slope" below
  C_L = 0.6, which agrees with the digitized values.

Caveat: the taper-1.0 wing reads 0.0449 here but 0.0463 in TN-1669 (wing 5), about 3% apart.
This may be the same model tested in a different run, or a different fitting range. Both
values are kept as published. Treat roughly ±3% as the realistic experimental scatter.

### 3. NACA TN-1422: AR 9, taper 0.4 unswept wings (3 wings, 1 planform)

Sivells, James C.: *Experimental and Calculated Characteristics of Three Wings of NACA 64-210
and 65-210 Airfoil Sections With and Without 2° Washout.* NACA TN-1422, August 1947.

- NTRS 19930082118: https://ntrs.nasa.gov/citations/19930082118
- PDF: https://ntrs.nasa.gov/api/citations/19930082118/downloads/19930082118.pdf

Test and model details:
- Facility: Langley 19-ft pressure tunnel, about 35 psi absolute.
- Flow conditions: R ≈ 4.4e6 based on the mean aerodynamic chord of 1.769 ft, M ≈ 0.17.
- Planform: AR 9, root/tip chord 2.5. Quarter-chord sweep 0°, dihedral 3°. Washout is 2° or
  0°, applied as linear twist about the 0.25c line (PDF p.5).
- Model: solid steel, held to contour within 0.003 in.
- Corrections: tare and interference, plus streamline-curvature and jet-boundary corrections
  to alpha and drag (PDF p.5).
- CL_alpha source: **stated in Table I, "Summary of results...", PDF p.10 (printed p.9)**.
  Experimental dCL/dalpha per deg:
  - 64-210 with 2° washout: 0.086 (calculated 0.083)
  - 65-210 with 2° washout: 0.085 (calculated 0.085)
  - 65-210 with 0° washout: 0.085 (calculated 0.085)
  
  The same table gives alpha(L=0), CLmax, (L/D)max, CDo,min, Cm(L=0) and dCm/dCL. These
  values were read from the rendered page image and are not in the CSV.
- Caveat: the sections are cambered, but camber and washout do not affect a linear CL_alpha.
  The values have only two significant figures.

### Sources examined and not used

- NACA TR-431 (Zimmerman 1932, Clark Y small aspect ratios; NTRS 19930091504). Its
  coefficient tables (Tables II–VIII) are for 20° yaw only. Table IX summarizes CLmax,
  CDmin and L/D but gives no CL_alpha, and the zero-yaw lift curves appear only in figures.
- NACA TN-1093 (Shortal & Maggin 1946; NTRS 19930087879). This is a compilation with data
  in figures only.
- NACA TN-2445 (Purser & Spearman 1951; NTRS 19930083063). Figure 42 compares measured and
  computed CL_alpha, but its symbols do not identify individual wings, so it was not usable.
- NACA TN-1862 (Polhamus 1949; NTRS 19930082534). Table I lists 29 swept planforms (with
  section a0) used in a CL_alpha correlation. The measured values appear only in figure 3,
  without per-wing labels, so it was not digitized. It is a good lead for references
  (TN-1669, TN-1671, TN-1703, TN-1046).
- NACA TR-921 (DeYoung & Harper). Its experimental comparisons (figs. 10–12) cite TR-431,
  TM-1176 and TN-1093 and were not digitized.


### Reproducing the digitization

Render the cited page at 300 dpi (`pdftoppm -r 300 -f <page> -l <page> -png <report>.pdf out`), locate the grid lines, and read each symbol centre by linear interpolation between the grid lines. Expect to agree with the tabulated values within the stated reading uncertainty.
