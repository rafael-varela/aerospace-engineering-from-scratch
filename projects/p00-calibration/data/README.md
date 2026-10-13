# Data folder: p00-calibration

- `prediction.md`: predicted uncertainties, readings, and error type, written before any measurement.
- `raw/YYYY-MM-DD_<tier>_scale.csv`: columns `sequence,stack_size,reference_mass_g,reading_g,position,direction` where `position` is `center|N|S|E|W` and `direction` is `up|down`.
- `raw/YYYY-MM-DD_<tier>_length.csv`: columns `sequence,object,instrument,reading_mm`.
- `calibration.csv`: the fitted offset, scale factor, their uncertainties, and the date, one row per instrument.
- `instruments.md`: the final uncertainty statement for each instrument, for reuse in later projects.
- `reduced/`: outputs of the analysis script.

Never edit a raw file. Fix mistakes in the analysis script and say so in the report.
