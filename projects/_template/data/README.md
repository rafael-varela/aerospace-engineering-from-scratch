# Data folder

- `prediction.md`: your prediction with uncertainty, written before any measurement.
- `calibration.csv`: instrument, reference value, reading, date.
- `raw/`: untouched measurements, one CSV per session, named `YYYY-MM-DD_<tier>_<run>.csv`, with a header row and SI units in column names (for example `load_N,deflection_mm`). Never edit a raw file; fix mistakes in the analysis script and say so.
- `reduced/`: outputs of the analysis script.
