# Data folder: p00-first-experiment

- `prediction.md`: heights, predicted times with uncertainty, the uncertainty budget, and the drag estimate, written before any drop.
- `raw/YYYY-MM-DD_<tier>_heights.csv`: columns `height_id,reading_m,reference` (reference = `floor|pad_top`).
- `raw/YYYY-MM-DD_<tier>_drops.csv`: columns `height_id,drop,video_file,frame_rate_fps,release_frame,impact_frame,note`. For the stopwatch variant: `height_id,drop,time_s,note`.
- `reduced/`: outputs of the analysis script, including the two plots.

Never edit a raw file. Fix mistakes in the analysis script and say so in the report. Keep the videos; they are your raw data too.
