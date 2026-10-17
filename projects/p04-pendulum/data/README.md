# Data folder: p04-pendulum

- `prediction.md`: lengths, predicted periods T₀ ± u, the frame term per period, the predicted amplitude excess at 10°, 20°, 30° with detectable / not-detectable verdicts, the source of g. Dated, written before any timing, never edited afterwards.
- `calibration.csv`: columns `date,item,reference,reading,unit,note`, for example the frame-rate check (`frame_rate,stopwatch 10.00 s,301,frames,30 fps setting`) and the bob's diameter and mass.
- `raw/YYYY-MM-DD_<tier>_lengths.csv`: columns `length_id,reading,L_to_bob_centre_m,note`, three readings per length.
- `raw/YYYY-MM-DD_<tier>_timings.csv`: columns `length_id,L_m,run,amplitude_start_deg,amplitude_end_deg,n_swings,start_frame,end_frame,fps_nominal,video_file,note`. The start and end frames are the frames in which the string crosses the 0° line in the agreed direction; n_swings full swings separate them. For the enhanced (break-beam) tier: `length_id,L_m,run,amplitude_start_deg,amplitude_end_deg,passage,time_s,note`, one row per beam interruption.
- `reduced/`: outputs of the analysis script (per-length T and g tables, the two plots).

Never edit a raw file. Fix mistakes in the analysis script and say so in the report. Keep the videos; they are raw data too.

If you have no kit and use `make_synthetic_pendulum_data`, save its output as `reduced/SYNTHETIC_seed<seed>.csv`, put **SYNTHETIC** in every file name, table caption and plot title derived from it, and do not place it in `raw/`.
