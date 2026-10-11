# References

`sources.yaml` is the machine-readable manifest of every external resource the curriculum cites: title, creator, version, canonical URL, license or usage terms, access tier (`required-free`, `supplementary-free`, `optional-paid`), mapped lessons, verification date, and caveats. `scripts/check_references.py` validates it in CI; the weekly `--links` run checks reachability (reachability is not content verification, so each source also gets a periodic human review).

`phase-XX.md` and `bridge-XX.md` map each lesson of a phase to exact chapters or sections of those sources, with a one-line note on what the learner should extract, and mark any `GAP` where no adequate free source exists and original instruction must be written.

Rules: every required reading is free to access; paid books appear only under "Optional textbooks for deeper study", never as prerequisites; we link to sources and write original explanations rather than reproducing protected text; free to read does not mean free to copy.
