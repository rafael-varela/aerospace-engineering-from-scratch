# P10.21  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer. Nothing here requires XFOIL, AVL or OpenVSP; where a problem says "if you have the tool", that part is optional practice and its output must be stored as a labelled CSV (part 4 of the lesson). Never present a number from your own code, or from a test fixture, as tool output.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: what a polar file must say)

A colleague sends `naca0012.csv` with the columns alpha, CL, CD, CM and nothing else.

(a) Show that `read_tool_polar` refuses it, and write the six header lines it would need for an XFOIL polar at Re = 6 × 10⁶, M = 0.15, Ncrit = 9, 160 panels, an iteration limit of 100, free transition, run with XFOIL 6.99.
(b) For each of Re, Mach, Ncrit, the panel count and the iteration limit, give one sentence on how changing it could change the polar.
(c) Which of those settings are written into XFOIL's own polar file, and which must the person who ran it add by hand?

Answer form: a refusal message; six lines; five sentences; two lists.

## 2. (Derivation: comparison metrics that mean something)

(a) For a least-squares line cl = m α + c over n points, derive α_L0 = −c/m and, from the covariance of (m, c), an expression for its standard uncertainty. Why is α_L0's uncertainty large when the fitted α range does not straddle the zero-lift angle?
(b) Two polars are tabulated at different angles. Explain why the RMS cl difference must be computed at common angles, why interpolating the second polar onto the first's angles is acceptable inside the linear range, and why it is not acceptable near stall.
(c) The lesson's own code gives a lift slope 12.7 % above Ladson's at Re ≈ 6 × 10⁶. Show that a constant ratio between two slopes produces an RMS cl difference that grows linearly with the largest |α| in the range. What does that imply about quoting "RMS difference" without the range?

Answer form: two derivations and a formula; a short argument; a derivation and a sentence.

## 3. (Unfamiliar scenario: a suspicious drag value)

A student report claims, from "a quick XFOIL run", that a NACA 0012 at α = 0 and Re = 6 × 10⁶ has cd = 0.0045, and uses it to size an aircraft's drag.

(a) Compare with Ladson's free-transition and wraparound-grit values at that Reynolds number (Tables I and VII). Is 0.0045 plausible for the wind-tunnel model? For a wing with insect debris and rivets?
(b) List the settings that, if unrecorded, could explain a low value (think about Ncrit and transition).
(c) Write the header the report's polar should have carried, and the sentence you would add to the report.

Answer form: two published values and a judgment; a list; a header and a sentence.

## 4. (Numerical, uses your code: your own Reynolds number)

Pick the Ladson Table I block with index s mod 4 in the list [2.00, 3.94, 5.97, 8.86] × 10⁶.

(a) Run `validate_against_ladson` at that Reynolds number. Report the published and computed cl_α (per degree), α_L0, cd_min and the RMS cl difference over −4.5° to 8.5°.
(b) Tabulate cd(own)/cd(published) at each published angle in that range. Where is the own code worst, and why?
(c) Give the own-code polar its header (`own_polar(...)["meta"]`) and write it with `write_tool_polar`; read it back with `read_tool_polar`. Why should your own code's output be labelled exactly like a tool's?

Answer form: eight numbers; a table and a sentence; a file and a sentence.

## 5. (Numerical, uses your code: which setting matters?)

For the own code at Re = 5.97 × 10⁶, α = 0° and 4°:

(a) Change the panel count (80, 160, 320 panels). How much do cl and cd change?
(b) Change the Squire–Young station (x/c = 0.98, 0.99, 0.995). How much does cd change?
(c) Force transition at 0.05c on both surfaces and compare with Ladson's Table VII value at α = 0, Re = 6.00 × 10⁶. Why is the grit-strip drag still higher than your forced-transition estimate?
(d) Rank the three settings by their effect on cd, and say which must be recorded in a polar header.

Answer form: a table of cl and cd; three values; two numbers and a reason; a ranking.

## 6. (Judgment: tool, code and tunnel disagree)

Suppose you have run XFOIL (if you have it) and your own code, and both disagree with Ladson's data in cl_α, in different amounts.

(a) Which of the three is the reference, and why? What corrections has the reference already had applied (see the data README)?
(b) Your one-way code overpredicts cl_α by about 13 %. What physical effect is missing, and which XFOIL formulation (from its user primer) includes it?
(c) XFOIL's predicted transition depends on Ncrit. What does Ncrit represent, and why might the best Ncrit for a low-turbulence pressure tunnel differ from that for a flight test?
(d) Write the two-sentence statement you would put in a report when no tool was available, as in this lesson.

Answer form: four short paragraphs.
