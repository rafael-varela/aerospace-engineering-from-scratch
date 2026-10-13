# P00.10  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Round the uncertainty first, to one significant figure (two if the leading digit is 1). Then round the value to the same decimal place. The concise form puts the uncertainty's digits in parentheses at the position of the last shown digit.
2. Write $\rho = m^{1} r^{-2} h^{-1}$. The exponents are $1$, $-2$, $-1$. Sign does not matter after squaring, but the 2 does.
3. Convert 20° to radians before squaring. Compare the correction you get with 2 %, and ask whether a cause that cannot produce the size of the effect can be "the" cause.
4. The linter checks the Results section only for uncertainties, and only recognises units in its list. Either fact gives you a fourth mistake it will miss.
5. The caption line already contains "Figure n", so exclude caption lines from the search. A regular expression such as `\bFig(?:ure|\.)\s*(\d+)` over the remaining lines gives you the referenced numbers.
6. The linter never compares the Prediction with the Results, and never checks whether the Conclusions follow. Put a prediction of 10 ± 1, a result of 20 ± 1, and conclude "agreement".
7. "Apparatus" for a computation is the software environment: versions, solver settings, random seeds, and the machine. "Calibration" becomes verification against a known case. "Safety" is honestly "none".
