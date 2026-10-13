# P00.11  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. Split the string into digit 1, digit 2, digits 3–4. 2. $m = d_1/100$, $p = d_2/10$, $t = d_{34}/100$. 3. Multiply $m$ and $t$ by 1200 mm. 4. Note that for 0009 the camber position is undefined.

## 2.
1. Write $y_c = ax^2 + bx + c$ and its derivative. 2. Impose $y_c'(p) = 0$ to get $b = -2ap$. 3. Impose $y_c(1) = 0$ to get $c$ in terms of $a$ and $p$. 4. Impose $y_c(p) = m$ to solve for $a$. 5. Substitute back and factor to match equation (2). 6. Differentiate both arcs; evaluate at $p$ (both zero), at 0 (front arc), and at 1 (rear arc) with $m = 0.02$, $p = 0.4$.

## 3.
1. Plot $C_d$ (horizontal) against $C_l$ (vertical), or the reverse as the report does; be consistent. 2. Lowest $C_d$ in the table and its $C_l$. 3. Compute $C_l/C_d$ for each row; the largest is the maximum L/D; mark the tangent from the origin on the sketch. 4. Run `drag_polar_summary(cl, cd)` and compare. 5. Explain: minimum drag is where the polar is lowest; maximum L/D is where the ratio is largest, which is further up the polar because lift grows faster than drag at first.

## 4.
1. Open the 2412 plot; identify the Re ≈ 6 million smooth curve by its symbol. 2. Choose about eight α values across −6° to +14° and read $C_l$; note the minor-division size as $u_{C_l}$. 3. Save the CSV. 4. Fit once; shift all $C_l$ by $+u$ and refit; by $-u$ and refit; drop the end points and refit. 5. Uncertainty on slope and $\alpha_{L0}$ = half the range over the four fits (or the max deviation from the first fit). 6. Convert $2\pi$ rad⁻¹ to deg⁻¹ and compare.

## 5.
1. Find the peak of each Reynolds-number curve; read $C_{l,\max}$ and $u$. 2. Compare the difference with $\sqrt{u_1^2 + u_2^2}$. 3. Fit the linear region of each with `lift_curve_fit`. 4. Compare slopes the same way. 5. Argue from where viscosity acts.

## 6.
1. $Re = \rho V c / \mu$ with the given values. 2. Compare with 3–9 million (ratio). 3. Decide using §2.6: stall and $C_{l,\max}$ depend strongly on Re; slope and $\alpha_{L0}$ weakly; $C_{d,\min}$ rises substantially at low Re. 4. Rank.

## 7.
1. Sentence 1: the question the report answers. 2. Sentence 2: tunnel, Reynolds numbers, surface conditions. 3. Sentence 3: a quantitative headline (slopes near 0.1/deg; minimum drag range). 4. Sentence 4: agreement with thin-airfoil expectation. 5. Sentence 5: two-dimensional, low-speed limitation. 6. Abstract: what the document contains, for whom. 7. One sentence on the difference.
