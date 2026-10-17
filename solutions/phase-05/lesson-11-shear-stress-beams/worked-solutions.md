# P05.11  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe script that calls it, given to the figures shown. Material values are typical, not design allowables.

## 1. A rectangle and a pin (s = 15: V = 3500 N)

**Assumptions.** Elementary beam theory; τ uniform across the width of each cut.

(a) A = 25 × 60 = 1500 mm². τ_avg = 3500/1500 = **2.333 MPa**; τ_max = 1.5 τ_avg = **3.500 MPa** at the neutral axis; at a quarter depth from the top, y = h/4 from the centroid, τ = 1.5(V/A)(1 − 4/16) = (9/8)(V/A) = **2.625 MPa**.

(b) A = π × 10² = 314.2 mm²: τ_avg = **11.14 MPa**, τ_max = 4/3 × 11.14 = **14.85 MPa**.

(c) ∫_{−h/2}^{h/2} (3V/2bh)(1 − 4y²/h²) b dy = (3V/2h)[h − (4/h²)(h³/12)] = (3V/2h)(2h/3) = **V** ✓.

**Independent check.** VQ/(It) at the neutral axis with Q = b(h/2)(h/4) = 25 × 30 × 15 = 11 250 mm³, I = 25 × 60³/12 = 450 000 mm⁴, t = 25 mm: 3500 × 11 250/(450 000 × 25) = 3.500 MPa ✓.

**Common mistakes.** Measuring y from the top instead of the centroid; using the average stress for a pin.

## 2. Derivations

(a) The part above the cut (area A′) has end forces ∫_{A′} σ dA at z and z + dz. With σ = −My/I the difference is −(dM/I)∫_{A′} y dA = −dM Q/I. The cut face (width t, length dz) carries τ t dz. Axial equilibrium: τ t dz = −dM Q/I; with dM/dz = −V, **τ = VQ/(It)** (positive in the sense that balances the part). The derivation is equilibrium applied twice (axial equilibrium of the part, and dM/dz = −V from the slice); the flexure stress it uses came from compatibility and the constitutive law in `P05.10`.

(b) I = t_w h_c³/12 + 2bt_f(h_c/2)². Half-flange, from the tip (y = ±h_c/2 constant): |q(s)| = (S_y/I) t_f (h_c/2) s, linear, reaching **q_f = S_y t_f b h_c/(4I)** at the junction. Web, from the top junction (incoming 2q_f), with y = h_c/2 − s: |q(s)| = 2q_f + (S_y/I) t_w (h_c s/2 − s²/2), a parabola with its maximum at s = h_c/2. Its integral over the web is 2q_f h_c + (S_y t_w/I)(h_c³/4 − h_c³/6) = (S_y/I)[2bt_f(h_c/2)² + t_w h_c³/12] = **S_y** exactly. The flanges are horizontal, so their flows have no vertical component; their horizontal forces cancel in pairs (the two half-flanges carry equal and opposite flows).

(c) The flange flow is linear from 0 to S_y t_f b h/(2I), so each flange carries F = S_y t_f b² h/(4I). The two flange forces are equal and opposite, h apart: a couple Fh. For the resultant S_y to have the same moment about the web, S_y e = Fh, e = t_f b² h²/(4I). With I = t_w h³/12 + bt_f h²/2: **e = 3b²t_f/(6bt_f + ht_w)**. As t_w → 0, **e → b/2**: all the bending stiffness is in the flanges and the flange forces are as large as they can be. As b → 0, **e → 0**: a flat plate's shear centre is on the plate.

**Common mistakes.** Starting the web integration at zero instead of at the flanges' incoming flow; using the overall depth h instead of the centreline depth h_c; forgetting the factor 2 for two half-flanges.

## 3. Angles and tees (s = 15: a = 55 mm)

(a) In thin-wall theory the flow in a straight leg acts along that leg. If all legs meet at one point, every flow's line of action passes through it, so the flows have zero moment about it, and their resultant passes through it: **the shear centre is at the leg junction** (the corner of the angle, the flange–stem junction of the T). `shear_centre` returns (0, 0) for both the 40 × 40 × 3 angle with the corner at the origin and the 60 × 4 / 50 × 3 T with the junction at the origin (to round-off).

(b) Centroid at **(a/4, a/4) = (13.75, 13.75) mm** from the corner. About centroidal axes parallel to the legs: I_xx = I_yy = **103 984 mm⁴**, I_xy = **−62 391 mm⁴** (= −ta³/8). A vertical 500 N load through the centroid is 13.75 mm from the shear centre (the corner): torque **6.875 N·m**.

(c) I_xy ≠ 0, so a vertical load bends the angle about an inclined neutral axis: it deflects sideways as well as down. Use the general (unsymmetric) formula of `P05.10`, or work in principal axes (the angle's axes of symmetry at ±45°, where I_xy = 0). The torque in (b) twists it as well, so a real angle stiffener should be loaded through its corner or restrained by the skin it is attached to.

**Common mistakes.** Putting the shear centre at the centroid; using the symmetric bending formula for an angle loaded parallel to a leg.

## 4. The web and the joint of your own spar (s = 15: h = 0.170 m, b_f = 46 mm, V = 13 601.8 N)

(a) I = **6.739 × 10⁶ mm⁴**. τ_NA = **30.43 MPa**, τ at the top of the web = **24.76 MPa**, V/(h_w t_w) = 13 601.8/(150 × 3) = **30.23 MPa**. The web carries **94.4 %** of V.

(b) The centreline model gives a peak web stress of **30.44 MPa**, **0.03 %** above (a).

(c) Q_cap = 46 × 10 × 80 = 36 800 mm³; q = VQ_cap/I = **74 272 N/m**; pitch for two rows at 4 kN with factor 1.5 = **71.8 mm**.

(d) MS = 331/(1.5 × 30.43) − 1 = **+6.25**: strength is not the issue for the web (Problem 6).

**Independent check.** The fraction carried by the caps, 1 − 0.944 = 5.6 %, is consistent with the lesson's 5.4 % for a slightly wider cap and deeper section.

**Common mistakes.** Using Q of half the section for the joint (that is the neutral-axis Q); using one row's pitch.

## 5. A control-surface channel spar (s = 15: h = 95 mm)

(a) Shear centre at **e = 9.818 mm** behind the web by `shear_centre` and by (6) (identical); centroid **5.806 mm** in front of the web.

(b) A 1000 N hinge reaction through the web is 9.818 mm from the shear centre: torque **9.818 N·m**.

(c) With 10 mm inward lips: **e = 13.839 mm**, centroid **8.571 mm**, I_xx from **330 917** to **388 983 mm⁴**. The lips move the shear centre **away** from the web: the flow in each lip enters the flange at its tip, so the flange flow no longer starts from zero, the flange forces and hence their couple grow faster than I does.

(d) Putting the hinge line on the shear-centre line (14 mm behind the web for the lipped section) would remove the torque, but hinge position is also set by aerodynamic balance (hinge moment), mass balance for flutter, and the space between spar and skin. In practice the spar is often closed into a D-cell or torque tube so that the torque is carried efficiently, rather than moving the hinge.

**Common mistakes.** Starting the lipped section's flow at the flange tip instead of the lip's free edge; expecting lips to pull the shear centre toward the web.

## 6. Can the web be thinner?

V = 13 974.5 N, web 160 mm deep, stiffeners every 200 mm: b/a = 0.8, k_s = 5.35 + 4 × 0.64 = **7.91**; E = 71.7 GPa, ν = 0.33.

| t (mm) | τ (MPa) | MS strength (F_su, 1.5) | τ_cr (MPa) | τ/τ_cr at limit |
|---|---|---|---|---|
| 3.0 | 29.1 | +6.58 | 184.0 | 0.16 |
| 2.0 | 43.7 | +4.05 | 81.8 | 0.53 |
| 1.0 | 87.3 | +1.53 | 20.4 | 4.27 |
| 0.5 | 174.7 | +0.26 | 5.1 | 34.2 |

(The strength margins use τ = V/(bt), the average web stress, and F_su = 331 MPa.)

(c) A tension-field web pulls diagonally on the caps (adding bending in the caps between stiffeners and an inward pull), pushes the stiffeners in compression, puts the web-to-cap fasteners in tension as well as shear, and leaves permanent wrinkles if it buckles at limit load. Other limits: fastener bearing and pull-through in thin sheet, minimum gauge for handling and corrosion, fatigue and crack growth from the wrinkle lines, damage tolerance after a strike, and the web's shear stiffness, which contributes to the wing's torsional stiffness (flutter).

(d) Setting τ_ult = τ_cr gives t³ = 1.5Vb/K with K = k_sπ²E/(12(1 − ν²)): **t = 1.86 mm** keeps the web unbuckled up to ultimate load, **1.62 mm** up to limit load. Recommendation: about **2 mm**, unless a tension-field analysis of web, stiffeners and caps (and a test) is done, which could justify thinner. 0.5 mm is far into the post-buckled regime (34 times τ_cr) and would need heavier stiffeners and caps that may cancel the saving.

**Common mistakes.** Treating strength margins as sufficient for thin webs; using b = a or the full spar depth in τ_cr.

## Check-yourself answers

1. Q is the first moment of the area on one side of the cut (between the cut and the free surface) about the neutral axis (the centroidal axis).
2. 1.5 V/A, at the neutral axis.
3. The same Q is divided by a different width: the stress just inside the flange is t_w/b_f times the stress in the web just below the junction, a jump by the factor b_f/t_w.
4. q = τt, the shear force per unit length of wall, in N/m. q = 0 at a free edge; at a junction the flows in sum to the flows out.
5. τ ≈ V/(h_w t_w) = 14 000/(160 × 3) = 29 MPa.
6. The point through which a shear force must act for the beam to bend without twisting. A channel loaded through its web is loaded off its shear centre (which is behind the web), so it twists as well as bends.
7. At the corner where the two legs meet: the flow in each leg acts along the leg, and both legs pass through the corner, so the flows have no moment about it.
8. p = P_allow/(1.5 q) = 4000/(1.5 × 50 000) = 53 mm.
