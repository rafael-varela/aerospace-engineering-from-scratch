"""P07.03 Properties of Pure Substances: reference solution.

A property-table lookup for two working fluids of aerospace interest: water (steam; spacecraft
evaporative coolers) and anhydrous ammonia (the fuel of the X-15's XLR99 engine, burned with liquid
oxygen; also the coolant of large spacecraft thermal-control loops):

* saturation properties at a temperature by interpolation, the saturation pressure by linear or
  logarithmic (Clausius-Clapeyron-shaped) interpolation, and its inverse, the saturation temperature;
* phase identification from (T, p) and from (T, v), quality x and two-phase mixture properties;
* superheated-vapour lookup by interpolation in temperature and then in pressure;
* the complete state from (T, v) in any region, and the rigid-tank heating problem (closed system, first law);
* the independent checks of a table: Clausius-Clapeyron dp/dT = h_fg / (T v_fg), the identity
  h = u + p v to table round-off, the compressibility factor Z = p v / (R T) -> 1 at low pressure, and
  superheated-steam enthalpy differences against the NASA ideal-gas fit for H2O.

The tables are the NIST Chemistry WebBook values (SRD 69, Lemmon, McLinden & Friend) as published in
Yan, Introduction to Engineering Thermodynamics (LibreTexts), Appendices A and B, embedded below as given
constants (a subset of each table).

Units: tables and the public functions use T in deg C, p in kPa, v in m^3/kg, u and h in kJ/kg,
s in kJ/(kg K). Kelvin (T + 273.15) is used inside every thermodynamic formula.

Run:  python solution.py          (prints the investigation tables, writes property_tables.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-07/lesson-03-properties-pure-substances
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.optimize import brentq

R_U = 8.31446261815324            # J/(mol K), CODATA 2018
T_OFFSET = 273.15                 # K = deg C + 273.15
M_WATER = 18.01528e-3             # kg/mol (NASA CEA thermo.inp)
M_AMMONIA = 17.03052e-3           # kg/mol, from the NASA CEA atomic masses N 14.0067 and H 1.00794

# Water, saturated (temperature table): Yan Appendix A Table A1 (NIST SRD 69 values)
# columns: T_C, p_MPa, vf, vg (m^3/kg), uf, ug, hf, hg (kJ/kg), sf, sg (kJ/(kg K))
WATER_SAT = np.array([
    [0.01, 0.00061, 0.001000, 205.991, 0.00, 2374.90, 0.00, 2500.92, 0.0000, 9.1555],
    [5, 0.00087, 0.001000, 147.011, 21.02, 2381.78, 21.02, 2510.06, 0.0763, 9.0248],
    [10, 0.00123, 0.001000, 106.303, 42.02, 2388.65, 42.02, 2519.21, 0.1511, 8.8998],
    [15, 0.00171, 0.001001, 77.8755, 62.98, 2395.49, 62.98, 2528.33, 0.2245, 8.7803],
    [20, 0.00234, 0.001002, 57.7567, 83.91, 2402.32, 83.91, 2537.43, 0.2965, 8.6660],
    [25, 0.00317, 0.001003, 43.3373, 104.83, 2409.13, 104.83, 2546.51, 0.3672, 8.5567],
    [30, 0.00425, 0.001004, 32.8783, 125.73, 2415.91, 125.73, 2555.55, 0.4368, 8.4520],
    [35, 0.00563, 0.001006, 25.2053, 146.63, 2422.67, 146.63, 2564.55, 0.5051, 8.3517],
    [40, 0.00738, 0.001008, 19.5151, 167.53, 2429.39, 167.53, 2573.51, 0.5724, 8.2556],
    [45, 0.00960, 0.001010, 15.2521, 188.43, 2436.08, 188.44, 2582.43, 0.6386, 8.1633],
    [50, 0.01235, 0.001012, 12.0269, 209.33, 2442.73, 209.34, 2591.29, 0.7038, 8.0748],
    [55, 0.01576, 0.001015, 9.56428, 230.24, 2449.34, 230.26, 2600.09, 0.7680, 7.9898],
    [60, 0.01995, 0.001017, 7.66723, 251.16, 2455.90, 251.18, 2608.83, 0.8313, 7.9082],
    [65, 0.02504, 0.001020, 6.19354, 272.09, 2462.41, 272.12, 2617.50, 0.8937, 7.8296],
    [70, 0.03120, 0.001023, 5.03954, 293.03, 2468.86, 293.07, 2626.10, 0.9551, 7.7540],
    [75, 0.03860, 0.001026, 4.12894, 313.99, 2475.24, 314.03, 2634.60, 1.0158, 7.6812],
    [80, 0.04741, 0.001029, 3.40516, 334.96, 2481.56, 335.01, 2643.02, 1.0756, 7.6111],
    [85, 0.05787, 0.001032, 2.82584, 355.96, 2487.81, 356.02, 2651.33, 1.1346, 7.5434],
    [90, 0.07018, 0.001036, 2.35906, 376.97, 2493.97, 377.04, 2659.53, 1.1929, 7.4781],
    [95, 0.08461, 0.001040, 1.98056, 398.00, 2500.04, 398.09, 2667.61, 1.2504, 7.4151],
    [100, 0.10142, 0.001043, 1.67177, 419.06, 2506.02, 419.17, 2675.57, 1.3072, 7.3541],
    [105, 0.12090, 0.001047, 1.41838, 440.15, 2511.90, 440.27, 2683.39, 1.3634, 7.2952],
    [110, 0.14338, 0.001052, 1.20929, 461.26, 2517.67, 461.42, 2691.06, 1.4188, 7.2381],
    [115, 0.16918, 0.001056, 1.03584, 482.42, 2523.33, 482.59, 2698.58, 1.4737, 7.1828],
    [120, 0.19867, 0.001060, 0.891212, 503.60, 2528.86, 503.81, 2705.93, 1.5279, 7.1291],
    [125, 0.23224, 0.001065, 0.770026, 524.83, 2534.27, 525.07, 2713.10, 1.5816, 7.0770],
    [130, 0.27028, 0.001070, 0.668004, 546.10, 2539.53, 546.38, 2720.08, 1.6347, 7.0264],
    [135, 0.31323, 0.001075, 0.581729, 567.41, 2544.65, 567.75, 2726.87, 1.6872, 6.9772],
    [140, 0.36154, 0.001080, 0.508454, 588.77, 2549.62, 589.16, 2733.44, 1.7392, 6.9293],
    [145, 0.41569, 0.001085, 0.445962, 610.19, 2554.42, 610.64, 2739.80, 1.7908, 6.8826],
    [150, 0.47617, 0.001091, 0.392453, 631.66, 2559.05, 632.18, 2745.93, 1.8418, 6.8371],
    [155, 0.54350, 0.001096, 0.346460, 653.19, 2563.51, 653.79, 2751.81, 1.8924, 6.7926],
    [160, 0.61824, 0.001102, 0.306782, 674.79, 2567.78, 675.47, 2757.44, 1.9426, 6.7491],
    [165, 0.70093, 0.001108, 0.272431, 696.46, 2571.85, 697.24, 2762.81, 1.9923, 6.7066],
    [170, 0.79219, 0.001114, 0.242589, 718.20, 2575.73, 719.08, 2767.90, 2.0417, 6.6650],
    [175, 0.89260, 0.001121, 0.216581, 740.02, 2579.39, 741.02, 2772.71, 2.0906, 6.6241],
    [180, 1.00281, 0.001127, 0.193842, 761.92, 2582.83, 763.05, 2777.21, 2.1393, 6.5841],
    [185, 1.12346, 0.001134, 0.173901, 783.91, 2586.04, 785.19, 2781.41, 2.1875, 6.5447],
    [190, 1.25524, 0.001141, 0.156362, 806.00, 2589.01, 807.43, 2785.28, 2.2355, 6.5059],
    [195, 1.39882, 0.001149, 0.140892, 828.18, 2591.74, 829.79, 2788.82, 2.2832, 6.4678],
    [200, 1.55493, 0.001157, 0.127210, 850.47, 2594.20, 852.27, 2792.01, 2.3305, 6.4302],
    [205, 1.72430, 0.001164, 0.115078, 872.87, 2596.40, 874.88, 2794.83, 2.3777, 6.3930],
    [210, 1.90767, 0.001173, 0.104292, 895.39, 2598.31, 897.63, 2797.27, 2.4245, 6.3563],
    [215, 2.10584, 0.001181, 0.094680, 918.04, 2599.94, 920.53, 2799.32, 2.4712, 6.3200],
    [220, 2.31959, 0.001190, 0.086092, 940.82, 2601.25, 943.58, 2800.95, 2.5177, 6.2840],
    [225, 2.54972, 0.001199, 0.078404, 963.74, 2602.24, 966.80, 2802.15, 2.5640, 6.2483],
    [230, 2.79709, 0.001209, 0.071504, 986.81, 2602.90, 990.19, 2802.90, 2.6101, 6.2128],
    [235, 3.06253, 0.001219, 0.065298, 1010.04, 2603.20, 1013.77, 2803.17, 2.6561, 6.1775],
    [240, 3.34693, 0.001229, 0.059705, 1033.44, 2603.13, 1037.55, 2802.96, 2.7020, 6.1423],
    [245, 3.65117, 0.001240, 0.054654, 1057.02, 2602.67, 1061.55, 2802.22, 2.7478, 6.1072],
    [250, 3.97617, 0.001252, 0.050083, 1080.79, 2601.79, 1085.77, 2800.93, 2.7935, 6.0721],
    [255, 4.32289, 0.001264, 0.045939, 1104.77, 2600.49, 1110.23, 2799.07, 2.8392, 6.0369],
    [260, 4.69226, 0.001276, 0.042173, 1128.97, 2598.72, 1134.96, 2796.60, 2.8850, 6.0016],
    [265, 5.08529, 0.001289, 0.038746, 1153.41, 2596.45, 1159.96, 2793.49, 2.9307, 5.9661],
    [270, 5.50299, 0.001303, 0.035621, 1178.10, 2593.67, 1185.27, 2789.69, 2.9765, 5.9304],
    [275, 5.94639, 0.001318, 0.032766, 1203.07, 2590.33, 1210.90, 2785.17, 3.0224, 5.8944],
    [280, 6.41658, 0.001333, 0.030153, 1228.33, 2586.39, 1236.88, 2779.87, 3.0685, 5.8579],
    [285, 6.91466, 0.001349, 0.027756, 1253.92, 2581.81, 1263.25, 2773.73, 3.1147, 5.8209],
    [290, 7.44178, 0.001366, 0.025555, 1279.86, 2576.53, 1290.03, 2766.70, 3.1612, 5.7834],
    [295, 7.99911, 0.001385, 0.023529, 1306.19, 2570.49, 1317.27, 2758.70, 3.2080, 5.7451],
    [300, 8.58790, 0.001404, 0.021660, 1332.95, 2563.62, 1345.01, 2749.64, 3.2552, 5.7059],
    [305, 9.20943, 0.001425, 0.019934, 1360.18, 2555.85, 1373.30, 2739.43, 3.3028, 5.6657],
    [310, 9.86505, 0.001448, 0.018335, 1387.93, 2547.07, 1402.22, 2727.95, 3.3510, 5.6244],
    [315, 10.55620, 0.001472, 0.016851, 1416.28, 2537.17, 1431.83, 2715.05, 3.3998, 5.5816],
    [320, 11.28430, 0.001499, 0.015471, 1445.31, 2526.01, 1462.22, 2700.59, 3.4494, 5.5372],
    [325, 12.05100, 0.001528, 0.014183, 1475.11, 2513.41, 1493.52, 2684.33, 3.5000, 5.4909],
    [330, 12.85810, 0.001561, 0.012979, 1505.80, 2499.15, 1525.87, 2666.03, 3.5518, 5.4422],
    [335, 13.70730, 0.001597, 0.011847, 1537.56, 2482.96, 1559.45, 2645.35, 3.6050, 5.3906],
    [340, 14.60070, 0.001638, 0.010781, 1570.62, 2464.44, 1594.53, 2621.85, 3.6601, 5.3356],
    [345, 15.54060, 0.001685, 0.009769, 1605.30, 2443.09, 1631.48, 2594.90, 3.7176, 5.2762],
    [350, 16.52940, 0.001740, 0.008802, 1642.13, 2418.14, 1670.89, 2563.64, 3.7784, 5.2110],
    [355, 17.57010, 0.001808, 0.007868, 1681.96, 2388.40, 1713.72, 2526.65, 3.8439, 5.1380],
    [360, 18.66600, 0.001895, 0.006949, 1726.28, 2351.78, 1761.66, 2481.49, 3.9167, 5.0536],
    [365, 19.82140, 0.002017, 0.006012, 1777.79, 2303.79, 1817.77, 2422.95, 4.0014, 4.9497],
    [370, 21.04360, 0.002215, 0.004954, 1844.07, 2230.26, 1890.69, 2334.52, 4.1112, 4.8012],
    [373.95, 22.06400, 0.003106, 0.003106, 2015.73, 2015.73, 2084.26, 2084.26, 4.4070, 4.4070],
])

# Ammonia, saturated (temperature table): Yan Appendix B Table B1 (NIST SRD 69 values)
# columns: T_C, p_MPa, vf, vg (m^3/kg), uf, ug, hf, hg (kJ/kg), sf, sg (kJ/(kg K))
AMMONIA_SAT = np.array([
    [-50, 0.04084, 0.001424, 2.627750, 118.37, 1427.03, 118.43, 1534.34, 0.5661, 6.9112],
    [-45, 0.05449, 0.001436, 2.007080, 140.23, 1433.38, 140.31, 1542.75, 0.6630, 6.8100],
    [-40, 0.07169, 0.001449, 1.553280, 162.22, 1439.56, 162.33, 1550.92, 0.7583, 6.7141],
    [-35, 0.09310, 0.001462, 1.216780, 184.34, 1445.55, 184.48, 1558.83, 0.8522, 6.6232],
    [-30, 0.11943, 0.001475, 0.963963, 206.58, 1451.34, 206.76, 1566.47, 0.9446, 6.5367],
    [-25, 0.15147, 0.001489, 0.771672, 228.94, 1456.92, 229.17, 1573.81, 1.0357, 6.4543],
    [-20, 0.19008, 0.001503, 0.623730, 251.42, 1462.27, 251.71, 1580.83, 1.1254, 6.3757],
    [-15, 0.23617, 0.001518, 0.508679, 274.01, 1467.39, 274.37, 1587.52, 1.2137, 6.3005],
    [-10, 0.29071, 0.001534, 0.418303, 296.72, 1472.25, 297.16, 1593.86, 1.3009, 6.2285],
    [-5, 0.35476, 0.001550, 0.346635, 319.54, 1476.85, 320.09, 1599.82, 1.3868, 6.1593],
    [0, 0.42939, 0.001566, 0.289297, 342.48, 1481.17, 343.16, 1605.39, 1.4716, 6.0926],
    [5, 0.51575, 0.001583, 0.243039, 365.55, 1485.20, 366.36, 1610.55, 1.5553, 6.0284],
    [10, 0.61505, 0.001601, 0.205426, 388.74, 1488.92, 389.72, 1615.27, 1.6380, 5.9662],
    [15, 0.72853, 0.001619, 0.174614, 412.06, 1492.32, 413.24, 1619.53, 1.7197, 5.9060],
    [20, 0.85748, 0.001639, 0.149198, 435.53, 1495.38, 436.94, 1623.31, 1.8005, 5.8475],
    [25, 1.00324, 0.001659, 0.128092, 459.16, 1498.08, 460.82, 1626.59, 1.8804, 5.7904],
    [30, 1.16720, 0.001680, 0.110457, 482.95, 1500.40, 484.91, 1629.32, 1.9597, 5.7347],
    [35, 1.35077, 0.001702, 0.095633, 506.93, 1502.32, 509.23, 1631.49, 2.0382, 5.6802],
    [40, 1.55542, 0.001726, 0.083101, 531.11, 1503.81, 533.79, 1633.06, 2.1161, 5.6265],
    [45, 1.78266, 0.001750, 0.072450, 555.51, 1504.84, 558.63, 1633.99, 2.1936, 5.5736],
    [50, 2.03403, 0.001777, 0.063350, 580.16, 1505.37, 583.77, 1634.23, 2.2706, 5.5213],
    [55, 2.31113, 0.001804, 0.055538, 605.09, 1505.37, 609.26, 1633.73, 2.3473, 5.4693],
    [60, 2.61560, 0.001834, 0.048797, 630.33, 1504.79, 635.12, 1632.42, 2.4239, 5.4174],
    [65, 2.94913, 0.001866, 0.042955, 655.91, 1503.56, 661.42, 1630.24, 2.5004, 5.3655],
    [70, 3.31347, 0.001900, 0.037868, 681.90, 1501.62, 688.20, 1627.09, 2.5770, 5.3131],
    [75, 3.71045, 0.001937, 0.033419, 708.34, 1498.88, 715.53, 1622.88, 2.6539, 5.2601],
    [80, 4.14197, 0.001978, 0.029509, 735.31, 1495.24, 743.50, 1617.46, 2.7312, 5.2060],
    [85, 4.61002, 0.002022, 0.026058, 762.88, 1490.56, 772.20, 1610.69, 2.8093, 5.1505],
    [90, 5.11672, 0.002071, 0.022997, 791.16, 1484.67, 801.76, 1602.34, 2.8884, 5.0929],
    [95, 5.66432, 0.002127, 0.020268, 820.29, 1477.36, 832.34, 1592.16, 2.9689, 5.0328],
    [100, 6.25527, 0.002190, 0.017820, 850.46, 1468.32, 864.16, 1579.79, 3.0513, 4.9691],
    [105, 6.89227, 0.002263, 0.015610, 881.91, 1457.14, 897.51, 1564.72, 3.1363, 4.9007],
    [110, 7.57834, 0.002350, 0.013596, 915.03, 1443.19, 932.84, 1546.23, 3.2249, 4.8258],
    [115, 8.31700, 0.002456, 0.011740, 950.47, 1425.51, 970.89, 1523.15, 3.3190, 4.7418],
    [120, 9.11249, 0.002594, 0.009993, 989.44, 1402.32, 1013.08, 1493.38, 3.4218, 4.6435],
    [125, 9.97022, 0.002795, 0.008283, 1034.97, 1369.69, 1062.84, 1452.28, 3.5417, 4.5199],
    [130, 10.89770, 0.003202, 0.006379, 1100.28, 1312.96, 1135.17, 1382.47, 3.7153, 4.3287],
    [132.25, 11.33930, 0.004444, 0.004444, 1211.98, 1211.98, 1262.38, 1262.38, 4.0258, 4.0258],
])

# Water, superheated vapour: Yan Appendix A Table A2 (NIST SRD 69 values), four of its pressures
# {p_kPa: rows of [T_C, v (m^3/kg), u (kJ/kg), h (kJ/kg), s (kJ/(kg K))]}; the 'Sat' rows are omitted
WATER_SUPERHEATED = {
    10: np.array([
        [50, 14.8670, 2443.30, 2591.97, 8.1741],
        [100, 17.1964, 2515.49, 2687.46, 8.4489],
        [150, 19.5132, 2587.91, 2783.04, 8.6892],
        [200, 21.8256, 2661.34, 2879.60, 8.9049],
        [250, 24.1361, 2736.09, 2977.45, 9.1015],
        [300, 26.4456, 2812.28, 3076.74, 9.2827],
        [350, 28.7545, 2890.00, 3177.55, 9.4513],
        [400, 31.0631, 2969.32, 3279.95, 9.6094],
        [450, 33.3714, 3050.27, 3383.98, 9.7584],
        [500, 35.6796, 3132.89, 3489.69, 9.8998],
        [550, 37.9876, 3217.24, 3597.11, 10.0344],
        [600, 40.2956, 3303.32, 3706.27, 10.1631],
        [650, 42.6035, 3391.17, 3817.20, 10.2866],
        [700, 44.9113, 3480.79, 3929.90, 10.4055],
        [750, 47.2191, 3572.18, 4044.37, 10.5202],
        [800, 49.5269, 3665.35, 4160.62, 10.6311],
        [850, 51.8347, 3760.27, 4278.61, 10.7386],
        [900, 54.1424, 3856.91, 4398.33, 10.8429],
        [950, 56.4501, 3955.25, 4519.75, 10.9442],
        [1000, 58.7578, 4055.24, 4642.82, 11.0428],
    ]),
    100: np.array([
        [150, 1.93665, 2582.94, 2776.60, 7.6148],
        [200, 2.17243, 2658.21, 2875.45, 7.8356],
        [250, 2.40615, 2733.90, 2974.51, 8.0346],
        [300, 2.63885, 2810.65, 3074.53, 8.2172],
        [350, 2.87095, 2888.73, 3175.83, 8.3866],
        [400, 3.10269, 2968.29, 3278.56, 8.5452],
        [450, 3.33420, 3049.42, 3382.84, 8.6946],
        [500, 3.56553, 3132.18, 3488.73, 8.8362],
        [550, 3.79675, 3216.62, 3596.30, 8.9710],
        [600, 4.02788, 3302.79, 3705.58, 9.0998],
        [650, 4.25895, 3390.70, 3816.60, 9.2235],
        [700, 4.48997, 3480.37, 3929.37, 9.3424],
        [750, 4.72094, 3571.82, 4043.91, 9.4572],
        [800, 4.95189, 3665.02, 4160.21, 9.5682],
        [850, 5.18280, 3759.97, 4278.25, 9.6757],
        [900, 5.41370, 3856.64, 4398.01, 9.7800],
        [950, 5.64458, 3955.00, 4519.46, 9.8813],
        [1000, 5.87545, 4055.01, 4642.56, 9.9800],
    ]),
    1000: np.array([
        [200, 0.206022, 2622.24, 2828.26, 6.6955],
        [250, 0.232747, 2710.38, 2943.12, 6.9265],
        [300, 0.257988, 2793.64, 3051.63, 7.1246],
        [350, 0.282504, 2875.68, 3158.18, 7.3029],
        [400, 0.306606, 2957.87, 3264.47, 7.4669],
        [450, 0.330447, 3040.85, 3371.30, 7.6200],
        [500, 0.354113, 3124.99, 3479.10, 7.7642],
        [550, 0.377655, 3210.48, 3588.14, 7.9008],
        [600, 0.401108, 3297.47, 3698.58, 8.0310],
        [650, 0.424494, 3386.04, 3810.53, 8.1557],
        [700, 0.447828, 3476.25, 3924.08, 8.2755],
        [750, 0.471121, 3568.14, 4039.26, 8.3909],
        [800, 0.494382, 3661.71, 4156.09, 8.5024],
        [850, 0.517618, 3756.98, 4274.60, 8.6103],
        [900, 0.540832, 3853.92, 4394.75, 8.7150],
        [950, 0.564028, 3952.52, 4516.54, 8.8166],
        [1000, 0.587210, 4052.74, 4639.95, 8.9155],
    ]),
    4000: np.array([
        [300, 0.058870, 2726.21, 2961.69, 6.3639],
        [350, 0.066473, 2827.43, 3093.32, 6.5843],
        [400, 0.073431, 2920.75, 3214.47, 6.7714],
        [450, 0.080043, 3011.00, 3331.17, 6.9386],
        [500, 0.086442, 3100.25, 3446.02, 7.0922],
        [550, 0.092700, 3189.54, 3560.34, 7.2355],
        [600, 0.098859, 3279.43, 3674.87, 7.3705],
        [650, 0.104945, 3370.30, 3790.08, 7.4988],
        [700, 0.110976, 3462.37, 3906.28, 7.6214],
        [750, 0.116966, 3555.78, 4023.64, 7.7390],
        [800, 0.122922, 3650.63, 4142.32, 7.8523],
        [850, 0.128852, 3746.97, 4262.38, 7.9616],
        [900, 0.134759, 3844.83, 4383.87, 8.0674],
        [950, 0.140649, 3944.22, 4506.82, 8.1701],
        [1000, 0.146525, 4045.13, 4631.23, 8.2697],
    ]),
}

# Ammonia, superheated vapour: Yan Appendix B Table B2 (NIST SRD 69 values), four of its pressures, T <= 200 C
# {p_kPa: rows of [T_C, v (m^3/kg), u (kJ/kg), h (kJ/kg), s (kJ/(kg K))]}; the 'Sat' rows are omitted
AMMONIA_SUPERHEATED = {
    100: np.array([
        [-30, 1.15730, 1453.48, 1569.21, 6.6323],
        [-20, 1.21016, 1470.67, 1591.69, 6.7229],
        [-10, 1.26220, 1487.57, 1613.79, 6.8085],
        [0, 1.31365, 1504.29, 1635.65, 6.8900],
        [10, 1.36463, 1520.91, 1657.37, 6.9681],
        [20, 1.41527, 1537.50, 1679.03, 7.0433],
        [30, 1.46562, 1554.10, 1700.66, 7.1159],
        [40, 1.51576, 1570.75, 1722.33, 7.1862],
        [50, 1.56571, 1587.48, 1744.05, 7.2544],
        [60, 1.61552, 1604.30, 1765.85, 7.3209],
        [70, 1.66521, 1621.23, 1787.76, 7.3857],
        [80, 1.71480, 1638.30, 1809.78, 7.4490],
        [90, 1.76431, 1655.51, 1831.94, 7.5108],
        [100, 1.81374, 1672.88, 1854.25, 7.5714],
        [110, 1.86311, 1690.40, 1876.71, 7.6308],
        [120, 1.91243, 1708.10, 1899.34, 7.6891],
        [130, 1.96170, 1725.97, 1922.14, 7.7464],
        [140, 2.01094, 1744.02, 1945.12, 7.8027],
        [150, 2.06013, 1762.27, 1968.28, 7.8581],
        [160, 2.10930, 1780.70, 1991.63, 7.9126],
        [170, 2.15844, 1799.34, 2015.18, 7.9664],
        [180, 2.20755, 1818.17, 2038.92, 8.0194],
        [190, 2.25664, 1837.21, 2062.87, 8.0716],
        [200, 2.30572, 1856.45, 2087.02, 8.1232],
    ]),
    500: np.array([
        [10, 0.257544, 1496.73, 1625.50, 6.0957],
        [20, 0.269459, 1516.75, 1651.48, 6.1859],
        [30, 0.280985, 1536.04, 1676.54, 6.2699],
        [40, 0.292219, 1554.85, 1700.96, 6.3492],
        [50, 0.303228, 1573.34, 1724.95, 6.4246],
        [60, 0.314060, 1591.62, 1748.65, 6.4968],
        [70, 0.324751, 1609.77, 1772.15, 6.5663],
        [80, 0.335327, 1627.88, 1795.54, 6.6335],
        [90, 0.345807, 1645.98, 1818.88, 6.6987],
        [100, 0.356208, 1664.10, 1842.21, 6.7621],
        [110, 0.366541, 1682.30, 1865.57, 6.8239],
        [120, 0.376818, 1700.58, 1888.99, 6.8842],
        [130, 0.387045, 1718.97, 1912.49, 6.9432],
        [140, 0.397230, 1737.48, 1936.10, 7.0011],
        [150, 0.407378, 1756.14, 1959.83, 7.0578],
        [160, 0.417494, 1774.94, 1983.69, 7.1136],
        [170, 0.427580, 1793.91, 2007.70, 7.1683],
        [180, 0.437642, 1813.04, 2031.86, 7.2223],
        [190, 0.447681, 1832.35, 2056.19, 7.2754],
        [200, 0.457700, 1851.85, 2080.70, 7.3277],
    ]),
    1000: np.array([
        [30, 0.132037, 1510.15, 1642.19, 5.8437],
        [40, 0.138657, 1532.67, 1671.32, 5.9383],
        [50, 0.144962, 1554.01, 1698.97, 6.0252],
        [60, 0.151035, 1574.55, 1725.59, 6.1063],
        [70, 0.156931, 1594.54, 1751.47, 6.1829],
        [80, 0.162688, 1614.15, 1776.84, 6.2558],
        [90, 0.168333, 1633.51, 1801.85, 6.3256],
        [100, 0.173887, 1652.71, 1826.60, 6.3929],
        [110, 0.179366, 1671.83, 1851.19, 6.4579],
        [120, 0.184782, 1690.91, 1875.69, 6.5210],
        [130, 0.190144, 1709.99, 1900.14, 6.5824],
        [140, 0.195460, 1729.12, 1924.58, 6.6423],
        [150, 0.200737, 1748.32, 1949.05, 6.7008],
        [160, 0.205979, 1767.61, 1973.59, 6.7581],
        [170, 0.211191, 1787.01, 1998.20, 6.8143],
        [180, 0.216376, 1806.54, 2022.91, 6.8695],
        [190, 0.221538, 1826.21, 2047.74, 6.9237],
        [200, 0.226678, 1846.02, 2072.70, 6.9770],
    ]),
    2000: np.array([
        [50, 0.064729, 1507.24, 1636.69, 5.5357],
        [60, 0.068754, 1534.87, 1672.38, 5.6445],
        [70, 0.072463, 1560.13, 1705.06, 5.7411],
        [80, 0.075952, 1583.81, 1735.72, 5.8292],
        [90, 0.079278, 1606.42, 1764.97, 5.9109],
        [100, 0.082480, 1628.26, 1793.22, 5.9876],
        [110, 0.085583, 1649.58, 1820.75, 6.0604],
        [120, 0.088608, 1670.52, 1847.73, 6.1300],
        [130, 0.091567, 1691.20, 1874.33, 6.1968],
        [140, 0.094473, 1711.71, 1900.65, 6.2613],
        [150, 0.097332, 1732.11, 1926.78, 6.3237],
        [160, 0.100152, 1752.47, 1952.77, 6.3845],
        [170, 0.102939, 1772.81, 1978.69, 6.4436],
        [180, 0.105696, 1793.19, 2004.58, 6.5014],
        [190, 0.108427, 1813.62, 2030.47, 6.5579],
        [200, 0.111136, 1834.13, 2056.40, 6.6133],
    ]),
}


# NASA Glenn 9-coefficient polynomial for H2O(g), from data/thermo.inp of NASA CEA (Apache-2.0;
# McBride, Zehe & Gordon, NASA/TP-2002-211556): (T_low, T_high, [a1..a7], [b1, b2]).
NASA_H2O = [
    (200.0, 1000.0, [-39479.6083, 575.573102, 0.931782653, 0.00722271286, -7.34255737e-06, 4.95504349e-09, -1.336933246e-12], [-33039.7431, 17.24205775]),
    (1000.0, 6000.0, [1034972.096, -2412.698562, 4.64611078, 0.002291998307, -6.83683048e-07, 9.42646893e-11, -4.82238053e-15], [-13842.86509, -7.97814851]),
]

FLUIDS = {
    "water": {"sat": WATER_SAT, "sup": WATER_SUPERHEATED, "M": M_WATER},
    "ammonia": {"sat": AMMONIA_SAT, "sup": AMMONIA_SUPERHEATED, "M": M_AMMONIA},
}
SAT_COLUMNS = ["T", "p", "vf", "vg", "uf", "ug", "hf", "hg", "sf", "sg"]


# ---------------------------------------------------------------- given helpers

def sat_table(fluid: str) -> dict[str, np.ndarray]:
    """The saturation table of a fluid as named columns, with p converted to kPa."""
    tab = FLUIDS[fluid]["sat"]
    cols = {name: tab[:, i].astype(float) for i, name in enumerate(SAT_COLUMNS)}
    cols["p"] = cols["p"] * 1000.0
    return cols


def critical_point(fluid: str) -> tuple[float, float]:
    """(T_c in deg C, p_c in kPa): the last row of the saturation table, where v_f = v_g."""
    t = sat_table(fluid)
    return float(t["T"][-1]), float(t["p"][-1])


def gas_constant(fluid: str) -> float:
    """Specific gas constant R = R_U / M in kJ/(kg K)."""
    return R_U / FLUIDS[fluid]["M"] / 1000.0


def interp(x: float, xp: np.ndarray, fp: np.ndarray) -> float:
    """Linear interpolation that refuses to extrapolate (raises ValueError outside [xp[0], xp[-1]])."""
    if x < xp[0] - 1e-12 or x > xp[-1] + 1e-12:
        raise ValueError(f"{x} is outside the table range [{xp[0]}, {xp[-1]}]")
    return float(np.interp(x, xp, fp))


def nasa_h2o_enthalpy(T_K: float) -> float:
    """Ideal-gas enthalpy of H2O [kJ/kg] from the NASA polynomial (formation-enthalpy datum)."""
    for lo, hi, a, b in NASA_H2O:
        if lo <= T_K <= hi:
            T = T_K
            hrt = (-a[0] / T ** 2 + a[1] * np.log(T) / T + a[2] + a[3] * T / 2 + a[4] * T ** 2 / 3
                   + a[5] * T ** 3 / 4 + a[6] * T ** 4 / 5 + b[0] / T)
            return float(R_U * T * hrt / M_WATER / 1000.0)
    raise ValueError("T outside the NASA range")


# ---------------------------------------------------------------- saturation

def saturation_pressure(fluid: str, T_C: float, method: str = "log") -> float:
    """Saturation pressure [kPa] at T_C. method "linear": p linear in T between table rows; "log":
    ln p linear in 1/T (the shape Clausius-Clapeyron predicts), much more accurate between rows."""
    t = sat_table(fluid)
    if method == "linear":
        return interp(T_C, t["T"], t["p"])
    if method == "log":
        interp(T_C, t["T"], t["T"])                                   # range check
        i = int(np.clip(np.searchsorted(t["T"], T_C) - 1, 0, len(t["T"]) - 2))
        x0, x1 = 1.0 / (t["T"][i] + T_OFFSET), 1.0 / (t["T"][i + 1] + T_OFFSET)
        y0, y1 = np.log(t["p"][i]), np.log(t["p"][i + 1])
        x = 1.0 / (T_C + T_OFFSET)
        return float(np.exp(y0 + (y1 - y0) * (x - x0) / (x1 - x0)))
    raise ValueError(method)


def saturation_temperature(fluid: str, p_kPa: float) -> float:
    """Saturation temperature [deg C] at p_kPa: the inverse of saturation_pressure(..., "log")."""
    t = sat_table(fluid)
    if not t["p"][0] <= p_kPa <= t["p"][-1]:
        raise ValueError(f"p = {p_kPa} kPa is outside the saturation table")
    return float(brentq(lambda T: saturation_pressure(fluid, T) - p_kPa, t["T"][0], t["T"][-1], xtol=1e-10))


def saturation_state(fluid: str, T_C: float) -> dict[str, float]:
    """All saturation properties at T_C: p (log interpolation) and vf, vg, uf, ug, hf, hg, sf, sg (linear
    in T), plus the differences v_fg, u_fg, h_fg, s_fg."""
    t = sat_table(fluid)
    st = {"T": float(T_C), "p": saturation_pressure(fluid, T_C)}
    for k in SAT_COLUMNS[2:]:
        st[k] = interp(T_C, t["T"], t[k])
    for y in "vuhs":
        st[y + "fg"] = st[y + "g"] - st[y + "f"]
    return st


def quality(fluid: str, T_C: float, prop: str, value: float) -> float:
    """Quality x = (y - y_f) / (y_g - y_f) for prop in "v", "u", "h", "s" at saturation temperature T_C.
    Values outside [0, 1] mean the state is not a two-phase mixture."""
    st = saturation_state(fluid, T_C)
    return (value - st[prop + "f"]) / (st[prop + "g"] - st[prop + "f"])


def mixture_property(fluid: str, T_C: float, x: float, prop: str) -> float:
    """y = y_f + x (y_g - y_f) for a two-phase mixture of quality x (0 <= x <= 1)."""
    if not 0.0 <= x <= 1.0:
        raise ValueError("quality must lie in [0, 1]")
    st = saturation_state(fluid, T_C)
    return st[prop + "f"] + x * (st[prop + "g"] - st[prop + "f"])


def phase_from_Tp(fluid: str, T_C: float, p_kPa: float, rel_tol: float = 1e-6) -> str:
    """'compressed liquid' if p > p_sat(T), 'superheated vapor' if p < p_sat(T), 'saturated' if equal within
    rel_tol (then T and p are not independent and do not fix the state), 'supercritical' above T_c."""
    Tc, pc = critical_point(fluid)
    if T_C > Tc:
        return "supercritical"
    ps = saturation_pressure(fluid, T_C)
    if abs(p_kPa - ps) <= rel_tol * ps:
        return "saturated"
    return "compressed liquid" if p_kPa > ps else "superheated vapor"


def phase_from_Tv(fluid: str, T_C: float, v: float) -> tuple[str, float]:
    """(phase, x) from temperature and specific volume: compressed liquid (v < v_f, x = NaN), saturated
    mixture (v_f <= v <= v_g, x from v), superheated vapor (v > v_g, x = NaN), or supercritical (T > T_c)."""
    if T_C > critical_point(fluid)[0]:
        return "supercritical", float("nan")
    st = saturation_state(fluid, T_C)
    if v < st["vf"]:
        return "compressed liquid", float("nan")
    if v > st["vg"]:
        return "superheated vapor", float("nan")
    return "saturated mixture", (v - st["vf"]) / st["vfg"]


# ---------------------------------------------------------------- superheated vapour

def _block(fluid: str, p_kPa: float, T_C: float) -> np.ndarray:
    """[v, u, h, s] at T_C on one tabulated pressure, linear in T."""
    rows = FLUIDS[fluid]["sup"][p_kPa]
    return np.array([interp(T_C, rows[:, 0], rows[:, j]) for j in range(1, 5)])


def superheated_state(fluid: str, p_kPa: float, T_C: float) -> dict[str, float]:
    """Superheated vapour properties at (p, T): linear interpolation in T on the two tabulated pressures
    that bracket p, then between them with p v linear in p (exact for an ideal gas), u and h linear in p,
    and s linear in ln p (exact for an ideal gas). Raises ValueError outside the table."""
    ps = np.array(sorted(FLUIDS[fluid]["sup"]))
    if not ps[0] <= p_kPa <= ps[-1]:
        raise ValueError(f"p = {p_kPa} kPa is outside the superheated table [{ps[0]}, {ps[-1]}]")
    if p_kPa in FLUIDS[fluid]["sup"]:                               # on a tabulated pressure: that block only
        a = _block(fluid, p_kPa, T_C)
        return {"T": float(T_C), "p": float(p_kPa), "v": float(a[0]), "u": float(a[1]), "h": float(a[2]), "s": float(a[3])}
    i = int(np.clip(np.searchsorted(ps, p_kPa) - 1, 0, len(ps) - 2))
    p0, p1 = ps[i], ps[i + 1]
    a, b = _block(fluid, p0, T_C), _block(fluid, p1, T_C)
    f = (p_kPa - p0) / (p1 - p0)
    pv = (1 - f) * p0 * a[0] + f * p1 * b[0]
    g = np.log(p_kPa / p0) / np.log(p1 / p0)
    return {"T": float(T_C), "p": float(p_kPa), "v": float(pv / p_kPa), "u": float((1 - f) * a[1] + f * b[1]),
            "h": float((1 - f) * a[2] + f * b[2]), "s": float((1 - g) * a[3] + g * b[3])}


def state_from_Tv(fluid: str, T_C: float, v: float) -> dict[str, float]:
    """The complete state (phase, x, p, v, u, h, s) from (T, v). In the dome, p = p_sat and the rest
    follow from the quality; in the superheated (or supercritical) vapour region p is found by root
    finding on v(p, T) within the superheated table. Raises ValueError for a compressed liquid or a state
    outside the table."""
    phase, x = phase_from_Tv(fluid, T_C, v)
    if phase == "saturated mixture":
        st = saturation_state(fluid, T_C)
        return {"phase": phase, "x": x, "T": float(T_C), "p": st["p"], "v": float(v),
                "u": st["uf"] + x * st["ufg"], "h": st["hf"] + x * st["hfg"], "s": st["sf"] + x * st["sfg"]}
    if phase in ("superheated vapor", "supercritical"):
        ps = np.array(sorted(FLUIDS[fluid]["sup"]), dtype=float)
        lo = ps[0]
        hi = ps[-1] if phase == "supercritical" else min(ps[-1], saturation_pressure(fluid, T_C))
        p = brentq(lambda pp: superheated_state(fluid, pp, T_C)["v"] - v, lo, hi, xtol=1e-10)
        st = superheated_state(fluid, p, T_C)
        return {"phase": phase, "x": float("nan"), **st}
    raise ValueError("compressed liquid: this lookup has no compressed-liquid table")


def rigid_tank_heating(fluid: str, V: float, m: float, T1_C: float, T2_C: float) -> dict[str, float]:
    """A sealed rigid tank (closed system, no work) of volume V [m^3] holding m [kg] is taken from T1 to T2.
    v = V/m is fixed, so both states follow from (T, v). First law: Q = m (u2 - u1) [kJ]. Returns the two
    states' pressures, qualities and phases, the liquid volume fractions, and Q."""
    v = V / m
    s1, s2 = state_from_Tv(fluid, T1_C, v), state_from_Tv(fluid, T2_C, v)
    def liquid_fraction(st, T):
        if st["phase"] != "saturated mixture":
            return 0.0
        return (1.0 - st["x"]) * saturation_state(fluid, T)["vf"] / v
    return {"v": v, "p1": s1["p"], "p2": s2["p"], "x1": s1["x"], "x2": s2["x"], "phase1": s1["phase"],
            "phase2": s2["phase"], "liquid_volume_fraction1": liquid_fraction(s1, T1_C),
            "liquid_volume_fraction2": liquid_fraction(s2, T2_C), "Q": m * (s2["u"] - s1["u"])}


def max_liquid_fill(fluid: str, T_fill_C: float, T_max_C: float) -> float:
    """Largest liquid volume fraction at filling temperature T_fill such that a sealed rigid tank is still
    two-phase (not liquid-full) when it warms to T_max: the specific volume must be at least v_f(T_max)."""
    v = saturation_state(fluid, T_max_C)["vf"]
    st = saturation_state(fluid, T_fill_C)
    x = (v - st["vf"]) / st["vfg"]
    return (1.0 - x) * st["vf"] / v


# ---------------------------------------------------------------- checks of the table itself

def clausius_clapeyron_hfg(fluid: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """At every interior row of the saturation table, h_fg predicted by Clausius-Clapeyron,
    h_fg = T v_fg dp/dT, with dp/dT = -(p / T^2) d(ln p)/d(1/T) from central differences of ln p
    against 1/T across the neighbouring rows. Returns (T_C, h_fg predicted, h_fg tabulated) in kJ/kg."""
    t = sat_table(fluid)
    TK = t["T"] + T_OFFSET
    x, y = 1.0 / TK, np.log(t["p"])
    dydx = (y[2:] - y[:-2]) / (x[2:] - x[:-2])
    T, p = TK[1:-1], t["p"][1:-1]
    dpdT = -p / T ** 2 * dydx
    vfg = (t["vg"] - t["vf"])[1:-1]
    return t["T"][1:-1], T * vfg * dpdT, (t["hg"] - t["hf"])[1:-1]


def enthalpy_identity_residual(fluid: str) -> tuple[np.ndarray, np.ndarray]:
    """Residuals h - (u + p v) [kJ/kg] for every saturated-liquid, saturated-vapour and superheated entry,
    with the bound that table round-off alone allows: half a unit in the last printed digit of h, u, p and
    v (h, u: 0.005 kJ/kg each; p_sat: 0.005 kPa; v: 5e-7 m^3/kg for v_f, 5e-6 relative for v_g and
    superheated v, which are printed to six significant figures). Returns (residuals, bounds)."""
    t = sat_table(fluid)
    res, bnd = [], []
    for v, u, h, dv in ((t["vf"], t["uf"], t["hf"], 5e-7 + 0 * t["vf"]), (t["vg"], t["ug"], t["hg"], 5e-6 * t["vg"])):
        res.append(h - u - t["p"] * v)
        bnd.append(0.01 + 0.005 * v + t["p"] * dv)
    for p, rows in FLUIDS[fluid]["sup"].items():
        v, u, h = rows[:, 1], rows[:, 2], rows[:, 3]
        res.append(h - u - p * v)
        bnd.append(0.01 + p * 5e-6 * v)
    return np.concatenate(res), np.concatenate(bnd)


def compressibility_factor(fluid: str, p_kPa: float, T_C: float, v: float) -> float:
    """Z = p v / (R T): 1 for an ideal gas."""
    return p_kPa * v / (gas_constant(fluid) * (T_C + T_OFFSET))


def superheated_dh_vs_ideal_gas(p_kPa: float, T1_C: float, T2_C: float) -> tuple[float, float]:
    """(table h(T2) - h(T1), NASA ideal-gas h(T2) - h(T1)) for superheated steam at p_kPa [kJ/kg]."""
    tab = superheated_state("water", p_kPa, T2_C)["h"] - superheated_state("water", p_kPa, T1_C)["h"]
    ig = nasa_h2o_enthalpy(T2_C + T_OFFSET) - nasa_h2o_enthalpy(T1_C + T_OFFSET)
    return float(tab), float(ig)


def leave_one_out_psat(fluid: str, method: str) -> tuple[np.ndarray, np.ndarray]:
    """Drop every odd interior row of the saturation table, interpolate p_sat there from the remaining
    rows with `method`, and return (T_C of the dropped rows, relative error against the dropped values)."""
    t = sat_table(fluid)
    keep = np.ones(len(t["T"]), bool)
    keep[1:-1:2] = False
    Tk, pk = t["T"][keep], t["p"][keep]
    errs = []
    for T, p in zip(t["T"][~keep], t["p"][~keep]):
        i = np.searchsorted(Tk, T) - 1
        if method == "linear":
            pe = pk[i] + (pk[i + 1] - pk[i]) * (T - Tk[i]) / (Tk[i + 1] - Tk[i])
        else:
            x0, x1, x = 1 / (Tk[i] + T_OFFSET), 1 / (Tk[i + 1] + T_OFFSET), 1 / (T + T_OFFSET)
            pe = np.exp(np.log(pk[i]) + (np.log(pk[i + 1]) - np.log(pk[i])) * (x - x0) / (x1 - x0))
        errs.append(pe / p - 1.0)
    return t["T"][~keep], np.array(errs)


# ---------------------------------------------------------------- driver

def main(out_dir: str | Path = ".") -> None:
    """Print the data-based investigation and the engineering answer; write property_tables.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    print("Data: NIST SRD 69 values via Yan (LibreTexts) Appendices A (water) and B (ammonia).")

    print("\n1. Clausius-Clapeyron: h_fg = T v_fg dp/dT from the table's own p_sat(T) against tabulated h_fg")
    for fl in ("water", "ammonia"):
        T, hcc, htab = clausius_clapeyron_hfg(fl)
        err = hcc / htab - 1
        mask = T <= (critical_point(fl)[0] - 30)
        print(f"  {fl}: {len(T)} interior rows; median |error| {100*np.median(np.abs(err)):.3f} %, "
              f"largest |error| below Tc - 30 C {100*np.max(np.abs(err[mask])):.2f} % (at {T[mask][np.argmax(np.abs(err[mask]))]:.0f} C), "
              f"nearest the critical point {100*err[-1]:+.1f} % (at {T[-1]:.2f} C)")
        for Tq in ((20.0, 100.0, 200.0, 300.0) if fl == "water" else (-30.0, 0.0, 40.0, 100.0)):
            j = int(np.argmin(np.abs(T - Tq)))
            print(f"    {T[j]:7.2f} C: CC {hcc[j]:8.2f}, table {htab[j]:8.2f} kJ/kg ({100*err[j]:+.3f} %)")

    print("\n2. h = u + p v, residual against what table round-off allows")
    for fl in ("water", "ammonia"):
        r, b = enthalpy_identity_residual(fl)
        print(f"  {fl}: {len(r)} entries; largest |residual| {np.max(np.abs(r)):.3f} kJ/kg; largest |residual|/bound "
              f"{np.max(np.abs(r)/b):.2f}; entries exceeding the bound: {int(np.sum(np.abs(r) > b))}")

    print("\n3. Compressibility factor Z = p v / (R T) of superheated vapour")
    for fl, Ts in (("water", (300.0, 500.0)), ("ammonia", (60.0, 200.0))):
        for T in Ts:
            zs = []
            for p in sorted(FLUIDS[fl]["sup"]):
                rows = FLUIDS[fl]["sup"][p]
                if rows[0, 0] <= T <= rows[-1, 0]:
                    zs.append(f"{p} kPa: {compressibility_factor(fl, p, T, superheated_state(fl, p, T)['v']):.4f}")
            print(f"  {fl} at {T:.0f} C: " + "; ".join(zs))

    print("\n4. Interpolating p_sat between rows: leave every other row out")
    for fl in ("water", "ammonia"):
        for meth in ("linear", "log"):
            T, e = leave_one_out_psat(fl, meth)
            print(f"  {fl}, {meth:>6}: largest |error| {100*np.max(np.abs(e)):.3f} % (at {T[np.argmax(np.abs(e))]:.0f} C), "
                  f"median {100*np.median(np.abs(e)):.3f} %")

    print("\n5. Superheated steam at 10 kPa against the NASA ideal-gas H2O fit: h(T2) - h(T1), kJ/kg")
    for T1, T2 in ((100.0, 300.0), (100.0, 600.0), (300.0, 1000.0)):
        tab, ig = superheated_dh_vs_ideal_gas(10, T1, T2)
        print(f"  {T1:.0f} -> {T2:.0f} C: table {tab:.2f}, NASA {ig:.2f}, difference {100*(tab/ig-1):+.3f} %")
    tab, ig = superheated_dh_vs_ideal_gas(4000, 300.0, 600.0)
    print(f"  at 4000 kPa, 300 -> 600 C: table {tab:.2f}, NASA {ig:.2f}, difference {100*(tab/ig-1):+.2f} % (non-ideal)")

    print("\nEngineering case (typical, illustrative): a sealed ammonia propellant tank filled at 15 C")
    Tf, Tm = 15.0, 40.0
    st = saturation_state("ammonia", Tf)
    print(f"  p_sat(15 C) = {st['p']:.1f} kPa: the tank must hold at least this to keep ammonia liquid; "
          f"p_sat(40 C) = {saturation_pressure('ammonia', Tm):.1f} kPa")
    fill = 0.90
    V = 1.0
    m = fill * V / st["vf"] + (1 - fill) * V / st["vg"]
    r = rigid_tank_heating("ammonia", V, m, Tf, Tm)
    print(f"  1 m^3 tank, 90 % liquid by volume at 15 C: m = {m:.2f} kg, x1 = {r['x1']:.5f} "
          f"(vapour is {100*r['x1']:.3f} % of the mass but {100*(1-r['liquid_volume_fraction1']):.0f} % of the volume)")
    print(f"  warmed to 40 C: p2 = {r['p2']:.1f} kPa, x2 = {r['x2']:.5f}, liquid {100*r['liquid_volume_fraction2']:.2f} % "
          f"of the volume; heat in Q = m (u2 - u1) = {r['Q']/1e3:.1f} MJ")
    print(f"  largest liquid fill at 15 C that stays two-phase up to 40 C: {100*max_liquid_fill('ammonia', Tf, Tm):.2f} %; "
          f"up to 50 C: {100*max_liquid_fill('ammonia', Tf, 50.0):.2f} %")
    stw = saturation_state("water", 5.0)
    print(f"  water evaporated at 5 C ({stw['p']:.3f} kPa) from a 20 C feed absorbs h_g(5 C) - h_f(20 C) = "
          f"{stw['hg'] - saturation_state('water', 20.0)['hf']:.2f} kJ/kg; ammonia boiled at -30 C from 15 C liquid: "
          f"{saturation_state('ammonia', -30.0)['hg'] - st['hf']:.2f} kJ/kg")

    fig, axs = plt.subplots(1, 3, figsize=(14, 4.2))
    for fl, c in (("water", "C0"), ("ammonia", "C1")):
        t = sat_table(fl)
        axs[0].semilogx(t["vf"], t["T"], c, label=f"{fl} saturated liquid")
        axs[0].semilogx(t["vg"], t["T"], c + "--", label=f"{fl} saturated vapour")
        T, hcc, htab = clausius_clapeyron_hfg(fl)
        axs[1].plot(T, 100 * (hcc / htab - 1), c + ".-", label=fl)
        for meth, ls in (("linear", "x"), ("log", "o")):
            Tl, e = leave_one_out_psat(fl, meth)
            axs[2].semilogy(Tl, np.abs(e) * 100 + 1e-6, c + ls, ms=4, label=f"{fl}, {meth}")
    axs[0].set_xlabel("v (m^3/kg)")
    axs[0].set_ylabel("T (C)")
    axs[0].set_title("saturation domes (NIST SRD 69 via Yan)")
    axs[0].legend(fontsize=7)
    axs[1].set_ylim(-3, 3)
    axs[1].set_xlabel("T (C)")
    axs[1].set_ylabel("Clausius-Clapeyron h_fg error (%)")
    axs[1].set_title("dp/dT = h_fg / (T v_fg) from the table")
    axs[1].legend(fontsize=8)
    axs[2].set_xlabel("T (C)")
    axs[2].set_ylabel("|p_sat interpolation error| (%)")
    axs[2].set_title("leave-one-out: linear vs ln p in 1/T")
    axs[2].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out_dir / "property_tables.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main()
