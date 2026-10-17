"""Public tests for P07.14 Heat Transfer: Radiation. They assert properties, not magic numbers:
Planck's law integrates to sigma T^4 and peaks at lambda T = 2898 um K; band fractions lie in [0, 1]
and depend on lambda T only; every closed-form view factor obeys reciprocity A_i F_ij = A_j F_ji and
its geometric limits; the six faces of a box built from two closed forms satisfy the summation rule
sum_j F_ij = 1; numerical double-area integration reproduces the closed forms, at second order for
parallel plates and first order when the surfaces share an edge; the crossed-string rule reproduces
the two-dimensional closed forms; the radiosity solver's net heat flows sum to zero, reduce to the
black-body exchange when every emittance is one and to the parallel-plate formula for two surfaces,
and give zero net flow in an isothermal enclosure; a reradiating surface settles between the hot and
cold temperatures; and the radiator's energy balance closes.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.integrate import quad

from conftest import load_impl

impl = load_impl(__file__)
SIG = 5.670374419e-8


# ---------------------------------------------------------------- black body

@pytest.mark.parametrize("T", [300.0, 1000.0, 5772.0])
def test_planck_integrates_to_sigma_T4_and_peaks_at_wien(T):
    lam_pk = 2.8978e-3 / T
    tot = quad(lambda lam: float(impl.planck_spectral(lam, T)), lam_pk / 50, lam_pk * 2000, limit=500,
               points=[lam_pk])[0]
    assert tot == pytest.approx(SIG * T ** 4, rel=1e-6)
    lam = np.linspace(0.5 * lam_pk, 2 * lam_pk, 20001)
    assert lam[np.argmax(impl.planck_spectral(lam, T))] * T == pytest.approx(2.8978e-3, rel=2e-4)
    assert impl.blackbody_emissive_power(T) == pytest.approx(SIG * T ** 4)


def test_band_fractions_are_fractions_and_scale_with_lambda_T():
    assert impl.blackbody_fraction(1e-9, np.inf, 1000.0) == pytest.approx(1.0, abs=1e-10)
    # the fraction below lambda T = 2898 um K (the Wien peak) is 25.0 %
    assert impl.blackbody_fraction(1e-9, 2.8978e-3 / 1000.0, 1000.0) == pytest.approx(0.2501, abs=2e-4)
    f1 = impl.blackbody_fraction(1e-6, 3e-6, 1000.0)
    f2 = impl.blackbody_fraction(2e-6, 6e-6, 500.0)                 # same lambda T band
    assert 0 < f1 < 1 and f1 == pytest.approx(f2, rel=1e-8)
    # sunlight lies mostly below 3 um, room-temperature emission almost entirely above
    assert impl.blackbody_fraction(1e-9, 3e-6, impl.T_SUN) > 0.95
    assert impl.blackbody_fraction(1e-9, 3e-6, 300.0) < 1e-3


# ---------------------------------------------------------------- closed forms

def test_parallel_rectangles_limits():
    assert impl.vf_parallel_rectangles(100.0, 100.0, 0.01) > 0.999            # large, close plates see only each other
    a, b, c = 0.01, 0.02, 10.0                                                 # small, far: F -> A2 / (pi c^2)
    assert impl.vf_parallel_rectangles(a, b, c) == pytest.approx(a * b / (np.pi * c ** 2), rel=1e-4)
    # moving the plates apart lowers the view factor
    Fs = [impl.vf_parallel_rectangles(1.0, 0.8, c) for c in (0.1, 0.5, 1.0, 3.0)]
    assert np.all(np.diff(Fs) < 0)


def test_perpendicular_rectangles_reciprocity_and_limits():
    l, w, h = 1.0, 0.8, 1.2
    assert l * w * impl.vf_perpendicular_rectangles(l, w, h) == pytest.approx(l * h * impl.vf_perpendicular_rectangles(l, h, w), rel=1e-12)
    # very long common edge: the 2-D strip formula, which tends to 1/2 for a tall wall beside a narrow strip
    for hh in (0.7, 50.0):
        assert impl.vf_perpendicular_rectangles(1e5, 1.0, hh) == pytest.approx(impl.vf_perpendicular_strips(1.0, hh), rel=1e-3)
    assert impl.vf_perpendicular_strips(1.0, 1e4) == pytest.approx(0.5, rel=1e-4)


def test_coaxial_disks_reciprocity_and_limits():
    r1, r2, h = 0.1, 0.2, 0.3
    assert np.pi * r1 ** 2 * impl.vf_coaxial_disks(r1, r2, h) == pytest.approx(np.pi * r2 ** 2 * impl.vf_coaxial_disks(r2, r1, h), rel=1e-12)
    assert impl.vf_coaxial_disks(0.1, 100.0, 0.1) == pytest.approx(1.0, abs=1e-5)       # an infinite plane catches all
    # small disk far below a disk: F -> r2^2 / (h^2 + r2^2) (the disk seen from a point on the axis)
    assert impl.vf_coaxial_disks(1e-4, 0.5, 1.0) == pytest.approx(0.25 / 1.25, rel=1e-6)


def test_box_view_factors_obey_summation_and_reciprocity():
    for dims in ((1.0, 0.8, 0.5), (2.0, 0.3, 0.7), (1.0, 1.0, 1.0)):
        A, F = impl.box_view_factors(*dims)
        npt.assert_allclose(F.sum(axis=1), 1.0, rtol=1e-12)         # two independent closed forms close the box
        AF = A[:, None] * F
        npt.assert_allclose(AF, AF.T, rtol=1e-12, atol=1e-15)
        assert np.all(np.diag(F) == 0)
    A, F = impl.box_view_factors(1.0, 1.0, 1.0)
    assert F[0, 1] == pytest.approx(0.19982, abs=1e-4)              # cube: opposite faces 0.2, adjacent 0.2
    assert F[0, 2] == pytest.approx((1 - F[0, 1]) / 4, rel=1e-12)


# ---------------------------------------------------------------- independent check: numerical integration

def test_numerical_view_factor_matches_parallel_closed_form_at_second_order():
    a, b, c = 1.0, 0.8, 0.5
    Fx = impl.vf_parallel_rectangles(a, b, c)
    errs = [abs(impl.vf_numerical([0, 0, 0], [a, 0, 0], [0, b, 0], [0, 0, c], [a, 0, 0], [0, b, 0], n) - Fx)
            for n in (8, 16, 32)]
    orders = np.log2(np.array(errs[:-1]) / np.array(errs[1:]))
    assert np.all(orders > 1.8) and np.all(orders < 2.2)
    assert errs[-1] < 5e-4 * Fx


def test_numerical_view_factor_matches_perpendicular_closed_form_and_is_reciprocal():
    l, w, h = 1.0, 0.8, 1.2
    Fx = impl.vf_perpendicular_rectangles(l, w, h)
    errs = []
    for n in (8, 16, 32):
        F12 = impl.vf_numerical([0, 0, 0], [0, l, 0], [w, 0, 0], [0, 0, 0], [0, l, 0], [0, 0, h], n)
        F21 = impl.vf_numerical([0, 0, 0], [0, l, 0], [0, 0, h], [0, 0, 0], [0, l, 0], [w, 0, 0], n)
        assert l * w * F12 == pytest.approx(l * h * F21, rel=1e-12)  # the discrete kernel is symmetric
        errs.append(abs(F12 - Fx))
    orders = np.log2(np.array(errs[:-1]) / np.array(errs[1:]))
    assert np.all(orders > 0.8) and np.all(orders < 1.3)            # kernel singular along the shared edge
    assert errs[-1] < 0.03 * Fx


def test_crossed_strings_reproduce_two_dimensional_closed_forms():
    for w, h in ((1.0, 0.5), (0.4, 2.0)):
        assert impl.vf_crossed_strings((0, 0), (w, 0), (0, h), (w, h)) == pytest.approx(impl.vf_parallel_strips(w, h), rel=1e-12)
        assert impl.vf_crossed_strings((w, 0), (0, 0), (0, h), (0, 0)) == pytest.approx(impl.vf_perpendicular_strips(w, h), rel=1e-12)
    # a triangle (three long flat surfaces): F12 = (L1 + L2 - L3) / (2 L1), and each row sums to one
    P = [(0.0, 0.0), (3.0, 0.0), (0.5, 2.0)]
    L = [np.hypot(*np.subtract(P[1], P[0])), np.hypot(*np.subtract(P[2], P[1])), np.hypot(*np.subtract(P[0], P[2]))]
    F12 = impl.vf_crossed_strings(P[0], P[1], P[2], P[1])
    F13 = impl.vf_crossed_strings(P[1], P[0], P[2], P[0])
    assert F12 == pytest.approx((L[0] + L[1] - L[2]) / (2 * L[0]), rel=1e-12)
    assert F12 + F13 == pytest.approx(1.0, rel=1e-12)


# ---------------------------------------------------------------- enclosure

def _box_enclosure():
    A, F = impl.box_view_factors(1.0, 0.8, 0.5)
    return A, F


def test_enclosure_net_heat_sums_to_zero_and_reproduces_specified_values():
    A, F = _box_enclosure()
    eps = [0.9, 0.3, 0.6, 0.8, 0.5, 0.95]
    T = [500.0, None, 350.0, None, 300.0, 420.0]
    Q = [None, 150.0, None, 0.0, None, None]
    s = impl.solve_enclosure(A, F, eps, T=T, Q=Q)
    assert abs(s["Q"].sum()) < 1e-9 * np.max(np.abs(s["Q"]))
    assert s["Q"][1] == pytest.approx(150.0) and abs(s["Q"][3]) < 1e-9
    # re-solving with the computed temperatures everywhere gives the same heat flows
    s2 = impl.solve_enclosure(A, F, eps, T=list(s["T"]))
    npt.assert_allclose(s2["Q"], s["Q"], rtol=1e-9, atol=1e-8)


def test_black_enclosure_equals_direct_exchange_sum():
    A, F = _box_enclosure()
    T = np.array([500.0, 300.0, 350.0, 400.0, 320.0, 450.0])
    s = impl.solve_enclosure(A, F, np.ones(6), T=list(T))
    direct = np.array([sum(A[i] * F[i, j] * SIG * (T[i] ** 4 - T[j] ** 4) for j in range(6)) for i in range(6)])
    npt.assert_allclose(s["Q"], direct, rtol=1e-10)
    npt.assert_allclose(s["J"], SIG * T ** 4, rtol=1e-12)                     # black: radiosity = emissive power


def test_two_parallel_plates_and_isothermal_enclosure():
    e1, e2, T1, T2 = 0.8, 0.6, 600.0, 300.0
    s = impl.solve_enclosure([1.0, 1.0], [[0, 1], [1, 0]], [e1, e2], T=[T1, T2])
    plates = SIG * (T1 ** 4 - T2 ** 4) / (1 / e1 + 1 / e2 - 1)
    assert s["Q"][0] == pytest.approx(plates, rel=1e-12)
    assert impl.two_surface_exchange(1.0, 1.0, 1.0, e1, e2, T1, T2) == pytest.approx(plates, rel=1e-12)
    # concentric spheres r1 < r2: F12 = 1; the two-surface formula reduces to the textbook form
    r1, r2 = 0.3, 0.5
    A1, A2 = 4 * np.pi * r1 ** 2, 4 * np.pi * r2 ** 2
    s = impl.solve_enclosure([A1, A2], [[0, 1], [A1 / A2, 1 - A1 / A2]], [e1, e2], T=[T1, T2])
    assert s["Q"][0] == pytest.approx(A1 * SIG * (T1 ** 4 - T2 ** 4) / (1 / e1 + (A1 / A2) * (1 / e2 - 1)), rel=1e-12)
    # isothermal enclosure: no net exchange whatever the emittances (second law)
    A, F = _box_enclosure()
    s = impl.solve_enclosure(A, F, [0.1, 0.9, 0.5, 0.3, 0.7, 0.2], T=[400.0] * 6)
    assert np.max(np.abs(s["Q"])) < 1e-9 * SIG * 400.0 ** 4


def test_reradiating_surface_floats_between_hot_and_cold():
    Fp = impl.vf_parallel_rectangles(1.0, 1.0, 1.0)
    F = np.array([[0, Fp, 1 - Fp], [Fp, 0, 1 - Fp], [(1 - Fp) / 4, (1 - Fp) / 4, 0]])
    F[2, 2] = 1 - F[2, 0] - F[2, 1]
    s = impl.solve_enclosure([1, 1, 4], F, [0.8, 0.6, 0.5], T=[600.0, 300.0, None], Q=[None, None, 0.0])
    assert 300.0 < s["T"][2] < 600.0
    assert s["Q"][0] == pytest.approx(-s["Q"][1], rel=1e-12)
    # closed form of the network: the reradiating wall is a parallel path 1/(A1 F1R) + 1/(A2 F2R) beside 1/(A1 F12)
    A1 = A2 = 1.0
    R1, R2 = (1 - 0.8) / (0.8 * A1), (1 - 0.6) / (0.6 * A2)
    Rmid = 1 / (A1 * Fp + 1 / (1 / (A1 * (1 - Fp)) + 1 / (A2 * (1 - Fp))))
    assert s["Q"][0] == pytest.approx(SIG * (600.0 ** 4 - 300.0 ** 4) / (R1 + Rmid + R2), rel=1e-10)
    # reradiating wall: its emittance does not matter
    s2 = impl.solve_enclosure([1, 1, 4], F, [0.8, 0.6, 0.05], T=[600.0, 300.0, None], Q=[None, None, 0.0])
    assert s2["Q"][0] == pytest.approx(s["Q"][0], rel=1e-10)


def test_linearised_radiation_coefficients():
    eps, T, Ts = 0.9, 340.0, 295.0
    assert impl.h_rad_linear(eps, T, Ts) * (T - Ts) == pytest.approx(eps * SIG * (T ** 4 - Ts ** 4), rel=1e-12)
    rel_err = 1 - impl.h_rad_mean(eps, T, Ts) / impl.h_rad_linear(eps, T, Ts)
    assert rel_err == pytest.approx(((T - Ts) / (T + Ts)) ** 2 / (1 + ((T - Ts) / (T + Ts)) ** 2), rel=1e-10)
    assert 0 < rel_err < 0.01


# ---------------------------------------------------------------- radiator

def test_radiator_energy_balance_and_limits():
    Q, eps, A = 200.0, 0.93, 0.6
    T = impl.radiator_temperature(Q, A, eps, alpha_s=0.15, cos_sun=0.5)
    assert eps * SIG * A * (T ** 4 - impl.T_SPACE ** 4) == pytest.approx(Q + 0.15 * impl.G_SOLAR * 0.5 * A, rel=1e-12)
    A_ecl = impl.radiator_area(Q, 293.15, eps)
    assert A_ecl == pytest.approx(Q / (eps * SIG * (293.15 ** 4 - impl.T_SPACE ** 4)), rel=1e-12)
    assert impl.radiator_temperature(Q, A_ecl, eps) == pytest.approx(293.15, rel=1e-12)       # inverse of each other
    assert impl.radiator_area(Q, 293.15, eps, alpha_s=0.15, cos_sun=1.0) > A_ecl                # the sun costs area
    assert impl.radiator_temperature(Q, A, eps, alpha_s=0.15, cos_sun=-1.0) == pytest.approx(impl.radiator_temperature(Q, A, eps))
    with pytest.raises(ValueError):
        impl.radiator_area(Q, 293.15, 0.83, alpha_s=0.94, cos_sun=1.0)                         # black paint in full sun


def test_radiator_beside_a_solar_array():
    Q, eps = 200.0, 0.93
    alone = impl.radiator_temperature(Q, 0.6, eps)
    tiny = impl.radiator_with_array(Q, 1.0, 0.6, 1e-6, eps, 0.83, 330.0)
    assert tiny["T_rad"] == pytest.approx(alone, rel=1e-5)                    # no array: the isolated radiator
    cold = impl.radiator_with_array(Q, 1.0, 0.6, 1.0, eps, 1.0, 3.0)
    hot = impl.radiator_with_array(Q, 1.0, 0.6, 1.0, eps, 0.83, 330.0)
    assert hot["T_rad"] > alone                                               # a warm neighbour heats the radiator
    assert cold["T_rad"] == pytest.approx(alone, abs=0.05)                    # a black neighbour at 3 K is just more space
    assert hot["Q_array"] + hot["Q_space"] + Q == pytest.approx(0.0, abs=1e-6 * Q)


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "radiation.png").stat().st_size > 10_000
