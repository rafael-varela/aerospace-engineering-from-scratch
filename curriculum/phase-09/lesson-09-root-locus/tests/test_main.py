"""Public tests for P09.09 Root Locus. They assert properties, not magic numbers: every computed locus
point satisfies the angle condition arg G = 180 deg and the magnitude condition K |G| = 1 (independent
check 1), and points off the locus do not; branches are continuous, start at the open-loop poles and,
for n - m >= 2, keep a constant root sum; far out they follow the asymptotes; the real-axis rule agrees
with the angle condition; breakaway points are double roots; the Routh count equals the count from
np.roots, including the epsilon and row-of-zeros cases; the j omega crossing from the Im/Re conditions
agrees with Routh (the RHP count jumps from 0 to 2 there, and the auxiliary polynomial gives the same
omega: independent check 2); and the gain for a damping target lands a root on the zeta line, equal
to the closed form where one exists.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy import signal

from conftest import load_impl

impl = load_impl(__file__)

NUM, DEN = impl.altitude_plant()                                     # 50 / (0.5 s^3 + 1.5 s^2 + s)
PLANTS = {
    "altitude P": (NUM, DEN),
    "altitude PD Td=2": (np.convolve(NUM, [2.0, 1.0]), DEN),
    "complex poles": (np.array([1.0, 1.0]), np.convolve([1.0, 0.0], [1.0, 2.0, 5.0])),
}


# ---------------------------------------------------------------- the locus and its two conditions

@pytest.mark.parametrize("name", list(PLANTS))
def test_every_locus_point_satisfies_the_angle_and_magnitude_conditions(name):
    num, den = PLANTS[name]
    gains = np.logspace(-4, 2, 300)
    R = impl.locus_roots(num, den, gains)
    assert R.shape == (gains.size, len(den) - 1)
    for K, roots in zip(gains, R):
        npt.assert_allclose(impl.angle_condition_error(num, den, roots), 0.0, atol=1e-6)
        npt.assert_allclose(K * np.abs(impl.tf_eval(num, den, roots)), 1.0, rtol=1e-6)
        npt.assert_allclose(impl.gain_at_point(num, den, roots), K, rtol=1e-6)


def test_points_off_the_locus_fail_the_angle_condition():
    for s in (-0.5 + 2.0j, 0.5 + 0.5j, -1.5 + 0.0j, -3.0 + 1.0j):
        assert abs(impl.angle_condition_error(NUM, DEN, s)) > 1.0


def test_locus_starts_at_the_open_loop_poles_and_branches_are_continuous():
    gains = np.concatenate([[0.0], np.logspace(-6, 0, 800)])
    R = impl.locus_roots(NUM, DEN, gains)
    npt.assert_allclose(np.sort_complex(R[0]), np.sort_complex(np.roots(DEN).astype(complex)), atol=1e-12)
    for K, roots in zip(gains[::50], R[::50]):                        # the same set as np.roots
        npt.assert_allclose(np.sort_complex(roots), np.sort_complex(np.roots(impl.char_poly(NUM, DEN, K))), atol=1e-9)
    assert np.max(np.abs(np.diff(R, axis=0))) < 0.15                   # no branch jumps


def test_root_sum_is_constant_when_the_pole_excess_is_at_least_two():
    R = impl.locus_roots(NUM, DEN, np.logspace(-4, 1, 100))
    npt.assert_allclose(np.sum(R, axis=1), -DEN[1] / DEN[0], atol=1e-9)


def test_far_branches_follow_the_asymptotes():
    sigma, angles = impl.asymptotes(NUM, DEN)
    assert sigma == pytest.approx(-1.0)
    npt.assert_allclose(np.sort(angles), [60.0, 180.0, 300.0])
    roots = np.roots(impl.char_poly(NUM, DEN, 1e6))
    dirs = np.sort(np.degrees(np.angle(roots - sigma)) % 360.0)
    npt.assert_allclose(dirs, [60.0, 180.0, 300.0], atol=0.5)
    s2, a2 = impl.asymptotes(*PLANTS["altitude PD Td=2"])               # one zero at -0.5: two asymptotes
    assert s2 == pytest.approx((-3.0 + 0.5) / 2) and np.allclose(np.sort(a2), [90.0, 270.0])


def test_real_axis_segments_agree_with_the_angle_condition():
    for name, (num, den) in PLANTS.items():
        segs = impl.real_axis_segments(num, den)
        for x in np.linspace(-6.0, 1.0, 141):
            if np.min(np.abs(np.concatenate([np.roots(den).real, np.roots(num).real if len(num) > 1 else []]) - x)) < 1e-6:
                continue
            on = any(lo < x < hi for lo, hi in segs)
            assert on == (abs(impl.angle_condition_error(num, den, x + 0j)) < 1e-6), (name, x)


def test_breakaway_point_is_a_double_root_and_the_off_locus_candidate_is_rejected():
    b = impl.breakaway_points(NUM, DEN)
    npt.assert_allclose(b, [-1.0 + 1.0 / np.sqrt(3.0)], atol=1e-9)    # -1.577 also solves dK/ds = 0 but needs K < 0
    K = -np.polyval(DEN, b[0]) / np.polyval(NUM, b[0])
    roots = np.roots(impl.char_poly(NUM, DEN, K))
    assert np.sum(np.abs(roots - b[0]) < 1e-4) == 2


# ---------------------------------------------------------------- Routh and the j omega crossing

@pytest.mark.parametrize("seed", range(6))
def test_routh_count_equals_the_count_from_np_roots(seed):
    rng = np.random.default_rng(seed)
    roots = rng.normal(0.0, 2.0, 3) + 0j
    roots = np.concatenate([roots, [complex(rng.normal(), rng.uniform(0.5, 3))]])
    roots = np.concatenate([roots, [np.conj(roots[-1])]])
    c = np.real(np.poly(roots)) * rng.uniform(0.5, 3.0)
    assert impl.routh_rhp_count(c) == int(np.sum(roots.real > 0))


def test_routh_special_cases():
    assert impl.routh_rhp_count([1.0, 1.0, 2.0, 2.0, 3.0]) == 2          # zero in the first column: epsilon
    assert impl.routh_rhp_count([1.0, 1.0, -4.0, -4.0]) == 1             # (s^2 - 4)(s + 1): row of zeros
    assert impl.routh_rhp_count([1.0, 1.0, 1.0, 1.0]) == 0               # (s^2 + 1)(s + 1): roots on the axis
    R = impl.routh_array([1.0, 1.0, -4.0, -4.0])
    npt.assert_allclose(R[2, 0], 2.0)                                   # d/ds of the auxiliary s^2 - 4


def test_jw_crossing_agrees_with_routh():
    (Kc, wc), = impl.jw_crossings(NUM, DEN)
    a3, a2, a1, _ = DEN
    assert Kc == pytest.approx(a2 * a1 / (a3 * NUM[0]), rel=1e-10)       # Routh's s^1 row vanishes
    assert impl.routh_rhp_count(impl.char_poly(NUM, DEN, 0.999 * Kc)) == 0
    assert impl.routh_rhp_count(impl.char_poly(NUM, DEN, 1.001 * Kc)) == 2
    R = impl.routh_array(impl.char_poly(NUM, DEN, Kc))
    w_aux = np.sqrt(R[1, 1] / R[1, 0])                                  # auxiliary a2 s^2 + K V = 0
    assert wc == pytest.approx(w_aux, rel=1e-8) and wc == pytest.approx(np.sqrt(2.0), rel=1e-10)
    assert np.min(np.abs(np.roots(impl.char_poly(NUM, DEN, Kc)) - 1j * wc)) < 1e-6


def test_rate_feedback_removes_the_jw_crossing():
    num, den = PLANTS["altitude PD Td=2"]
    assert impl.jw_crossings(num, den) == []
    for K in (1e-3, 1e-1, 10.0, 1e3):
        assert impl.routh_rhp_count(impl.char_poly(num, den, K)) == 0


# ---------------------------------------------------------------- damping target

@pytest.mark.parametrize("zeta", [0.3, 0.5, 0.7])
def test_gain_for_damping_puts_a_closed_loop_root_on_the_zeta_line(zeta):
    K, s = impl.gain_for_damping(NUM, DEN, zeta)
    assert -s.real / abs(s) == pytest.approx(zeta, rel=1e-6)
    assert np.min(np.abs(np.roots(impl.char_poly(NUM, DEN, K)) - s)) < 1e-6
    assert 0 < K < impl.jw_crossings(NUM, DEN)[0][0]                    # a stable design


def test_gain_for_damping_matches_the_second_order_closed_form():
    # Td = 1 s puts the zero on the pole at -1: L = 50 K / (s (0.5 s + 1)), s^2 + 2 s + 100 K = 0
    num = np.convolve(NUM, [1.0, 1.0])
    K, s = impl.gain_for_damping(num, DEN, 0.5)
    assert K == pytest.approx(0.04, rel=1e-8)                          # wn = 1/zeta = 2, wn^2 = 100 K
    assert abs(s) == pytest.approx(2.0, rel=1e-8)


def test_given_step_response_agrees_with_scipy_and_the_dominant_pair_where_it_should():
    K = 0.04
    num = np.convolve(NUM, [1.0, 1.0])
    cl = impl.char_poly(num, DEN, K)
    t = np.arange(0.0, 20.0, 0.01)
    y = impl.step_response(K * num, cl, t)
    _, ys = signal.step((K * num, cl), T=t)
    npt.assert_allclose(y, ys, atol=1e-6)
    os_, _ = impl.overshoot_and_settling(y, t)
    assert os_ == pytest.approx(100 * np.exp(-np.pi * 0.5 / np.sqrt(0.75)), abs=0.1)   # exact second order: 16.3 %


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "root_locus.png").exists()
