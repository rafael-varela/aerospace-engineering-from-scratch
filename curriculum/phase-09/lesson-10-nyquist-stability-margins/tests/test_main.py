"""Public tests for P09.10 Nyquist & Stability Margins. They assert properties, not magic numbers: the
D-contour is closed, avoids the imaginary-axis poles and turns clockwise; the winding number counts
known curves correctly; Z = N + P equals the number of closed-loop right-half-plane roots from the
characteristic polynomial for stable and unstable gains, for a loop with a double pole at the origin,
for an open-loop unstable airframe (where stability needs N = -1), and for loops with a pure delay
(checked against a Pade approximant; independent check 1); the gain and phase margins equal the
crossings read from ``scipy.signal.bode`` (independent check 2); the gain margin equals the Routh
critical gain divided by the gain; the delay margin puts the curve through -1; and the limits K -> 0
and K at the critical gain behave as derived.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy import signal

from conftest import load_impl

impl = load_impl(__file__)

K_CRIT_SERVO = (impl.TAU_M + impl.TAU_F) / (impl.TAU_M * impl.TAU_F)     # Routh: 120 1/s


# ---------------------------------------------------------------- contour and winding number

def test_contour_is_closed_clockwise_and_avoids_the_imaginary_axis_poles():
    num, den = impl.satellite_loop()                                    # double pole at s = 0
    s = impl.nyquist_contour(den, R=1e3, eps=1e-3)
    assert np.min(np.abs(s)) >= 1e-3 * (1 - 1e-9)                         # never on the pole
    assert abs(s[0] - (-1e3j)) < 1e-6 and abs(s[-1] - (-1e3j)) < 1e-6    # starts and ends at -jR
    # the contour turns clockwise about a point inside the right half plane
    assert impl.winding_number(s, 5.0 + 1.0j) == pytest.approx(-1.0, abs=1e-9)
    assert impl.winding_number(s, -5.0 + 1.0j) == pytest.approx(0.0, abs=1e-9)
    assert impl.winding_number(s, 0.0 + 0.0j) == pytest.approx(0.0, abs=1e-9)  # the indented pole is outside


def test_winding_number_of_known_curves():
    th = np.linspace(0, 2 * np.pi, 400, endpoint=False)
    circle = 0.5 * np.exp(1j * th) - 1.0                                # around -1, counterclockwise
    assert impl.winding_number(circle) == pytest.approx(1.0)
    assert impl.winding_number(circle[::-1]) == pytest.approx(-1.0)
    assert impl.winding_number(-1.0 + 0.5 * np.exp(2j * th)) == pytest.approx(2.0)
    assert impl.winding_number(3.0 + 0.5 * np.exp(1j * th)) == pytest.approx(0.0, abs=1e-12)


def test_imag_axis_poles_are_found():
    npt.assert_allclose(impl.imag_axis_poles(np.convolve([1.0, 0.0], [1.0, 0.0, 4.0])), [-2.0, 0.0, 2.0], atol=1e-9)
    assert impl.imag_axis_poles(np.array([1.0, 3.0, 2.0])).size == 0


# ---------------------------------------------------------------- Z = N + P against the closed-loop poles

@pytest.mark.parametrize("K", [5.0, 20.0, 100.0, 119.0, 121.0, 300.0])
def test_servo_nyquist_count_agrees_with_closed_loop_poles(K):
    num, den = impl.servo_loop(K)
    N, P, Z = impl.nyquist_count(num, den)
    assert P == 0
    assert Z == int(np.sum(np.roots(impl.poly_add(den, num)).real > 0))
    assert Z == (0 if K < K_CRIT_SERVO else 2)


@pytest.mark.parametrize("K", [3.0, 6.0, 7.0, 20.0, 100.0])
def test_open_loop_unstable_airframe_needs_one_counterclockwise_encirclement(K):
    num, den = impl.unstable_airframe_loop(K)
    N, P, Z = impl.nyquist_count(num, den)
    assert P == 1
    assert Z == int(np.sum(np.roots(impl.poly_add(den, num)).real > 0))
    assert (N == -1) == (K > 20.0 / 3.0)                                # stable only above K = 20/3


@pytest.mark.parametrize("frac", [0.0, 0.5, 0.9, 1.1, 2.0])
def test_delayed_satellite_loop_nyquist_agrees_with_pade_closed_loop_poles(frac):
    num, den = impl.satellite_loop()
    dm = impl.margins(num, den)["DM"]
    T = frac * dm
    N, P, Z = impl.nyquist_count(num, den, delay=T)
    assert Z == impl.closed_loop_rhp_roots(num, den, T, pade_order=8)
    assert Z == (0 if frac < 1 else 2)


def test_pade_approximant_matches_the_delay_at_low_frequency():
    pn, pd = impl.pade_delay(0.7, 6)
    w = np.linspace(0.0, 1.5, 40)                                       # omega T up to about 1
    npt.assert_allclose(impl.tf_eval(pn, pd, 1j * w), np.exp(-0.7j * w), atol=1e-9)
    npt.assert_allclose(np.abs(impl.tf_eval(pn, pd, 1j * np.logspace(-2, 3, 50))), 1.0, rtol=1e-9)   # all-pass


# ---------------------------------------------------------------- margins against Bode and Routh

def bode_crossings(num, den):
    w = np.geomspace(1e-3, 1e3, 200001)
    _, mag, ph = signal.bode((num, den), w=w)
    i = np.flatnonzero(np.sign(mag[:-1]) != np.sign(mag[1:]))[0]
    w_gc = np.interp(0.0, [mag[i + 1], mag[i]], [w[i + 1], w[i]])
    pm = 180.0 + np.interp(w_gc, w[i:i + 2], ph[i:i + 2])
    pm = (pm + 180.0) % 360.0 - 180.0
    j = np.flatnonzero(np.diff(np.floor((ph + 180.0) / 360.0)) != 0)
    out = {"w_gc": w_gc, "PM": pm}
    if j.size:
        k = j[0]
        out["w_pc"] = np.interp(0.0, [(ph[k] + 180) % 360 - 360 if (ph[k] + 180) % 360 > 180 else (ph[k] + 180) % 360,
                                      (ph[k + 1] + 180) % 360 - 360 if (ph[k + 1] + 180) % 360 > 180 else (ph[k + 1] + 180) % 360],
                                [w[k], w[k + 1]])
        out["GM_db"] = -np.interp(out["w_pc"], w[k:k + 2], mag[k:k + 2])
    return out


@pytest.mark.parametrize("K", [5.0, 20.0, 60.0])
def test_servo_margins_agree_with_the_bode_crossings(K):
    num, den = impl.servo_loop(K)
    m = impl.margins(num, den)
    b = bode_crossings(num, den)
    assert m["w_gc"] == pytest.approx(b["w_gc"], rel=1e-4)
    assert m["PM"] == pytest.approx(b["PM"], abs=0.01)
    assert m["w_pc"] == pytest.approx(b["w_pc"], rel=1e-4)
    assert m["GM_db"] == pytest.approx(b["GM_db"], abs=0.01)


def test_servo_gain_margin_equals_routh_critical_gain_over_k():
    for K in (5.0, 20.0, 60.0):
        gm, wpc = impl.gain_margin(*impl.servo_loop(K))
        assert gm == pytest.approx(K_CRIT_SERVO / K, rel=1e-8)
        assert wpc == pytest.approx(1.0 / np.sqrt(impl.TAU_M * impl.TAU_F), rel=1e-8)


def test_satellite_margins_agree_with_bode_and_gain_margin_is_infinite():
    num, den = impl.satellite_loop()
    m = impl.margins(num, den)
    assert np.isinf(m["GM"])                                            # the phase never reaches -180 for w > 0
    b = bode_crossings(num, den)
    assert m["PM"] == pytest.approx(b["PM"], abs=0.01)
    assert m["w_gc"] == pytest.approx(b["w_gc"], rel=1e-4)
    # analytic: arg L = -180 + atan(Td w) - atan(tau w)
    w = m["w_gc"]
    assert m["PM"] == pytest.approx(np.degrees(np.arctan(impl.TD_SAT * w) - np.arctan(impl.TAU_SAT * w)), abs=1e-6)


def test_delay_margin_puts_the_curve_through_minus_one():
    num, den = impl.satellite_loop()
    m = impl.margins(num, den)
    assert m["DM"] == pytest.approx(np.radians(m["PM"]) / m["w_gc"])
    s = 1j * m["w_gc"]
    assert abs(1.0 + impl.tf_eval(num, den, s) * np.exp(-s * m["DM"])) < 1e-9
    pm_half, _ = impl.phase_margin(num, den, delay=0.5 * m["DM"])
    assert pm_half == pytest.approx(0.5 * m["PM"], rel=1e-6)           # a delay costs w_gc T of phase, not gain


def test_limits_small_gain_and_critical_gain():
    gm_small, _ = impl.gain_margin(*impl.servo_loop(0.5))
    assert gm_small == pytest.approx(K_CRIT_SERVO / 0.5, rel=1e-6)      # small K: large margin, GM ~ 1/K
    pm_small, w = impl.phase_margin(*impl.servo_loop(0.5))
    assert w == pytest.approx(0.5, rel=1e-3)                            # crossover where K/w = 1, far below both lags
    assert pm_small == pytest.approx(90.0 - np.degrees(np.arctan(impl.TAU_M * w) + np.arctan(impl.TAU_F * w)), abs=1e-6)
    assert pm_small > 88.0                                              # K -> 0: PM -> 90 deg (the integrator alone)
    num, den = impl.servo_loop(K_CRIT_SERVO)
    m = impl.margins(num, den)
    assert m["GM"] == pytest.approx(1.0, rel=1e-8) and m["PM"] == pytest.approx(0.0, abs=1e-6)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "nyquist_margins.png").exists()
