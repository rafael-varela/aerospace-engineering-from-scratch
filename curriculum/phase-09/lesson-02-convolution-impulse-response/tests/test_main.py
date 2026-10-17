"""Public tests for P09.02 Convolution & Impulse Response.

The independent checks: the direct convolution sum against FFT multiplication with zero padding (done
here with ``numpy.fft``, not with your code) and against ``numpy.convolve``; continuous convolution by
Riemann sum against closed forms, with the convergence order measured; the first-order step response by
convolution against the closed form and against ``scipy.signal`` (oracle only). Also: the algebra of
convolution (commutative, associative, distributive, the impulse as identity, delay by a shifted
impulse), the wrap-around of unpadded circular convolution, and superposition of shifted steps.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
from __future__ import annotations

import numpy as np
import numpy.testing as npt
import pytest
from scipy import signal

from conftest import load_impl

impl = load_impl(__file__)
RNG = np.random.default_rng(20260902)


def fft_product_oracle(x: np.ndarray, h: np.ndarray) -> np.ndarray:
    """Linear convolution by FFT multiplication with zero padding to N + M - 1 (written here, independent of the impl)."""
    L = len(x) + len(h) - 1
    return np.fft.irfft(np.fft.rfft(x, L) * np.fft.rfft(h, L), L)


# ---------------------------------------------------------------- discrete convolution

@pytest.mark.parametrize("N,M", [(1, 1), (7, 3), (64, 64), (500, 37)])
def test_direct_sum_matches_fft_multiplication_with_zero_padding(N, M):
    x, h = RNG.standard_normal(N), RNG.standard_normal(M)
    y = impl.convolve_direct(x, h)
    assert y.shape == (N + M - 1,)
    npt.assert_allclose(y, fft_product_oracle(x, h), atol=1e-11)
    npt.assert_allclose(y, np.convolve(x, h), atol=1e-12)


@pytest.mark.parametrize("N,M", [(10, 4), (300, 300), (1000, 77)])
def test_learner_fft_convolution_equals_the_direct_sum(N, M):
    x, h = RNG.standard_normal(N), RNG.standard_normal(M)
    y = impl.convolve_fft(x, h)
    assert y.shape == (N + M - 1,) and np.isrealobj(y)
    npt.assert_allclose(y, impl.convolve_direct(x, h), atol=1e-10)


def test_circular_convolution_wraps_around_unless_padded():
    x, h = RNG.standard_normal(50), RNG.standard_normal(12)
    lin = impl.convolve_direct(x, h)
    # unpadded: the circular result equals the DFT product of length L = 50
    c50 = impl.circular_convolve(x, h, 50)
    npt.assert_allclose(c50, np.fft.ifft(np.fft.fft(x) * np.fft.fft(h, 50)).real, atol=1e-12)
    # the last M - 1 linear samples have wrapped onto the first M - 1; beyond them it is exact
    npt.assert_allclose(c50[:11], lin[:11] + lin[50:], atol=1e-12)
    npt.assert_allclose(c50[11:], lin[11:50], atol=1e-12)
    assert np.max(np.abs(c50[:11] - lin[:11])) > 1e-3
    # padded to N + M - 1 or more, circular == linear
    npt.assert_allclose(impl.circular_convolve(x, h, 61), lin, atol=1e-12)
    npt.assert_allclose(impl.circular_convolve(x, h, 64)[:61], lin, atol=1e-12)
    npt.assert_allclose(impl.circular_convolve(x, h, 64)[61:], 0.0, atol=1e-12)


def test_convolution_algebra():
    a, b, c = RNG.standard_normal(20), RNG.standard_normal(9), RNG.standard_normal(14)
    cd = impl.convolve_direct
    npt.assert_allclose(cd(a, b), cd(b, a), atol=1e-12)                                  # commutative
    npt.assert_allclose(cd(cd(a, b), c), cd(a, cd(b, c)), atol=1e-11)                    # associative
    b14 = np.concatenate([b, np.zeros(5)])
    npt.assert_allclose(cd(a, b14 + c), cd(a, b14) + cd(a, c), atol=1e-12)               # distributive
    npt.assert_allclose(cd(a, [1.0]), a, atol=0)                                          # impulse is the identity
    npt.assert_allclose(cd(a, [0.0, 0.0, 0.0, 1.0]), np.concatenate([np.zeros(3), a]))  # shifted impulse delays
    assert np.sum(cd(a, b)) == pytest.approx(np.sum(a) * np.sum(b), rel=1e-12)           # areas multiply


def test_causal_lti_output_matches_the_difference_equation():
    x = RNG.standard_normal(400)
    a, b = 0.95, 0.05
    h = b * a ** np.arange(400)
    y_ref = signal.lfilter([b], [1.0, -a], x)                      # oracle: the recursion itself
    npt.assert_allclose(impl.causal_lti_output(h, x), y_ref, atol=1e-12)


# ---------------------------------------------------------------- continuous convolution

@pytest.mark.parametrize("a,b", [(1.0, 3.0), (2.0, 4.0), (0.5, 5.0)])
def test_riemann_sum_converges_to_the_closed_form_at_first_and_second_order(a, b):
    errs = {"rect": [], "trap": []}
    for dt in (0.02, 0.01, 0.005):
        t = np.arange(0.0, 10.0, dt)
        exact = impl.exp_conv_exp(t, a, b)
        for rule in errs:
            y = impl.convolve_continuous(np.exp(-a * t), np.exp(-b * t), dt, rule)
            errs[rule].append(np.max(np.abs(y - exact)))
    order = lambda e: np.log2(np.array(e[:-1]) / np.array(e[1:]))
    npt.assert_allclose(order(errs["rect"]), 1.0, atol=0.1)
    npt.assert_allclose(order(errs["trap"]), 2.0, atol=0.1)
    assert errs["trap"][-1] < 1e-4


def test_closed_form_exp_conv_exp_is_symmetric_and_has_the_right_limit():
    t = np.linspace(0.0, 6.0, 601)
    npt.assert_allclose(impl.exp_conv_exp(t, 1.0, 3.0), impl.exp_conv_exp(t, 3.0, 1.0), atol=1e-15)
    npt.assert_allclose(impl.exp_conv_exp(t, 2.0, 2.0 + 1e-7), impl.exp_conv_exp(t, 2.0, 2.0), atol=1e-6)
    npt.assert_array_equal(impl.exp_conv_exp(np.array([-1.0, -0.1]), 1.0, 2.0), 0.0)
    # a == b: the integrand e^{-a tau} e^{-a (t - tau)} = e^{-a t} is constant in tau, so the trapezoid rule is exact
    dt = 0.05
    tt = np.arange(0.0, 5.0, dt)
    npt.assert_allclose(impl.convolve_continuous(np.exp(-2 * tt), np.exp(-2 * tt), dt, "trap"), impl.exp_conv_exp(tt, 2.0, 2.0), atol=1e-13)


def test_rect_convolved_with_itself_is_a_triangle():
    w, dt = 0.5, 1e-3
    t = np.arange(0.0, 2.0, dt)
    r = ((t >= 0) & (t < w)).astype(float)
    y = impl.convolve_continuous(r, r, dt, "rect")
    tri = impl.rect_conv_rect(t, w)
    npt.assert_allclose(y, tri, atol=2 * dt)
    assert tri.max() == pytest.approx(w) and t[np.argmax(tri)] == pytest.approx(w, abs=dt)
    assert np.sum(tri) * dt == pytest.approx(w * w, rel=1e-3)                      # area = product of areas


# ---------------------------------------------------------------- first-order systems

@pytest.mark.parametrize("K,tau", [(1.0, 1.6e-3), (2.5, 0.05)])
def test_first_order_step_by_convolution_matches_closed_form_and_scipy(K, tau):
    dt = tau / 50
    t = np.arange(0.0, 8 * tau, dt)
    h = impl.first_order_impulse_response(t, K, tau)
    s_conv = impl.convolve_continuous(h, np.ones_like(t), dt, "trap")
    s_cf = impl.first_order_step_response(t, K, tau)
    npt.assert_allclose(s_conv, s_cf, atol=1e-3 * K)                               # O(dt^2) ~ K (dt/tau)^2 / 12
    _, s_sp = signal.step(signal.lti([K], [tau, 1.0]), T=t)
    npt.assert_allclose(s_cf, s_sp, atol=1e-9 * K)
    # the area of h is the DC gain
    assert np.sum(h) * dt == pytest.approx(K, rel=0.02)


def test_step_and_impulse_responses_are_integral_and_derivative_of_each_other():
    tau, dt = 0.05, 1e-4
    t = np.arange(0.0, 0.5, dt)
    h = impl.first_order_impulse_response(t, 1.0, tau)
    s = impl.first_order_step_response(t, 1.0, tau)
    npt.assert_allclose(impl.step_from_impulse(h, dt), s, atol=1e-6)
    npt.assert_allclose(impl.impulse_from_step(s, dt)[1:-1], h[1:-1], rtol=1e-5)


def test_mean_delay_of_the_impulse_response_is_tau():
    for tau in (1.6e-3, 0.05):
        dt = tau / 200
        h = impl.first_order_impulse_response(np.arange(0.0, 30 * tau, dt), 1.0, tau)
        assert impl.impulse_response_centroid(h, dt) == pytest.approx(tau, abs=dt)  # sampled centroid is tau - dt/2


def test_sampled_first_order_system_has_unit_dc_gain_and_matches_sampling():
    K, tau, dt = 1.7, 0.05, 1e-3
    h = impl.sampled_first_order(K, tau, dt, 2000)
    assert np.sum(h) == pytest.approx(K, rel=1e-9)                                 # unit DC gain times K
    a = np.exp(-dt / tau)
    npt.assert_allclose(h[1:] / h[:-1], a, rtol=1e-12)                             # pole at e^{-dt/tau}
    # for dt << tau it approaches dt times the continuous impulse response
    npt.assert_allclose(h[:100], dt * impl.first_order_impulse_response(np.arange(100) * dt, K, tau), rtol=0.011)


def test_doublet_response_by_convolution_equals_superposed_shifted_steps():
    K, tau, dt = 1.0, 0.05, 2.5e-4
    t = np.arange(0.0, 2.0, dt)
    cmd = impl.doublet(t, 2.0, 0.5, start=0.1)
    y = impl.convolve_continuous(cmd, impl.first_order_impulse_response(t, K, tau), dt, "trap")
    yc = impl.doublet_response_closed_form(t, K, tau, 2.0, 0.5, start=0.1)
    # the command jumps, so even the trapezoid rule is first order there: error ~ 2 A dt / tau
    assert np.max(np.abs(y - yc)) < 2 * 2.0 * dt / tau + 1e-12
    _, y_sp, _ = signal.lsim(signal.lti([K], [tau, 1.0]), cmd, t)
    assert np.max(np.abs(yc - y_sp)) < 2 * 2.0 * dt / tau


# ---------------------------------------------------------------- figure

def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "convolution.png").stat().st_size > 1000
