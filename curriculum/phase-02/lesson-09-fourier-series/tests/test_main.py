"""Public tests for P02.09 Fourier Series. They assert properties, with scipy.integrate.quad, solve_ivp and SymPy as oracles.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-02/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.integrate import quad, solve_ivp

from conftest import load_impl

impl = load_impl(__file__)

CASE = impl.thruster_boom_case()
T = CASE["T"]
WN = 2.0 * np.pi * CASE["f_mode"]
ZETA = CASE["zeta"]


# ---------------------------------------------------------------- waveforms and closed-form coefficients

def test_waveforms_have_the_stated_symmetry_period_and_values():
    t = np.linspace(-2.0, 2.0, 4001)
    sq, tr = impl.square_wave(t, T), impl.triangle_wave(t, T)
    npt.assert_allclose(impl.square_wave(t + T, T), sq)                       # periodic
    npt.assert_allclose(impl.triangle_wave(t + T, T), tr)
    npt.assert_allclose(impl.square_wave(-t, T), -sq)                         # odd
    npt.assert_allclose(impl.triangle_wave(-t, T), tr)                        # even
    assert impl.square_wave(np.array([0.1 * T]), T)[0] == 1.0 and impl.square_wave(np.array([0.6 * T]), T)[0] == -1.0
    assert impl.triangle_wave(np.array([0.0]), T)[0] == pytest.approx(1.0) and impl.triangle_wave(np.array([0.5 * T]), T)[0] == pytest.approx(-1.0)
    assert impl.triangle_wave(np.array([0.25 * T]), T)[0] == pytest.approx(0.0)
    assert np.max(np.abs(sq)) == 1.0 and np.max(np.abs(tr)) == pytest.approx(1.0)
    npt.assert_allclose(impl.square_wave(t, T, 3.0), 3.0 * sq)


def test_closed_form_coefficients_match_quadrature_and_sympy():
    a0, an, bn = impl.square_wave_coefficients(8)
    assert a0 == 0.0 and np.all(an == 0.0)
    assert bn[0] == pytest.approx(4 / np.pi) and bn[1] == 0.0 and bn[4] == pytest.approx(4 / (5 * np.pi))
    for n in (1, 2, 3, 4, 7):
        q, _ = quad(lambda t: impl.square_wave(np.array([t]), T)[0] * np.sin(2 * np.pi * n * t / T), 0.0, T, points=[T / 2], limit=200)
        assert 2 * q / T == pytest.approx(bn[n - 1], abs=1e-10)
    import sympy as sp
    assert float(impl.sympy_square_wave_coefficient(3, T)) == pytest.approx(bn[2])
    assert impl.sympy_square_wave_coefficient(4, T) == 0
    a0t, ant, bnt = impl.triangle_wave_coefficients(8)
    assert a0t == 0.0 and np.all(bnt == 0.0) and ant[0] == pytest.approx(8 / np.pi ** 2) and ant[1] == 0.0
    for n in (1, 3, 5):
        q, _ = quad(lambda t: impl.triangle_wave(np.array([t]), T)[0] * np.cos(2 * np.pi * n * t / T), 0.0, T, points=[T / 2], limit=200)
        assert 2 * q / T == pytest.approx(ant[n - 1], abs=1e-10)


# ---------------------------------------------------------------- coefficients from scratch

def test_trapezoid_coefficients_match_closed_forms_with_first_order_error_at_a_jump():
    sq = lambda t: impl.square_wave(t, T)
    _, _, bn_exact = impl.square_wave_coefficients(10)
    errs = []
    for ns in (256, 512, 1024, 2048):
        a0, an, bn = impl.fourier_coefficients(sq, T, 10, n_samples=ns)
        assert abs(a0) < 1e-12 and np.max(np.abs(an)) < 1e-12                  # odd function: no cosines
        errs.append(np.max(np.abs(bn - bn_exact)))
    rates = np.log2(np.array(errs[:-1]) / np.array(errs[1:]))
    assert np.all(rates > 1.5), f"expected at least ~2nd-order error decay with midpoints at a jump, got {rates}"
    assert errs[-1] < 1e-5
    # A smooth periodic function: the trapezoid rule is exact to rounding once it resolves the harmonics.
    f = lambda t: 0.3 + 2.0 * np.cos(2 * np.pi * t / T) - 0.7 * np.sin(6 * np.pi * t / T)
    a0, an, bn = impl.fourier_coefficients(f, T, 5, n_samples=64)
    assert a0 == pytest.approx(0.6, abs=1e-13)
    npt.assert_allclose(an, [2.0, 0.0, 0.0, 0.0, 0.0], atol=1e-13)
    npt.assert_allclose(bn, [0.0, 0.0, -0.7, 0.0, 0.0], atol=1e-13)
    with pytest.raises(ValueError):
        impl.fourier_coefficients(f, T, 10, n_samples=20)                      # aliasing guard


def test_partial_sum_reproduces_a_finite_series_and_converges_in_mean_square():
    t = np.linspace(0.0, 2 * T, 801)
    f = lambda tt: 0.3 + 2.0 * np.cos(2 * np.pi * tt / T) - 0.7 * np.sin(6 * np.pi * tt / T)
    npt.assert_allclose(impl.partial_sum(t, 0.6, np.array([2.0, 0, 0]), np.array([0, 0, -0.7]), T), f(t), atol=1e-12)
    npt.assert_allclose(impl.partial_sum(t, 0.6, np.array([2.0, 0, 0]), np.array([0, 0, -0.7]), T, n_terms=1), 0.3 + 2 * np.cos(2 * np.pi * t / T), atol=1e-12)
    sq = lambda tt: impl.square_wave(tt, T)
    a0, an, bn = impl.square_wave_coefficients(201)
    e = [impl.truncation_error_rms(sq, T, a0, an, bn, N) for N in (5, 21, 81, 201)]
    assert e[0] > e[1] > e[2] > e[3]
    assert e[3] < 0.06
    # Parseval form of the same error: sqrt(1 - (1/2) sum_{n<=N} b_n^2).
    assert e[1] == pytest.approx(np.sqrt(1.0 - impl.parseval_power(a0, an[:21], bn[:21])), rel=0.02)
    with pytest.raises(ValueError):
        impl.partial_sum(t, a0, an, bn, T, n_terms=500)


def test_convergence_rate_follows_the_smoothness_of_the_waveform():
    Ns = np.array([5, 11, 21, 41, 81, 161])
    sq = lambda t: impl.square_wave(t, T)
    tr = lambda t: impl.triangle_wave(t, T)
    e_sq = np.array([impl.truncation_error_rms(sq, T, *impl.square_wave_coefficients(N), N) for N in Ns])
    e_tr = np.array([impl.truncation_error_rms(tr, T, *impl.triangle_wave_coefficients(N), N) for N in Ns])
    assert impl.convergence_rate(e_sq, Ns) == pytest.approx(0.5, abs=0.06)        # coefficients ~ 1/n: rms error ~ N^-1/2
    assert impl.convergence_rate(e_tr, Ns) == pytest.approx(1.5, abs=0.12)        # coefficients ~ 1/n^2: rms error ~ N^-3/2
    assert impl.convergence_rate(np.array([1.0, 0.25, 1 / 16]), np.array([1, 2, 4])) == pytest.approx(2.0)


def test_gibbs_overshoot_tends_to_the_gibbs_constant_and_does_not_shrink():
    g = np.array([impl.gibbs_overshoot(N, T) for N in (3, 9, 29, 99, 199)])
    assert np.all(np.diff(g) < 0.0)                                                # approaches from above
    assert g[-1] == pytest.approx(1.17898, abs=3e-4)                               # (2/pi) Si(pi): 8.95 % of the jump
    assert g[-1] - 1.0 > 0.17                                                      # the ear never goes away
    assert impl.gibbs_overshoot(1, T) == pytest.approx(4 / np.pi)                   # one term: the fundamental's peak


def test_parseval_identity_holds_for_each_waveform():
    sq = lambda t: impl.square_wave(t, T)
    tr = lambda t: impl.triangle_wave(t, T)
    assert impl.mean_square(sq, T) == pytest.approx(1.0)
    assert impl.mean_square(tr, T) == pytest.approx(1.0 / 3.0, rel=1e-6)
    assert impl.parseval_power(*impl.square_wave_coefficients(4001)) == pytest.approx(1.0, abs=2e-4)     # tail ~ 8/(pi^2 N) ... slow
    assert impl.parseval_power(*impl.triangle_wave_coefficients(201)) == pytest.approx(1.0 / 3.0, rel=1e-6)
    # Coefficients from the trapezoid rule satisfy Parseval against the sampled mean square too.
    f = lambda t: 1.0 + np.cos(2 * np.pi * t / T) ** 3
    a0, an, bn = impl.fourier_coefficients(f, T, 6, n_samples=256)
    assert impl.parseval_power(a0, an, bn) == pytest.approx(impl.mean_square(f, T), rel=1e-12)


def test_complex_form_round_trips_and_has_conjugate_symmetry():
    rng = np.random.default_rng(209)
    a0, an, bn = rng.normal(), rng.normal(size=6), rng.normal(size=6)
    c = impl.complex_coefficients(a0, an, bn)
    assert c.shape == (13,)
    npt.assert_allclose(c[::-1], np.conj(c))                                       # c_{-n} = conj(c_n)
    assert c[6] == pytest.approx(a0 / 2)
    a0r, anr, bnr = impl.from_complex_coefficients(c)
    assert a0r == pytest.approx(a0) and np.allclose(anr, an) and np.allclose(bnr, bn)
    # sum c_n e^{i n w0 t} over -N..N equals the real series.
    t = np.linspace(0.0, T, 101)
    n = np.arange(-6, 7)
    s_complex = (c[:, None] * np.exp(1j * np.outer(n, 2 * np.pi * t / T))).sum(axis=0)
    npt.assert_allclose(s_complex.imag, 0.0, atol=1e-12)
    npt.assert_allclose(s_complex.real, impl.partial_sum(t, a0, an, bn, T), atol=1e-12)
    with pytest.raises(ValueError):
        impl.from_complex_coefficients(np.zeros(4))


# ---------------------------------------------------------------- response to periodic forcing

def test_response_coefficients_reduce_to_the_static_and_resonant_limits():
    a0, an, bn = impl.square_wave_coefficients(15)
    A0, An, Bn = impl.response_coefficients(a0, an, bn, T, WN, ZETA)
    amp_f = impl.harmonic_amplitudes(an, bn)
    amp_q = impl.harmonic_amplitudes(An, Bn) * WN ** 2
    r = np.arange(1, 16) * 2 * np.pi / T / WN
    npt.assert_allclose(amp_q, amp_f * np.abs(impl.complex_gain(r, ZETA)), rtol=1e-12)
    assert impl.dominant_harmonic(An, Bn) == 5                                     # the mode sits on the 5th harmonic
    assert amp_q[4] == pytest.approx(amp_f[4] / (2 * ZETA), rel=1e-12)             # resonant gain 1 / (2 zeta)
    # A very stiff mode responds statically to every harmonic: q = f / wn^2.
    A0s, Ans, Bns = impl.response_coefficients(2.0, an, bn, T, 1e4, ZETA)
    assert A0s == pytest.approx(2.0 / 1e8)
    npt.assert_allclose(Ans * 1e8, an, atol=1e-3)                                # 2 zeta r ~ 1e-5 n leaks a little into A_n
    npt.assert_allclose(Bns * 1e8, bn, atol=1e-3)
    # Phase: at exact resonance the response of the 5th harmonic lags by 90 degrees (sin forcing -> -cos response).
    assert An[4] == pytest.approx(-bn[4] / (2 * ZETA) / WN ** 2, rel=1e-12)
    assert abs(Bn[4]) < 1e-12 * abs(An[4])


def test_harmonic_sum_matches_solve_ivp_in_steady_state():
    a0, an, bn = impl.square_wave_coefficients(51)
    A0, An, Bn = impl.response_coefficients(a0, an, bn, T, WN, ZETA)
    t_end = 40.0                                                                   # transient e^{-zeta wn t} = 3e-6
    t = np.linspace(t_end - 3.0, t_end, 3001)
    f = lambda tt: float(impl.square_wave(np.array([tt]), T)[0])
    sol = solve_ivp(impl.forced_mode_rhs, (0.0, t_end), [0.0, 0.0], t_eval=t, args=(f, WN, ZETA), rtol=1e-9, atol=1e-12, max_step=2e-3)
    q_series = impl.partial_sum(t, A0, An, Bn, T)
    scale = np.max(np.abs(q_series))
    # The solver sees the jumps only to within max_step and the series is truncated at n = 51: agree to 0.1 %.
    npt.assert_allclose(sol.y[0], q_series, atol=1e-3 * scale)
    # The square wave excites the mode about 14 times more than a unit sinusoid at the same fundamental would.
    assert scale * WN ** 2 > 13.0 and abs(impl.complex_gain(2 * np.pi / T / WN, ZETA)) < 1.1


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "fourier_series.png").stat().st_size > 1000
