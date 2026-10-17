"""Public tests for P09.13 State-Space Models.

They assert realisation properties, not magic numbers: every realisation (companion, observable,
modal, any similarity transform) has the same transfer function, checked pointwise as C(sI - A)^-1 B + D
against the polynomial TF; Faddeev-LeVerrier reproduces the characteristic polynomial and satisfies
Cayley-Hamilton; the physical pitch-altitude ODE model has exactly the hand-derived TF (the lesson's
independent check); simulation by the matrix exponential agrees with ``scipy.signal.lti`` (step and
``lsim`` with a held input); the ZOH matrices equal ``scipy.linalg.expm`` and a quadrature; superposition,
time invariance and the DC gain hold.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy import signal
from scipy.linalg import expm

from conftest import load_impl

impl = load_impl(__file__)

RNG = np.random.default_rng(91301)
S_TEST = np.array([0.3 + 0.7j, 2.0j, -0.4 + 3.0j, 5.0 + 0.0j, 0.1 - 9.0j])

TFS = [
    (np.array([2.0, 1.0]), np.array([1.0, 3.0, 2.0])),                     # strictly proper, 2nd order
    (np.array([1.0, 0.5, 4.0]), np.array([2.0, 1.0, 3.0, 1.0])),           # non-monic den, 3rd order
    (np.array([3.0, 1.0, 0.0, 2.0]), np.array([1.0, 2.0, 5.0, 4.0])),      # biproper: D = 3
]


def random_stable_system(n):
    A = RNG.standard_normal((n, n))
    A = A - (np.max(np.linalg.eigvals(A).real) + 0.5) * np.eye(n)     # shift so every pole has Re < -0.5
    return A, RNG.standard_normal((n, 1)), RNG.standard_normal((1, n)), np.array([[0.3]])


# ---------------------------------------------------------------- TF <-> SS

@pytest.mark.parametrize("num,den", TFS)
def test_companion_and_observable_forms_reproduce_the_tf_pointwise(num, den):
    for form in (impl.tf_to_controllable, impl.tf_to_observable):
        A, B, C, D = form(num, den)
        assert A.shape == (len(den) - 1, len(den) - 1)
        npt.assert_allclose(impl.ss_eval(A, B, C, D, S_TEST), impl.tf_eval(num, den, S_TEST), rtol=1e-12)


def test_companion_form_is_the_scipy_tf2ss_form():
    num, den = TFS[1]
    A, B, C, D = impl.tf_to_controllable(num, den)
    ref = signal.lti(num, den).to_ss()
    npt.assert_allclose(A, ref.A, atol=1e-14)
    npt.assert_allclose(B, ref.B, atol=1e-14)
    npt.assert_allclose(C, ref.C, atol=1e-14)
    npt.assert_allclose(D, ref.D, atol=1e-14)
    with pytest.raises(ValueError):
        impl.tf_to_controllable(np.array([1.0, 0.0, 0.0]), np.array([1.0, 1.0]))   # improper


@pytest.mark.parametrize("n", [3, 5])
def test_faddeev_leverrier_gives_the_characteristic_polynomial_and_cayley_hamilton(n):
    A = RNG.standard_normal((n, n))
    c, Ns = impl.faddeev_leverrier(A)
    npt.assert_allclose(c, np.poly(A), rtol=1e-10, atol=1e-10)
    # Cayley-Hamilton residual: A N_n + c_n I = 0
    npt.assert_allclose(A @ Ns[-1] + c[-1] * np.eye(n), 0.0, atol=1e-10)
    # adj(sI - A) = det(sI - A) (sI - A)^{-1} at a test point
    s = 0.7 + 1.3j
    adj = sum(N * s ** (n - k) for k, N in enumerate(Ns, start=1))
    npt.assert_allclose(adj, np.polyval(c, s) * np.linalg.inv(s * np.eye(n) - A), rtol=1e-9, atol=1e-10)


def test_ss_to_tf_matches_scipy_and_round_trips():
    A, B, C, D = random_stable_system(4)
    num, den = impl.ss_to_tf(A, B, C, D)
    ref = signal.lti(A, B, C, D).to_tf()
    npt.assert_allclose(den, ref.den, rtol=1e-9, atol=1e-10)
    npt.assert_allclose(num, ref.num, rtol=1e-9, atol=1e-10)
    npt.assert_allclose(impl.tf_eval(num, den, S_TEST), impl.ss_eval(A, B, C, D, S_TEST), rtol=1e-10)
    # round trip TF -> companion -> TF
    n2, d2 = impl.ss_to_tf(*impl.tf_to_controllable(*TFS[2]))
    npt.assert_allclose(n2, TFS[2][0] / TFS[2][1][0], atol=1e-12)
    npt.assert_allclose(d2, TFS[2][1] / TFS[2][1][0], atol=1e-12)


# ---------------------------------------------------------------- coordinates

def test_similarity_transform_changes_nothing_observable():
    A, B, C, D = random_stable_system(4)
    T = RNG.standard_normal((4, 4)) + 2.0 * np.eye(4)
    At, Bt, Ct, Dt = impl.similarity(A, B, C, D, T)
    npt.assert_allclose(np.sort_complex(np.linalg.eigvals(At)), np.sort_complex(np.linalg.eigvals(A)), atol=1e-9)
    npt.assert_allclose(impl.ss_eval(At, Bt, Ct, Dt, S_TEST), impl.ss_eval(A, B, C, D, S_TEST), rtol=1e-9)
    npt.assert_allclose(At, np.linalg.solve(T, A @ T), atol=1e-10)


def test_modal_form_is_real_block_diagonal_with_the_same_tf():
    A, B, C, D = random_stable_system(5)
    Am, Bm, Cm, Dm, T = impl.modal_form(A, B, C, D)
    assert np.isrealobj(Am) and np.isrealobj(T)
    # zero outside the 1x1 / 2x2 diagonal blocks
    mask = np.ones_like(Am, dtype=bool)
    i = 0
    while i < 5:
        if i + 1 < 5 and abs(Am[i, i + 1]) > 1e-8:
            mask[i:i + 2, i:i + 2] = False
            assert Am[i, i] == pytest.approx(Am[i + 1, i + 1], abs=1e-9)          # sigma on both
            assert Am[i, i + 1] == pytest.approx(-Am[i + 1, i], abs=1e-9)         # +-omega
            i += 2
        else:
            mask[i, i] = False
            i += 1
    npt.assert_allclose(Am[mask], 0.0, atol=1e-9)
    npt.assert_allclose(impl.ss_eval(Am, Bm, Cm, Dm, S_TEST), impl.ss_eval(A, B, C, D, S_TEST), rtol=1e-9)


# ---------------------------------------------------------------- the aircraft: ODE model vs hand TF

def test_pitch_altitude_ode_model_equals_the_hand_transfer_function():
    p = impl.pitch_altitude_params()
    A, B, C, D = impl.pitch_altitude_ss(p)
    num, den = impl.pitch_altitude_tf_hand(p)
    npt.assert_allclose(impl.ss_eval(A, B, C, D, S_TEST), impl.tf_eval(num, den, S_TEST), rtol=1e-12)
    nf, df = impl.ss_to_tf(A, B, C, D)
    npt.assert_allclose(impl.tf_eval(nf, df, S_TEST), impl.tf_eval(num, den, S_TEST), rtol=1e-10)
    # two poles at the origin (theta and h are integrals), two at the short period
    lam = np.sort_complex(np.linalg.eigvals(A))
    npt.assert_allclose(np.abs(lam[np.abs(lam) < 1e-6]), 0.0, atol=1e-8)
    assert np.sum(np.abs(lam) < 1e-6) == 2


def test_pitch_altitude_tf_scales_with_airspeed():
    p = impl.pitch_altitude_params()
    q = dict(p, V=2 * p["V"])
    n1, d1 = impl.pitch_altitude_tf_hand(p)
    n2, d2 = impl.pitch_altitude_tf_hand(q)
    s = 1.0 + 1.0j
    assert impl.tf_eval(n2, d2, s) == pytest.approx(2.0 * impl.tf_eval(n1, d1, s))


# ---------------------------------------------------------------- simulation against scipy.signal.lti

def test_step_simulation_matches_scipy_lti():
    p = impl.pitch_altitude_params()
    A, B, C, D = impl.pitch_altitude_ss(p)
    dt = 0.01
    t = np.arange(0.0, 5.0 + dt / 2, dt)
    y, _ = impl.simulate(A, B, C, D, np.ones_like(t), dt)
    _, y_ref = signal.step(signal.lti(A, B, C, D), T=t)
    npt.assert_allclose(y, y_ref, rtol=1e-9, atol=1e-9 * np.abs(y_ref).max())
    # the companion realisation of the hand TF gives the same output
    yc, _ = impl.simulate(*impl.tf_to_controllable(*impl.pitch_altitude_tf_hand(p)), np.ones_like(t), dt)
    npt.assert_allclose(yc, y, rtol=1e-9, atol=1e-9 * np.abs(y).max())


def test_held_random_input_matches_scipy_lsim_without_interpolation():
    A, B, C, D = random_stable_system(4)
    dt = 0.05
    t = np.arange(200) * dt
    u = RNG.standard_normal(200)
    x0 = RNG.standard_normal(4)
    y, x = impl.simulate(A, B, C, D, u, dt, x0)
    _, y_ref, x_ref = signal.lsim(signal.lti(A, B, C, D), u, t, X0=x0, interp=False)
    npt.assert_allclose(y, y_ref, atol=1e-9)
    npt.assert_allclose(x, x_ref, atol=1e-9)


def test_discretize_matches_expm_and_quadrature():
    A, B, _, _ = random_stable_system(3)
    dt = 0.2
    Ad, Bd = impl.discretize_zoh(A, B, dt)
    npt.assert_allclose(Ad, expm(A * dt), rtol=1e-10, atol=1e-12)
    s = np.linspace(0.0, dt, 4001)
    Bq = np.trapezoid(np.array([expm(A * si) @ B for si in s]), s, axis=0)
    npt.assert_allclose(Bd, Bq, rtol=1e-6, atol=1e-9)


def test_superposition_time_invariance_and_free_response():
    A, B, C, D = random_stable_system(3)
    dt = 0.05
    u1, u2 = RNG.standard_normal(100), RNG.standard_normal(100)
    y1, _ = impl.simulate(A, B, C, D, u1, dt)
    y2, _ = impl.simulate(A, B, C, D, u2, dt)
    y12, _ = impl.simulate(A, B, C, D, 2.0 * u1 - 3.0 * u2, dt)
    npt.assert_allclose(y12, 2.0 * y1 - 3.0 * y2, atol=1e-12)
    # shifting the input by 10 samples shifts the output by 10 samples
    ys, _ = impl.simulate(A, B, C, D, np.concatenate([np.zeros(10), u1[:-10]]), dt)
    npt.assert_allclose(ys[10:], y1[:-10], atol=1e-12)
    npt.assert_allclose(ys[:10], 0.0, atol=0.0)
    # zero input: x_k = e^{A k dt} x0
    x0 = np.array([1.0, -0.5, 0.2])
    _, x = impl.simulate(A, B, C, D, np.zeros(40), dt, x0)
    npt.assert_allclose(x[39], expm(A * 39 * dt) @ x0, rtol=1e-9, atol=1e-12)


def test_dc_gain_equals_tf_at_zero_and_the_step_final_value():
    A, B, C, D = random_stable_system(3)
    num, den = impl.ss_to_tf(A, B, C, D)
    g0 = impl.dc_gain(A, B, C, D)
    assert g0 == pytest.approx(impl.tf_eval(num, den, 0.0).real, rel=1e-10)
    dt = 0.05
    y, _ = impl.simulate(A, B, C, D, np.ones(2000), dt)                 # 100 s: slowest pole Re < -0.5 -> e^-50
    assert y[-1] == pytest.approx(g0, rel=1e-8, abs=1e-10)


def test_short_period_modal_block_and_observable_duality():
    p = impl.pitch_altitude_params()
    A, B, _, _ = impl.pitch_altitude_ss(p)
    Asp, Bsp, Csp, Dsp = A[:2, :2], B[:2], np.array([[1.0, 0.0]]), np.zeros((1, 1))
    Am, _, _, _, _ = impl.modal_form(Asp, Bsp, Csp, Dsp)
    za, ma, mq = p["Z_alpha"], p["M_alpha"], p["M_q"]
    sigma = (za + mq) / 2
    omega = np.sqrt((za * mq - ma) - sigma ** 2)
    npt.assert_allclose(np.diag(Am), [sigma, sigma], atol=1e-12)
    assert abs(Am[0, 1]) == pytest.approx(omega, rel=1e-12) and Am[0, 1] == pytest.approx(-Am[1, 0], rel=1e-12)
    Ac, Bc, Cc, Dc = impl.tf_to_controllable(*TFS[1])
    Ao, Bo, Co, Do = impl.tf_to_observable(*TFS[1])
    npt.assert_allclose(Ao, Ac.T)
    npt.assert_allclose(Bo, Cc.T)
    npt.assert_allclose(Co, Bc.T)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "state_space.png").stat().st_size > 1000
