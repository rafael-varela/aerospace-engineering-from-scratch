"""Tests for P09.17 Digital Control.

Run against the starter (default) or the reference solution (AEFS_IMPL=solution).
The independent checks: discretisation against ``scipy.signal.cont2discrete``; ZOH against the matrix
exponential (``scipy.linalg.expm``); step invariance of ZOH at the sample instants.
"""
from __future__ import annotations

import numpy as np
import numpy.testing as npt
import pytest
from conftest import load_impl
from scipy.linalg import expm
from scipy.signal import cont2discrete, dlti, dimpulse, lti, step

impl = load_impl(__file__)

RNG = np.random.default_rng(20250917)
SP = np.array([[-0.8, 1.0], [-4.0, -1.2]])
B2 = np.array([[0.0], [1.0]])
C2 = np.array([[1.0, 0.0]])
D2 = np.array([[0.0]])


def random_stable(n: int, m: int) -> tuple[np.ndarray, np.ndarray]:
    A = RNG.standard_normal((n, n))
    A = A - (np.max(np.linalg.eigvals(A).real) + 0.5) * np.eye(n)
    return A, RNG.standard_normal((n, m))


# ---------------------------------------------------------------- ZOH: the independent checks

@pytest.mark.parametrize("n,m,T", [(2, 1, 0.1), (3, 2, 0.05), (4, 1, 0.5), (5, 2, 0.01)])
def test_zoh_matches_cont2discrete_and_the_matrix_exponential(n, m, T):
    A, B = random_stable(n, m)
    Ad, Bd = impl.c2d_zoh(A, B, T)
    ref = cont2discrete((A, B, np.eye(n), np.zeros((n, m))), T, method="zoh")
    npt.assert_allclose(Ad, ref[0], rtol=1e-10, atol=1e-12)
    npt.assert_allclose(Bd, ref[1], rtol=1e-10, atol=1e-12)
    npt.assert_allclose(Ad, expm(A * T), rtol=1e-10, atol=1e-12)
    # Bd = A^{-1}(Ad - I) B for invertible A: a second closed form
    npt.assert_allclose(Bd, np.linalg.solve(A, (expm(A * T) - np.eye(n)) @ B), rtol=1e-9, atol=1e-12)


def test_zoh_of_an_integrator_and_the_small_T_limit():
    A = np.array([[0.0, 1.0], [0.0, 0.0]])                # double integrator
    B = np.array([[0.0], [1.0]])
    T = 0.3
    Ad, Bd = impl.c2d_zoh(A, B, T)
    npt.assert_allclose(Ad, [[1.0, T], [0.0, 1.0]], atol=1e-15)
    npt.assert_allclose(Bd, [[T ** 2 / 2], [T]], atol=1e-15)
    Ad, Bd = impl.c2d_zoh(SP, B2, 1e-6)
    npt.assert_allclose(Ad, np.eye(2) + 1e-6 * SP, atol=1e-11)
    npt.assert_allclose(Bd, 1e-6 * B2, rtol=1e-5, atol=1e-12)          # first row is T^2/2 = 5e-13


def test_step_invariance_at_the_sample_instants():
    """ZOH is exact for a step: the discrete step response equals the continuous one at t = kT."""
    A, B = impl.pitch_model()
    C = np.array([[0.0, 0.0, 1.0, 0.0]])
    T = 0.05
    Ad, Bd = impl.c2d_zoh(A - B @ impl.pitch_gains(), B * 2.0, T)   # closed loop driven by Kp * step
    y, _ = impl.dlsim(Ad, Bd, C, np.zeros((1, 1)), np.ones(81))
    t, yc = step(lti(A - B @ impl.pitch_gains(), B * 2.0, C, np.zeros((1, 1))), T=np.arange(81) * T)
    npt.assert_allclose(y[:, 0], yc, atol=1e-9)
    # Tustin is not step invariant: it differs at the samples
    At, Bt, Ct, Dt = impl.c2d_tustin(A - B @ impl.pitch_gains(), B * 2.0, C, np.zeros((1, 1)), T)
    yt, _ = impl.dlsim(At, Bt, Ct, Dt, np.ones(81))
    assert np.abs(yt[:, 0] - yc).max() > 1e-4


# ---------------------------------------------------------------- Tustin and matched pole-zero

@pytest.mark.parametrize("T", [0.01, 0.1, 0.4])
def test_tustin_state_space_matches_cont2discrete(T):
    A, B = random_stable(3, 1)
    C = RNG.standard_normal((2, 3))
    D = RNG.standard_normal((2, 1))
    got = impl.c2d_tustin(A, B, C, D, T)
    ref = cont2discrete((A, B, C, D), T, method="bilinear")
    for g, r in zip(got, ref[:4]):
        npt.assert_allclose(g, r, rtol=1e-10, atol=1e-12)


def test_tustin_tf_matches_cont2discrete_and_maps_the_axis_to_the_circle():
    num, den = np.array([0.5, 2.0]), np.array([0.05, 1.0])            # the lesson's lead
    for T in (0.02, 0.05, 0.2):
        nz, dz = impl.tustin_tf(num, den, T)
        rn, rd, _ = cont2discrete((num, den), T, method="bilinear")
        npt.assert_allclose(nz, np.ravel(rn) / rd[0], rtol=1e-10, atol=1e-12)
        npt.assert_allclose(dz, rd / rd[0], rtol=1e-10, atol=1e-12)
    # the s-plane frequency w_a maps to w_d = (2/T) atan(w_a T / 2): G_d(e^{j w_d T}) = G(j w_a)
    T = 0.05
    nz, dz = impl.tustin_tf(num, den, T)
    wa = np.array([0.5, 5.0, 30.0])
    wd = 2.0 / T * np.arctan(wa * T / 2)
    npt.assert_allclose(impl.freqresp_z(nz, dz, wd, T), impl.freqresp_s(num, den, wa), rtol=1e-10)


def test_prewarped_tustin_is_exact_at_the_prewarp_frequency():
    num, den = np.array([1.0]), np.array([1.0, 1.4, 1.0])            # second-order low pass, w_n = 1
    T, w0 = 0.5, 2.0
    nz, dz = impl.tustin_tf(num, den, T, w_prewarp=w0)
    npt.assert_allclose(impl.freqresp_z(nz, dz, [w0], T), impl.freqresp_s(num, den, [w0]), rtol=1e-12)
    nz0, dz0 = impl.tustin_tf(num, den, T)
    assert abs(impl.freqresp_z(nz0, dz0, [w0], T)[0] - impl.freqresp_s(num, den, [w0])[0]) > 1e-3


def test_matched_pole_zero_maps_poles_and_matches_dc():
    num, den = np.array([0.5, 2.0]), np.array([0.05, 1.0])
    T = 0.05
    nz, dz = impl.matched_tf(num, den, T)
    npt.assert_allclose(np.roots(dz), np.exp(np.roots(den) * T), rtol=1e-12)
    npt.assert_allclose(np.roots(nz), np.exp(np.roots(num) * T), rtol=1e-12)
    npt.assert_allclose(np.polyval(nz, 1.0) / np.polyval(dz, 1.0), np.polyval(num, 0.0) / np.polyval(den, 0.0), rtol=1e-12)
    # a strictly proper filter gets its zero at infinity mapped to z = -1
    nz2, dz2 = impl.matched_tf([4.0], [1.0, 3.0, 4.0], 0.1)
    npt.assert_allclose(np.sort(np.roots(nz2).real), [-1.0, -1.0], atol=1e-7)
    with pytest.raises(ValueError):
        impl.matched_tf([1.0], [1.0, 0.0], 0.1)                        # integrator: DC gain infinite


# ---------------------------------------------------------------- z-domain tools

def test_long_division_equals_the_impulse_response():
    num, den = np.array([0.5, 2.0]), np.array([0.05, 1.0])
    nz, dz = impl.tustin_tf(num, den, 0.05)
    h = impl.impulse_by_long_division(nz, dz, 30)
    _, (ref,) = dimpulse(dlti(nz, dz, dt=0.05), n=30)
    npt.assert_allclose(h, ref[:, 0], atol=1e-12)
    # first-order: 1/(z - a) = z^-1 + a z^-2 + ...
    npt.assert_allclose(impl.impulse_by_long_division([1.0], [1.0, -0.6], 6), [0, 1, 0.6, 0.36, 0.216, 0.1296], atol=1e-15)
    with pytest.raises(ValueError):
        impl.impulse_by_long_division([1.0, 0.0, 0.0], [1.0, 0.5], 4)


def test_dlsim_superposition_and_zero_input():
    Ad, Bd = impl.c2d_zoh(SP, B2, 0.1)
    u1, u2 = RNG.standard_normal(50), RNG.standard_normal(50)
    y1, _ = impl.dlsim(Ad, Bd, C2, D2, u1)
    y2, _ = impl.dlsim(Ad, Bd, C2, D2, u2)
    y12, _ = impl.dlsim(Ad, Bd, C2, D2, 2.0 * u1 - 3.0 * u2)
    npt.assert_allclose(y12, 2.0 * y1 - 3.0 * y2, atol=1e-12)
    y0, x0 = impl.dlsim(Ad, Bd, C2, D2, np.zeros(10))
    assert np.abs(y0).max() == 0.0
    # time invariance: delaying the input by 5 samples delays the output by 5 samples
    yd, _ = impl.dlsim(Ad, Bd, C2, D2, np.concatenate([np.zeros(5), u1[:-5]]))
    npt.assert_allclose(yd[5:], y1[:-5], atol=1e-12)


def test_stability_maps_from_the_left_half_plane_to_the_unit_disk():
    for T in (0.01, 0.1, 1.0):
        Ad, _ = impl.c2d_zoh(SP, B2, T)
        assert impl.spectral_radius(Ad) < 1.0
        npt.assert_allclose(np.sort_complex(impl.z_to_s(np.linalg.eigvals(Ad), T)), np.sort_complex(np.linalg.eigvals(SP)), atol=1e-9)
        At, *_ = impl.c2d_tustin(SP, B2, C2, D2, T)
        assert impl.spectral_radius(At) < 1.0                       # Tustin preserves stability at any T
    unstable = np.array([[0.2, 1.0], [-1.0, 0.2]])
    assert impl.spectral_radius(impl.c2d_zoh(unstable, B2, 0.1)[0]) > 1.0


# ---------------------------------------------------------------- computational delay

def test_delay_discretisation_limits():
    A, B = impl.pitch_model()
    T = 0.05
    Ad, Bd = impl.c2d_zoh(A, B, T)
    Phi, Gam = impl.c2d_zoh_delay(A, B, T, 0.0)
    npt.assert_allclose(Phi[:4, :4], Ad, atol=1e-14)
    npt.assert_allclose(Phi[:4, 4:], 0.0, atol=1e-14)
    npt.assert_allclose(Gam[:4], Bd, atol=1e-14)
    Phi, Gam = impl.c2d_zoh_delay(A, B, T, T)                       # a whole sample: u_{k-1} acts all period
    npt.assert_allclose(Phi[:4, 4:], Bd, atol=1e-14)
    npt.assert_allclose(Gam[:4], 0.0, atol=1e-14)
    # any tau: G0 + G1 = Bd (a constant input does not care when it switched)
    Phi, Gam = impl.c2d_zoh_delay(A, B, T, 0.3 * T)
    npt.assert_allclose(Phi[:4, 4:] + Gam[:4], Bd, atol=1e-13)
    with pytest.raises(ValueError):
        impl.c2d_zoh_delay(A, B, T, 1.5 * T)


@pytest.mark.parametrize("frac", [0.0, 0.25, 0.5, 1.0])
def test_delay_model_matches_a_simulation_resolved_between_samples(frac):
    T = 0.05
    th, xh = impl.hybrid_simulation(T, frac * T, 3.0, n_sub=20)
    Acl, Bcl = impl.discrete_pitch_loop(T, frac * T)
    _, zs = impl.dlsim(Acl, Bcl, np.eye(5), np.zeros((5, 1)), np.ones(60))
    npt.assert_allclose(xh[::20][:60], zs[:, :4], atol=1e-12)


# ---------------------------------------------------------------- the pitch loop

def test_fast_sampling_recovers_the_continuous_poles():
    A, B = impl.pitch_model()
    lam_c = np.sort_complex(np.linalg.eigvals(A - B @ impl.pitch_gains()))
    errs = []
    for T in (4e-3, 1e-3):
        Acl, _ = impl.discrete_pitch_loop(T, 0.0)
        lam = np.linalg.eigvals(Acl)
        s = np.sort_complex(impl.z_to_s(lam[np.abs(lam) > 1e-9], T))
        errs.append(np.max(np.abs(s - lam_c)))
    assert errs[1] < errs[0] / 3.0                                  # error shrinks roughly like T
    assert errs[1] < 0.15


def test_slower_sampling_and_longer_delay_cost_damping_and_margin():
    rows = impl.sample_rate_study(np.array([50.0, 20.0, 10.0]), delay_fracs=(0.0, 1.0))
    by = {(r["fs"], r["frac"]): r for r in rows}
    for fs in (50.0, 20.0, 10.0):
        assert by[(fs, 1.0)]["zeta"] < by[(fs, 0.0)]["zeta"]
        assert by[(fs, 1.0)]["pm"] < by[(fs, 0.0)]["pm"]
    assert by[(50.0, 1.0)]["pm"] > by[(20.0, 1.0)]["pm"] > by[(10.0, 1.0)]["pm"]
    assert all(r["rho"] < 1.0 for r in rows)
    # the rule PM_d ~ PM_c - wc (T/2 + tau) [rad -> deg] within 2 deg at 20 Hz and faster
    for fs in (50.0, 20.0):
        T = 1.0 / fs
        for frac in (0.0, 1.0):
            rule = 65.75 - np.degrees(7.12 * T * (0.5 + frac))
            assert abs(by[(fs, frac)]["pm"] - rule) < 2.0


def test_too_slow_a_flight_computer_destabilises_the_loop():
    Acl, _ = impl.discrete_pitch_loop(1.0 / 7.0, 1.0 / 7.0)
    assert impl.spectral_radius(Acl) > 1.0
    Acl, _ = impl.discrete_pitch_loop(1.0 / 50.0, 1.0 / 50.0)
    assert impl.spectral_radius(Acl) < 1.0


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "digital_pitch.png").stat().st_size > 1000
