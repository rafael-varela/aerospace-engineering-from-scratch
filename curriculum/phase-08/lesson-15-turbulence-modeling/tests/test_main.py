"""Public tests for P08.15 Turbulence Modeling.

They assert properties of the closures, not one magic number: the total stress is the exact linear
profile of the momentum balance; the velocity gradient solves the mixing-length quadratic and reduces
to the laminar value when the mixing length vanishes; the kappa -> 0 limit is plane Poiseuille flow
with c_f = 12/Re_b; the model has a viscous sublayer u+ = y+ and a log layer whose slope is 1/kappa;
the bulk velocity converges at second order with the grid; the friction law brackets Dean's
correlation within the stated band; the results depend only on Re_b (unit invariance); the k-epsilon
log layer has production = dissipation, nu_t = kappa u_tau y and -<u'v'>/k = sqrt(C_mu), and its
epsilon equation closes only for kappa^2 = sigma_eps sqrt(C_mu)(C_eps2 - C_eps1).

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)


# ---------------------------------------------------------------- mixing-length building blocks

def test_total_stress_is_the_linear_momentum_balance():
    yp = np.linspace(0, 550, 12)
    tau = impl.total_stress(yp, 550.0)
    assert tau[0] == 1.0 and abs(tau[-1]) < 1e-15
    npt.assert_allclose(np.diff(tau, 2), 0.0, atol=1e-14)            # straight line


def test_velocity_gradient_solves_the_quadratic_and_has_both_limits():
    rng = np.random.default_rng(15)
    tau = rng.uniform(0.0, 1.0, 200)
    lp = rng.uniform(0.0, 400.0, 200)
    g = impl.velocity_gradient(tau, lp)
    npt.assert_allclose(lp ** 2 * g ** 2 + g, tau, rtol=1e-12, atol=1e-15)  # (1 + l+^2 g) g = tau+
    assert np.all(g >= 0)
    npt.assert_allclose(impl.velocity_gradient(tau, np.zeros_like(tau)), tau, rtol=1e-15)   # laminar
    big = impl.velocity_gradient(np.array([0.64]), np.array([1e6]))
    assert big[0] == pytest.approx(np.sqrt(0.64) / 1e6, rel=1e-5)     # fully turbulent: g = sqrt(tau)/l


def test_van_driest_and_nikuradse_mixing_length_limits():
    yp = np.array([0.5, 1.0, 2.0])
    lp = impl.mixing_length(yp, 1000.0)
    npt.assert_allclose(lp, impl.KAPPA * yp ** 2 / impl.A_PLUS, rtol=0.05)   # damped: kappa y+^2/A+ near wall
    far = impl.mixing_length(np.array([500.0]), 1000.0)
    assert far[0] == pytest.approx(impl.KAPPA * 500.0, rel=1e-6)
    ln = impl.mixing_length(np.array([1e-3, 1000.0]), 1000.0, A_plus=0.0, outer="nikuradse")
    assert ln[0] == pytest.approx(0.4e-3, rel=1e-3) and ln[1] == pytest.approx(0.14 * 1000.0)


def test_kappa_to_zero_is_plane_poiseuille_flow_with_cf_12_over_re():
    prof = impl.mixing_length_profile(800.0, kappa=1e-12)
    npt.assert_allclose(prof["u_plus"], impl.laminar_profile(prof["y_plus"], 800.0), atol=1e-8)
    re_b, cf = impl.friction(prof)
    assert cf * re_b == pytest.approx(12.0, rel=1e-5)


# ---------------------------------------------------------------- the turbulent profile

def test_viscous_sublayer_and_log_layer_structure():
    prof = impl.mixing_length_profile(5200.0)
    yp, up, g, tau = prof["y_plus"], prof["u_plus"], prof["dudy"], prof["tau"]
    m = (yp > 0) & (yp < 3)
    npt.assert_allclose(up[m], yp[m], rtol=0.01)                       # u+ = y+ in the sublayer
    m = (yp > 200) & (yp < 1000)
    npt.assert_allclose(yp[m] * g[m] * impl.KAPPA / np.sqrt(tau[m]), 1.0, atol=0.01)   # slope 1/kappa
    # the eddy viscosity is kappa y+ sqrt(tau+) there, far larger than the molecular one
    npt.assert_allclose(prof["nu_t"][m] / (impl.KAPPA * yp[m] * np.sqrt(tau[m])), 1.0, atol=0.015)
    assert np.all(prof["nu_t"][m] > 50)
    # turbulent + viscous stress add to the total stress everywhere
    npt.assert_allclose(prof["reynolds_stress"] + g, tau, atol=1e-12)


def test_log_law_fit_recovers_kappa_and_damping_sets_B():
    Bs = []
    for A in (20.0, 26.0, 30.0):
        prof = impl.mixing_length_profile(2000.0, A_plus=A)
        k, B = impl.fit_log_law(prof["y_plus"], prof["u_plus"])
        assert k == pytest.approx(impl.KAPPA, abs=0.015)
        Bs.append(B)
    assert Bs[0] < Bs[1] < Bs[2]                       # a thicker damped layer shifts the log law up
    assert 4.5 < Bs[1] < 5.8                           # A+ = 26 gives B near the typical 5.0-5.5


def test_bulk_velocity_converges_at_second_order():
    ub = [impl.bulk_velocity(impl.mixing_length_profile(1000.0, n=n)) for n in (251, 501, 1001, 2001)]
    d = np.abs(np.diff(ub))
    orders = np.log2(d[:-1] / d[1:])
    assert np.all((orders > 1.8) & (orders < 2.2)), orders


def test_friction_brackets_the_dean_correlation():
    # Dean's c_f = 0.073 Re_b^-1/4 is a CORRELATION; the model is a two-constant closure. The band
    # (within 20 %) is the honest claim; the outer-region closure alone moves c_f by about 15 %.
    cfs = []
    for rt in (180.0, 550.0, 2000.0, 5200.0):
        re_b, cf = impl.friction(impl.mixing_length_profile(rt))
        assert abs(cf / impl.dean_cf(re_b) - 1) < 0.20
        cfs.append(cf)
    assert np.all(np.diff(cfs) < 0)
    re_b, cf_n = impl.friction(impl.mixing_length_profile(1000.0, outer="nikuradse"))
    re_b2, cf_p = impl.friction(impl.mixing_length_profile(1000.0))
    assert cf_n < impl.dean_cf(re_b) < cf_p                          # the two outer closures bracket it
    for rt in (500.0, 1000.0, 2000.0):                               # log-law friction law vs Dean
        re_b, cf = impl.log_law_friction(rt)
        assert abs(cf / impl.dean_cf(re_b) - 1) < 0.06


def test_re_tau_inversion_and_unit_invariance():
    rt = impl.re_tau_for_re_b(30000.0)
    assert impl.friction(impl.mixing_length_profile(rt))[0] == pytest.approx(30000.0, rel=1e-9)
    air = impl.duct_flow(5.0, 0.05)                                    # air
    water = impl.duct_flow(5.0 * 0.05 / 0.02 * (1.002e-3 / 998.2) / (1.81e-5 / 1.204), 0.02,
                           rho=998.2, mu=1.002e-3)                      # water, same Re_b
    assert water["re_b"] == pytest.approx(air["re_b"], rel=1e-12)
    assert water["cf"] == pytest.approx(air["cf"], rel=1e-9)
    assert water["re_tau"] == pytest.approx(air["re_tau"], rel=1e-9)
    # momentum balance on the duct: dp/dx times the half-height balances the wall shear
    assert air["dpdx"] * 0.025 == pytest.approx(-air["tau_w"], rel=1e-12)
    assert air["y1"] * air["u_tau"] / (1.81e-5 / 1.204) == pytest.approx(1.0)


# ---------------------------------------------------------------- k-epsilon in the log layer

def test_keps_log_layer_equilibrium_relations():
    y = np.linspace(0.01, 0.2, 50)
    u_tau = 0.3
    p = impl.keps_log_layer(y, u_tau)
    P = p["nu_t"] * p["dudy"] ** 2
    npt.assert_allclose(P, p["eps"], rtol=1e-12)                       # production = dissipation
    npt.assert_allclose(p["nu_t"], impl.KAPPA * u_tau * y, rtol=1e-12)  # same nu_t as the mixing length
    uv = p["nu_t"] * p["dudy"]
    npt.assert_allclose(uv, u_tau ** 2, rtol=1e-12)                     # constant stress
    npt.assert_allclose(uv / p["k"], np.sqrt(impl.C_MU), rtol=1e-12)     # -<u'v'>/k = sqrt(C_mu)
    assert impl.c_mu_from_stress_ratio(0.3) == pytest.approx(0.09)


def test_keps_constants_imply_kappa_0433():
    k = impl.kappa_from_keps()
    assert k == pytest.approx(np.sqrt(1.3 * 0.3 * 0.48), rel=1e-12)
    assert 0.43 < k < 0.435                                             # not the 0.41 of the log law


def test_epsilon_equation_closes_only_with_the_consistent_kappa():
    kc = impl.kappa_from_keps()
    res = []
    for n in (101, 201, 401, 801):
        y = np.linspace(0.01, 0.1, n)
        rk, re = impl.keps_residuals(y, 0.5, kappa=kc)
        assert np.max(np.abs(rk)) < 1e-12                               # k is constant: P = eps exactly
        res.append(np.max(np.abs(re)))
    orders = np.log2(np.array(res[:-1]) / np.array(res[1:]))
    assert np.all(orders > 1.8)                                         # residual is truncation, O(h^2)
    _, re = impl.keps_residuals(np.geomspace(0.01, 0.1, 401), 0.5, kappa=0.41)
    expected = (0.41 ** 2 / (1.3 * 0.3) - 0.48) / 1.92                  # analytic mismatch, about -2.6 %
    npt.assert_allclose(re, expected, atol=1e-4)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "turbulence_models.png").stat().st_size > 10_000
