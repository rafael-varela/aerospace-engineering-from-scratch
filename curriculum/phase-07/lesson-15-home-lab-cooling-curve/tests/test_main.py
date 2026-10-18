"""Public tests for P07.15 Home Lab: Cooling Curve. They assert properties: the bar's geometry and
lumped time constant follow their definitions; the convection correlation vanishes at zero excess
temperature and grows with it; the Biot number of a small aluminium bar is far below 0.1; first-order
uncertainty propagation agrees with a seeded Monte Carlo; the log-linear, nonlinear, lagged and
variable-h fits recover planted parameters from noise-free data and give honest uncertainties
(chi^2/dof near one) on seeded noise; a probe's lag does not bias tau once its transient has died,
but a wrong ambient temperature does; the series eigenvalues solve their equations and tend to the
lumped limits; the series solutions conserve energy (d mean/dFo = -2 Bi theta_surface) and reduce to
the lumped model at small Bi; the nonlinear lumped model closes its energy balance; and the SYNTHETIC
files are labelled, seeded, quantized, physically plausible and read back.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.special import j0, j1

from conftest import load_impl

impl = load_impl(__file__)

D, L = impl.D_BAR, impl.L_BAR
TINF = 294.15


def _bar():
    g = impl.bar_geometry(D, L)
    return g, impl.RHO_AL * impl.C_AL * g["V"]


# ---------------------------------------------------------------- the object and the prediction

def test_geometry_and_time_constant_definitions():
    g = impl.bar_geometry(D, L)
    assert g["V"] == pytest.approx(np.pi * D ** 2 / 4 * L)
    assert g["A"] == pytest.approx(np.pi * D * L + np.pi * D ** 2 / 2)
    assert g["Lc"] == pytest.approx(g["V"] / g["A"])
    assert impl.bar_geometry(D, 1e4)["Lc"] == pytest.approx(D / 4, rel=1e-5)        # long cylinder: R/2
    assert impl.time_constant(2700.0, 896.0, g["V"], 7.0, g["A"]) == pytest.approx(2700 * 896 * g["V"] / (7 * g["A"]))
    assert impl.biot_number(10.0, 167.0, 0.005) == pytest.approx(10 * 0.005 / 167)


def test_air_props_reproduce_table_nodes():
    for row in impl.AIR_A6:
        p = impl.air_props(row[0])
        assert p["nu"] == pytest.approx(row[1] * 1e-5) and p["k"] == pytest.approx(row[2]) and p["Pr"] == pytest.approx(row[4])
    with pytest.raises(ValueError):
        impl.air_props(400.0)


def test_convection_correlation_behaviour():
    assert impl.h_conv_bar(TINF, TINF, D, L)["h"] == 0.0
    hs = [impl.h_conv_bar(TINF + th, TINF, D, L)["h"] for th in (2.0, 8.0, 16.0, 32.0)]
    assert np.all(np.diff(hs) > 0)
    # laminar natural convection: (Nu - const) ~ Ra^(1/4) ~ dT^(1/4): from 8 K to 32 K (x4) h rises by
    # less than 4^(1/4) = 1.41 (the constant term and the film-temperature shift reduce it)
    assert 1.15 < hs[3] / hs[1] < 1.41
    r = impl.h_conv_bar(TINF + 14.0, TINF, D, L)
    assert r["h_side"] < r["h_end"] and r["h_side"] < r["h"] < r["h_end"]
    assert 3.0 < r["h"] < 10.0                                     # still air: a few W/m^2 K
    assert impl.h_rad_linear(0.1, 310.0, TINF) * 15.85 == pytest.approx(0.1 * impl.SIGMA * (310.0 ** 4 - TINF ** 4), rel=1e-12)


def test_convection_correlation_coefficients():
    # the correlations written out here, with film properties interpolated from Table A.6 (nu, k, alpha
    # are columns 1-3 in 1e-5 m^2/s, W/m K, 1e-5 m^2/s; Pr column 4), pin every coefficient
    Ts = TINF + 14.0
    Tf = 0.5 * (Ts + TINF)
    a6 = np.asarray(impl.AIR_A6)
    nu, k, alpha, Pr = (float(np.interp(Tf, a6[:, 0], a6[:, j])) for j in (1, 2, 3, 4))
    Ra = 9.80665 * (1 / Tf) * 14.0 * D ** 3 / (nu * 1e-5 * alpha * 1e-5)
    Nu_side = 0.36 + 0.518 * Ra ** 0.25 / (1 + (0.559 / Pr) ** (9 / 16)) ** (4 / 9)    # Churchill-Chu cylinder
    Nu_end = 0.68 + 0.670 * Ra ** 0.25 / (1 + (0.492 / Pr) ** (9 / 16)) ** (4 / 9)     # Churchill-Chu plate
    A_side, A_end = np.pi * D * L, np.pi * D ** 2 / 2
    r = impl.h_conv_bar(Ts, TINF, D, L)
    assert r["h_side"] == pytest.approx(Nu_side * k / D, rel=1e-3)
    assert r["h_end"] == pytest.approx(Nu_end * k / D, rel=1e-3)
    assert r["h"] == pytest.approx((Nu_side * A_side + Nu_end * A_end) * k / D / (A_side + A_end), rel=1e-3)


def test_prediction_is_lumped_and_responds_to_emittance():
    p = impl.predict_tau(14.0, TINF)
    g, C = _bar()
    assert p["tau"] == pytest.approx(C / (p["h"] * g["A"]))
    assert p["h"] == pytest.approx(p["h_conv"] + p["h_rad"])
    assert p["Bi"] < 1e-3                                           # far below 0.1: lumped is justified
    assert impl.predict_tau(14.0, TINF, eps=0.9)["tau"] < p["tau"] < impl.predict_tau(14.0, TINF, eps=0.0)["tau"]
    # a steel bar of the same size has a larger rho c, so it cools more slowly
    assert impl.predict_tau(14.0, TINF, rho=impl.RHO_1010, c=impl.C_1010)["tau"] > p["tau"]


def test_gum_and_monte_carlo_agree():
    u = {"D": 0.05e-3, "L": 0.1e-3, "eps": 0.04, "rho": 13.5, "c": 9.0, "Tinf": 0.5, "f_conv": 0.10}
    r = impl.predict_tau_uncertainty(14.0, TINF, u, n=4000, seed=7)
    # 4000 draws: the sample sd has a relative standard error of about 1/sqrt(2 n) = 1.1 %; the model is
    # mildly nonlinear in f_conv, so allow 10 % between first-order and Monte Carlo
    assert r["mc_sd"] == pytest.approx(r["u_gum"], rel=0.10)
    assert r["mc_mean"] == pytest.approx(r["tau"], rel=0.03)
    assert max(r["terms"], key=lambda k: abs(r["terms"][k])) == "f_conv"   # the correlation dominates
    assert np.sqrt(sum(v ** 2 for v in r["terms"].values())) == pytest.approx(r["u_gum"])


# ---------------------------------------------------------------- fits

def _exp_data(tau=1800.0, theta0=28.0, Tinf=TINF, dt=15.0, t_end=3600.0):
    t = np.arange(0.0, t_end + 1e-9, dt)
    return t, Tinf + theta0 * np.exp(-t / tau)


def test_fits_recover_noise_free_exponentials():
    t, T = _exp_data()
    ll = impl.fit_loglinear(t, T, TINF, 0.05, theta_min=0.5)
    assert ll["tau"] == pytest.approx(1800.0, rel=1e-9) and ll["theta0"] == pytest.approx(28.0, rel=1e-9)
    ne = impl.fit_exponential(t, T, 0.05)
    assert ne["tau"] == pytest.approx(1800.0, rel=1e-6) and ne["Tinf"] == pytest.approx(TINF, abs=1e-5)
    nf = impl.fit_exponential(t, T, 0.05, Tinf=TINF)
    assert nf["tau"] == pytest.approx(1800.0, rel=1e-8)


def test_loglinear_fit_uncertainty_is_honest_on_seeded_noise():
    rng = np.random.default_rng(11)
    t, T = _exp_data()
    u = 0.05
    pulls, chis = [], []
    for _ in range(40):
        Tn = T + rng.normal(0, u, T.size)
        r = impl.fit_loglinear(t, Tn, TINF, u, theta_min=2.0)
        pulls.append((r["tau"] - 1800.0) / r["u_tau"])
        chis.append(r["chi2"] / r["dof"])
    # pulls should have sd ~ 1 (40 samples: sd known to about 11 %); chi2/dof ~ 1 (sd sqrt(2/dof) ~ 0.09)
    assert 0.7 < np.std(pulls) < 1.35 and abs(np.mean(pulls)) < 0.6
    assert 0.85 < np.mean(chis) < 1.15


def test_probe_lag_model_and_fit():
    tau, tau_p, th0 = 1800.0, 20.0, 28.0
    t = np.arange(0.0, 3600.0 + 1e-9, 5.0)
    Tm = impl.lagged_model(t, th0, tau, tau_p, TINF)
    assert Tm[0] == pytest.approx(TINF)                                   # probe starts at room temperature
    # tau_p dTm/dt = T_body - Tm, checked by central differences on a fine grid
    tf = np.linspace(0.0, 600.0, 60001)
    Tf = impl.lagged_model(tf, th0, tau, tau_p, TINF)
    body = TINF + th0 * np.exp(-tf / tau)
    npt.assert_allclose((tau_p * np.gradient(Tf, tf) + Tf)[1:-1], body[1:-1], atol=1e-5)
    lg = impl.fit_lagged(t, Tm, 0.05, TINF)
    assert lg["tau"] == pytest.approx(tau, rel=1e-5) and lg["tau_p"] == pytest.approx(tau_p, rel=1e-4)
    # once the lag transient has died (t > 5 tau_p) a plain log-linear fit is unbiased in tau
    ll = impl.fit_loglinear(t, Tm, TINF, 0.05, t_min=5 * tau_p * 3)
    assert ll["tau"] == pytest.approx(tau, rel=1e-5)
    assert ll["theta0"] == pytest.approx(th0 * tau / (tau - tau_p), rel=1e-5)   # but the amplitude is scaled


def test_wrong_ambient_biases_tau():
    t, T = _exp_data()
    hi = impl.fit_loglinear(t, T, TINF + 0.5, 0.05, theta_min=1.0)
    lo = impl.fit_loglinear(t, T, TINF - 0.5, 0.05, theta_min=1.0)
    assert hi["tau"] < 1800.0 < lo["tau"]                         # ambient set too high: theta falls too fast
    assert hi["chi2"] > 1.5 * hi["dof"]                           # and chi^2 (zero for the right Tinf) says so


def test_variable_h_fit_recovers_the_convection_factor():
    t = np.arange(0.0, 3600.0 + 1e-9, 15.0)
    sim = impl.simulate_cooling(t, 323.0, TINF, f_conv=1.15)
    vh = impl.fit_variable_h(t, sim["T"], 0.05, TINF, t_min=120.0)
    assert vh["f_conv"] == pytest.approx(1.15, rel=1e-4)
    assert vh["chi2"] < 1e-3 * vh["dof"]


# ---------------------------------------------------------------- distributed model

def test_eigenvalues_solve_their_equations_and_limits():
    for Bi in (1e-3, 0.1, 1.0, 10.0):
        lc = impl.cylinder_eigenvalues(Bi, 6)
        npt.assert_allclose(lc * j1(lc) - Bi * j0(lc), 0.0, atol=1e-10)
        ls = impl.slab_eigenvalues(Bi, 6)
        npt.assert_allclose(ls * np.tan(ls), Bi, rtol=1e-8)
        assert np.all(np.diff(lc) > 0) and np.all(np.diff(ls) > 0)
    assert impl.cylinder_eigenvalues(1e-4, 1)[0] ** 2 == pytest.approx(2e-4, rel=1e-4)     # lumped: lambda^2 = 2 Bi
    assert impl.slab_eigenvalues(1e-4, 1)[0] ** 2 == pytest.approx(1e-4, rel=1e-4)         # lumped: lambda^2 = Bi
    assert impl.cylinder_eigenvalues(1e6, 1)[0] == pytest.approx(2.404826, rel=1e-5)       # Bi -> inf: first zero of J0


def test_series_solutions_conserve_energy_and_start_uniform():
    Bi = 0.5
    Fo = np.array([0.05, 0.2, 0.6, 1.5])
    dF = 1e-5
    dmean = (impl.cylinder_series_mean(Bi, Fo + dF) - impl.cylinder_series_mean(Bi, Fo - dF)) / (2 * dF)
    npt.assert_allclose(dmean, -2 * Bi * impl.cylinder_series(1.0, Bi, Fo), rtol=1e-5)   # rho c V dtheta/dt = -h A theta_s
    # slab: mean = integral of the series over xi in [0, 1]; d mean / dFo = -Bi theta(1)
    xi = np.linspace(0, 1, 2001)
    mean = lambda F: np.trapezoid(impl.slab_series(xi, Bi, np.full_like(xi, F)), xi)
    for F in (0.1, 0.5):
        d = (mean(F + dF) - mean(F - dF)) / (2 * dF)
        assert d == pytest.approx(-Bi * float(impl.slab_series(1.0, Bi, F)), rel=1e-4)
    assert float(impl.cylinder_series(0.0, Bi, 1e-3, n_terms=80)) == pytest.approx(1.0, abs=1e-6)   # centre untouched at first
    assert float(impl.slab_series(0.0, Bi, 1e-3, n_terms=80)) == pytest.approx(1.0, abs=1e-6)


def test_series_reduce_to_lumped_at_small_Bi_and_differ_at_large_Bi():
    Fo = np.linspace(0.5, 200.0, 50)
    Bi = 1e-3
    npt.assert_allclose(impl.cylinder_series_mean(Bi, Fo), np.exp(-2 * Bi * Fo), rtol=2e-3)
    surf, cen = impl.cylinder_series(1.0, 1.0, 1.0), impl.cylinder_series(0.0, 1.0, 1.0)
    assert cen > float(impl.cylinder_series_mean(1.0, 1.0)) > surf          # Bi = 1: strong gradient
    assert abs(float(impl.cylinder_series_mean(1.0, 1.0)) - np.exp(-2.0)) > 0.05


def test_short_cylinder_bar_is_lumped():
    g, C = _bar()
    h = 7.0
    alpha = impl.K_AL / (impl.RHO_AL * impl.C_AL)
    t = np.array([300.0, 1800.0, 3600.0])
    lump = impl.lumped_theta(t, h, impl.RHO_AL, impl.C_AL, g["V"], g["A"])
    cen = impl.short_cylinder(0.0, 0.0, t, h, D, L, impl.K_AL, alpha)
    edge = impl.short_cylinder(1.0, 1.0, t, h, D, L, impl.K_AL, alpha)
    npt.assert_allclose(cen, lump, rtol=1e-3)                              # Bi ~ 2e-4: lumped to 0.1 %
    assert np.all(cen > edge) and np.all(cen - edge < 1e-3)


# ---------------------------------------------------------------- nonlinear lumped model and energy balances

def test_nonlinear_model_energy_balance_and_limits():
    t = np.linspace(0, 3600, 121)
    s = impl.simulate_cooling(t, 323.0, TINF, eps=0.12, f_conv=1.1)
    assert s["C"] * (s["T"][0] - s["T"][-1]) == pytest.approx(s["Q_conv"][-1] + s["Q_rad"][-1], rel=1e-7)
    assert np.all(np.diff(s["T"]) < 0) and s["T"][-1] > TINF
    still = impl.simulate_cooling(t, 323.0, TINF, eps=0.0, f_conv=0.0)
    npt.assert_allclose(still["T"], 323.0, atol=1e-9)                       # no heat path: no cooling
    # near ambient the local time constant equals the linearised prediction
    th = s["T"] - TINF
    i = 100
    tau_local = -th[i] / ((s["T"][i + 1] - s["T"][i - 1]) / (t[i + 1] - t[i - 1]))
    assert tau_local == pytest.approx(impl.predict_tau(th[i], TINF, eps=0.12, f_conv=1.1)["tau"], rel=0.03)


def test_fitted_model_energy_balance_closes():
    _, C = _bar()
    stored, lost = impl.fitted_model_energy_balance(27.5, 1772.0, 3600.0, C)
    assert stored == pytest.approx(lost, rel=1e-12)
    assert impl.en_number(1.0, 0.1, 1.0, 0.2) == 0.0
    assert impl.en_number(1.6, 0.1, 1.0, 0.2) == pytest.approx(0.6 / np.hypot(0.2, 0.4))


# ---------------------------------------------------------------- SYNTHETIC data

def test_synthetic_dataset_is_labelled_seeded_and_plausible(tmp_path):
    a, b = impl.make_synthetic_cooling_data(715), impl.make_synthetic_cooling_data(715)
    assert a["cooling"] == b["cooling"]                                     # seeded
    assert impl.make_synthetic_cooling_data(716)["cooling"] != a["cooling"]
    T = np.array(a["cooling"])[:, 1]
    assert T.max() < impl.T_WATER_MAX_C                                      # never hotter than the water limit
    npt.assert_allclose(T * 10, np.round(T * 10), atol=1e-9)                 # 0.1 C quantization
    tr = a["truth"]
    assert tr["f_conv"] != 1.0 and tr["eps"] != impl.EPS_AL_NOMINAL         # hidden truth differs from nominal
    paths = impl.write_synthetic_csv(tmp_path, seed=715)
    for p in paths:
        assert "SYNTHETIC" in p.name and p.parent == tmp_path
        assert "SYNTHETIC" in p.read_text().splitlines()[0]
    arr = impl.read_csv(paths[0])
    npt.assert_allclose(arr, np.array(a["cooling"]), atol=1e-9)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "cooling_lab.png").stat().st_size > 10_000
