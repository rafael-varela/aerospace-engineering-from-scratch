"""Public tests for P08.19 Home Lab: Flow & Pressure Measurement.

They assert properties, not one magic number: the Torricelli level satisfies the mass-and-Bernoulli
ODE and its initial condition; the draining time scales as 1/C_d, (D/d)^2 and 1/sqrt(g); first-order
uncertainty propagation agrees with the power-law rule and with a seeded Monte Carlo; the weighted fit
recovers an exact line and has a chi^2 consistent with seeded noise; C_d is recovered from noise-free
levels and, within its uncertainty, from the SYNTHETIC data; the air density reproduces the ISA
sea-level value; manometer, tilt and airspeed formulas invert each other; the resolution limit gives
exactly the target uncertainty; and the SYNTHETIC files shipped in projects/p08-flow-pressure match the
seeded generator, say SYNTHETIC, and read back with the documented columns.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
from pathlib import Path

import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

G = impl.G
D, d, CD = 0.090, 0.004, 0.61
KIT = Path(__file__).resolve().parents[4] / "projects" / "p08-flow-pressure" / "data" / "synthetic"


# ---------------------------------------------------------------- experiment 1: Torricelli

def test_torricelli_level_satisfies_the_draining_ode():
    t = np.linspace(0.0, 120.0, 2001)
    h = impl.torricelli_level(t, 0.2, CD, d, D)
    assert h[0] == pytest.approx(0.2)
    dhdt = np.gradient(h, t)
    A, a = np.pi * D ** 2 / 4, impl.orifice_area(d)
    npt.assert_allclose(A * dhdt[1:-1], -CD * a * np.sqrt(2 * G * h[1:-1]), rtol=1e-4)   # mass + Bernoulli
    assert impl.orifice_area(d) == pytest.approx(np.pi * 4e-6)
    T_empty = impl.drain_time(0.2, 0.0, CD, d, D)
    assert np.all(impl.torricelli_level(np.array([T_empty, T_empty + 10]), 0.2, CD, d, D) == 0.0)


def test_drain_time_consistency_and_scaling():
    T = impl.drain_time(0.2, 0.03, CD, d, D)
    assert impl.torricelli_level(T, 0.2, CD, d, D) == pytest.approx(0.03, rel=1e-10)
    assert impl.drain_time(0.2, 0.03, 2 * CD, d, D) == pytest.approx(T / 2)
    assert impl.drain_time(0.2, 0.03, CD, d, 2 * D) == pytest.approx(4 * T)
    assert impl.drain_time(0.2, 0.03, CD, 2 * d, D) == pytest.approx(T / 4)
    assert impl.drain_time(0.2, 0.03, CD, d, D, g=4 * G) == pytest.approx(T / 2)     # time ~ 1/sqrt(g)
    assert impl.drain_time(0.2, 0.2, CD, d, D) == 0.0


def test_prediction_budget_follows_the_power_law_rule_and_monte_carlo():
    u_cd, u_d, u_D = 0.03, 0.1e-3, 0.5e-3
    pr = impl.predict_drain_time(0.2, 0.03, CD, u_cd, d, u_d, D, u_D)
    rel = np.sqrt((u_cd / CD) ** 2 + (2 * u_d / d) ** 2 + (2 * u_D / D) ** 2)
    assert pr["u_T"] / pr["T"] == pytest.approx(rel, rel=1e-6)
    assert pr["budget"]["d"] == pytest.approx(pr["T"] * 2 * u_d / d, rel=1e-6)
    rng = np.random.default_rng(1919)
    n = 40_000                     # MC standard error of the sd is about 1/sqrt(2n) = 0.35 %; allow 3 %
    Ts = impl.drain_time(0.2, 0.03, rng.normal(CD, u_cd, n), rng.normal(d, u_d, n), rng.normal(D, u_D, n))
    assert np.std(Ts) == pytest.approx(pr["u_T"], rel=0.03)


def test_weighted_line_fit_exact_and_chi2_on_seeded_noise():
    x = np.linspace(0, 10, 11)
    f = impl.fit_line(x, 3.0 - 0.25 * x, 0.1)
    assert f["slope"] == pytest.approx(-0.25) and f["intercept"] == pytest.approx(3.0)
    assert f["chi2"] < 1e-20 and f["dof"] == 9
    rng = np.random.default_rng(7)
    chi2 = [impl.fit_line(x, 3.0 - 0.25 * x + rng.normal(0, 0.1, x.size), 0.1)["chi2"] for _ in range(400)]
    assert np.mean(chi2) == pytest.approx(9.0, rel=0.1)   # E[chi2] = dof; sd of the mean ~ sqrt(18/400) = 0.21


def test_discharge_coefficient_recovered_from_noise_free_levels():
    t = np.linspace(2.0, 100.0, 18)
    h = impl.torricelli_level(t, 0.205, 0.63, d, D)
    r = impl.fit_discharge_coefficient(t, h, 0.5e-3, 0.08, d, 0.1e-3, D, 0.5e-3)
    assert r["cd"] == pytest.approx(0.63, rel=1e-9)
    assert r["h0"] == pytest.approx(0.205, rel=1e-9)
    assert r["rel_budget"]["d"] == pytest.approx(2 * 0.1e-3 / d)
    assert r["u_cd"] / r["cd"] == pytest.approx(np.sqrt(sum(v ** 2 for v in r["rel_budget"].values())))


def test_discharge_coefficient_from_the_synthetic_data_is_consistent_with_its_truth():
    # SYNTHETIC truth: C_d = 0.62 with d = 4.05 mm and D = 90.4 mm; read through the nominal 4.0 and 90.0 mm
    # the apparent C_d is 0.62 (4.05/4.0)^2 (90.0/90.4)^2 = 0.630; the fit must agree within 2 u
    data = impl.make_synthetic_flow_data(819)["drain"]
    apparent = 0.62 * (4.05 / 4.0) ** 2 * (90.0 / 90.4) ** 2
    for run in (1, 2, 3):
        m = data[:, 0] == run
        r = impl.fit_discharge_coefficient(data[m, 2], data[m, 1] * 1e-3, 0.5e-3, 0.08, d, 0.1e-3, D, 0.5e-3)
        assert abs(r["cd"] - apparent) < 2 * r["fit"]["u_slope"] / r["k"] * r["cd"] + 0.003
        assert 0.3 < r["fit"]["chi2"] / r["fit"]["dof"] < 2.5


def test_jet_reynolds_number():
    re = impl.jet_reynolds(np.array([0.01, 0.04, 0.2]), d)
    assert re[1] == pytest.approx(2 * re[0])                               # V ~ sqrt(h)
    assert re[2] == pytest.approx(impl.RHO_WATER * np.sqrt(2 * G * 0.2) * d / impl.MU_WATER)
    assert np.all(re > 1500)


# ---------------------------------------------------------------- experiment 2: manometer and airspeed

def test_air_density_reproduces_isa_sea_level():
    assert impl.air_density(101325.0, 288.15) == pytest.approx(1.225, abs=5e-4)


def test_manometer_tilt_and_airspeed_invert_each_other():
    rho = 1.2
    V = np.array([2.0, 5.0, 12.0])
    dp = 0.5 * rho * V ** 2
    dh = dp / ((impl.RHO_WATER - rho) * G)
    npt.assert_allclose(impl.manometer_dp(dh, rho), dp, rtol=1e-12)
    npt.assert_allclose(impl.airspeed(impl.manometer_dp(dh, rho), rho), V, rtol=1e-12)
    assert impl.vertical_head(0.010, 90.0) == pytest.approx(0.010)
    assert impl.vertical_head(0.010, 30.0) == pytest.approx(0.005)
    assert np.isnan(impl.airspeed(-1.0, rho))
    assert impl.airspeed(0.0, rho) == 0.0


def test_airspeed_uncertainty_and_the_tilted_tube():
    rho = 1.2
    V, uV = impl.airspeed_with_uncertainty(0.009, 0.0006, 90.0, 0.0, rho, 0.0)
    assert uV / V == pytest.approx(0.5 * 0.0006 / 0.009)                  # V ~ sqrt(dh): half the relative u
    # the same flow read on a tube tilted at 10 deg: 1/sin(10) times the reading, same absolute u_L
    L10 = 0.009 / np.sin(np.radians(10))
    V10, uV10 = impl.airspeed_with_uncertainty(L10, 0.0006, 10.0, 0.0, rho, 0.0)
    assert V10 == pytest.approx(V) and uV10 == pytest.approx(uV * np.sin(np.radians(10)))
    # an uncertain tilt matters more at shallow angles: u_theta/tan(theta)
    _, uA = impl.airspeed_with_uncertainty(L10, 0.0, 10.0, 0.5, rho, 0.0)
    assert uA / V10 == pytest.approx(0.5 * np.radians(0.5) / np.tan(np.radians(10)), rel=1e-9)


def test_resolution_limit_gives_exactly_the_target_uncertainty():
    rho = 1.196
    for u_dh, target in ((0.6e-3, 0.10), (0.1e-3, 0.05)):
        Vmin = impl.min_speed_for_uncertainty(u_dh, target, rho)
        dh = 0.5 * rho * Vmin ** 2 / ((impl.RHO_WATER - rho) * G)
        assert u_dh / (2 * dh) == pytest.approx(target, rel=1e-10)
    v1 = impl.min_speed_for_uncertainty(0.4e-3, 0.1, rho)
    v4 = impl.min_speed_for_uncertainty(1.6e-3, 0.1, rho)
    assert v4 == pytest.approx(2 * v1)                                     # V_min ~ sqrt(u_dh)


def test_reynolds_fan_flow_and_en():
    assert impl.reynolds(1.2, 3.0, 0.006, 1.8e-5) == pytest.approx(1200.0)
    assert impl.fan_mean_velocity(3600.0 * np.pi / 4, 1.0) == pytest.approx(1.0)
    assert impl.fan_mean_velocity(1100.0, 0.28, 0.08) > impl.fan_mean_velocity(1100.0, 0.28)
    assert impl.en_number(1.0, 0.1, 1.0, 0.2) == 0.0
    assert impl.en_number(1.0, 0.3, 0.0, 0.4) == pytest.approx(1.0)        # |1| / (2 * 0.5)


# ---------------------------------------------------------------- SYNTHETIC data

def test_shipped_synthetic_files_match_the_generator_and_are_labelled(tmp_path):
    shipped = sorted(KIT.glob("SYNTHETIC_seed819_*.csv"))
    assert [p.name for p in shipped] == ["SYNTHETIC_seed819_drain.csv", "SYNTHETIC_seed819_manometer.csv"]
    fresh = impl.write_synthetic_csv(tmp_path, 819)
    for p in fresh:
        assert "SYNTHETIC" in p.name and "SYNTHETIC" in p.read_text().splitlines()[0]
        assert p.read_text() == (KIT / p.name).read_text()                  # regenerated byte for byte
    dr = impl.read_csv(KIT / "SYNTHETIC_seed819_drain.csv")
    assert list(dr) == ["run", "mark_mm", "t_s"] and dr["t_s"].size == 54
    ma = impl.read_csv(KIT / "SYNTHETIC_seed819_manometer.csv")
    assert list(ma) == ["source", "setting", "distance_mm", "angle_deg", "repeat", "zero_left_mm",
                        "zero_right_mm", "left_mm", "right_mm", "anemometer_mps"]
    assert set(ma["angle_deg"]) == {90.0, 10.0} and ma["left_mm"].size == 16
    a, b = impl.make_synthetic_flow_data(5), impl.make_synthetic_flow_data(5)
    npt.assert_array_equal(a["drain"], b["drain"])                          # seeded: reproducible


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "flow_lab.png").stat().st_size > 10_000
