"""Public tests for P10.23 Home Lab: Airfoil Lift & Drag on a Bench Rig. They assert properties, not
magic numbers: room-air density and viscosity reproduce the standard room values; the manometer and
anemometer relations invert and propagate uncertainty correctly; the balance reduction inverts its
forward model (sign of the lift scale, the arm ratio, tares); coefficient uncertainties behave at zero
load; predictions have the right limits (finite span -> 2-D, thickness raises the panel slope, laminar
drag falls as Re^-1/2, Thwaites' separation point does not depend on Re but the transition does);
and the SYNTHETIC dataset is labelled, seeded, matches the shipped files, and reduces back to the
generator's stated truth within the computed uncertainties. Synthetic data are not measurements.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)


# ---------------------------------------------------------------- air and airspeed

def test_room_air_properties():
    assert impl.air_density(20.0, 101.3e3) == pytest.approx(1.204, rel=2e-3)
    assert impl.air_density(20.0, 101.3e3, RH=0.6) < impl.air_density(20.0, 101.3e3)   # water vapour is lighter
    assert impl.air_viscosity(20.0) == pytest.approx(1.81e-5, rel=0.01)
    assert impl.air_viscosity(40.0) > impl.air_viscosity(20.0)


def test_manometer_relation_and_its_uncertainty():
    rho = 1.204
    q, u = impl.manometer_q(1.0, 90.0, rho)                   # vertical tube: q = (rho_w - rho) g h
    assert q == pytest.approx((998.0 - rho) * 9.81 * 1e-3) and u == 0.0
    q10, _ = impl.manometer_q(1.0 / np.sin(np.radians(10)), 10.0, rho)
    assert q10 == pytest.approx(q)                            # inclining stretches the reading by 1/sin
    _, u1 = impl.manometer_q(5.6, 10.0, rho, u_dl_mm=0.3)
    _, u2 = impl.manometer_q(5.6, 10.0, rho, u_dl_mm=0.6)
    assert u2 == pytest.approx(2 * u1)
    _, ua = impl.manometer_q(5.6, 10.0, rho, u_incline_deg=0.2)
    qq, _ = impl.manometer_q(5.6, 10.0, rho)
    assert ua / qq == pytest.approx(np.radians(0.2) / np.tan(np.radians(10.0)), rel=1e-9)


def test_airspeed_round_trip_and_half_relative_uncertainty():
    V, u = impl.airspeed(9.63, 1.204, u_q=0.0963)
    assert 0.5 * 1.204 * V ** 2 == pytest.approx(9.63)
    assert u / V == pytest.approx(0.005)
    q, uq = impl.anemometer_q(4.0, 1.204)
    assert q == pytest.approx(0.5 * 1.204 * 16) and uq / q == pytest.approx(2 * (0.12 + 0.1) / np.sqrt(3) / 4.0)
    assert impl.reynolds_number(4.0, 0.1, 1.204, 1.81e-5) == pytest.approx(1.204 * 0.4 / 1.81e-5)


# ---------------------------------------------------------------- balance

def test_line_fit_and_calibrations_recover_exact_data():
    x = np.array([0.0, 1, 2, 5, 10])
    f = impl.fit_line(x, 18.4 + 2.04 * x)
    assert f["slope"] == pytest.approx(2.04) and f["intercept"] == pytest.approx(18.4) and f["u_slope"] < 1e-10
    R, uR = impl.pivot_calibration(x, 18.4 + 2.04 * x)
    assert R == pytest.approx(2.04)
    m = np.array([0.0, 50, 100, 200])
    k, uk = impl.spring_calibration(m, m / 1000 * 9.81 / 25.0 * 1000)    # 25 N/m spring
    assert k == pytest.approx(25.0)


def test_line_fit_uncertainty_on_seeded_noise():
    # 200 repeats: the spread of fitted slopes should match the reported u_slope within ~15 %
    rng = np.random.default_rng(3)
    x = np.linspace(0, 10, 8)
    fits = [impl.fit_line(x, 2.0 * x + 1 + rng.normal(0, 0.05, x.size)) for _ in range(200)]
    s = np.std([f["slope"] for f in fits], ddof=1)
    assert np.mean([f["u_slope"] for f in fits]) == pytest.approx(s, rel=0.15)


def test_balance_reduction_inverts_the_forward_model():
    g, R = 9.81, 2.04
    L_true, D_true, tare_l, tare_d = 0.05, 0.008, -0.02, 1.05
    off_l, off_d = 352.6, 18.40
    on_l = off_l - L_true / g * 1000 + tare_l
    on_d = off_d + R * D_true / g * 1000 + tare_d
    assert impl.lift_from_scale(off_l, on_l, tare_l) == pytest.approx(L_true)
    assert impl.drag_from_pivot(off_d, on_d, R, tare_d) == pytest.approx(D_true)
    assert impl.lift_from_scale(off_l, off_l - 1.0) > 0                   # the reading falls when the wing lifts
    assert impl.drag_from_spring(10.0, 12.0, 25.0) == pytest.approx(0.05)


def test_coefficient_uncertainty_behaves_at_zero_and_scales():
    C, u = impl.coefficient(0.0, 0.001, 9.6, 0.5, 0.02, 1e-4)
    assert C == 0.0 and u == pytest.approx(0.001 / (9.6 * 0.02))
    C, u = impl.coefficient(0.1, 0.0, 10.0, 0.5, 0.02)
    assert C == pytest.approx(0.5) and u / C == pytest.approx(0.05)


def test_blockage_estimate_is_small_and_grows_with_alpha():
    e = [impl.blockage_fraction(a) for a in (0, 5, 10, 15)]
    assert np.all(np.diff(e) > 0) and 0 < e[0] < e[-1] < 0.05
    assert impl.blockage_fraction(-10) == pytest.approx(impl.blockage_fraction(10))


# ---------------------------------------------------------------- predictions

def test_lift_predictions_have_the_right_limits():
    assert impl.thin_airfoil_cl(5.0) == pytest.approx(2 * np.pi * np.radians(5.0))
    cl_p = impl.panel_cl("0012", [0.0, 4.0])
    assert abs(cl_p[0]) < 1e-12
    assert cl_p[1] / impl.thin_airfoil_cl(4.0) == pytest.approx(1 + 0.77 * 0.12, abs=0.02)
    assert impl.finite_span_slope(2 * np.pi, 1e9) == pytest.approx(2 * np.pi)
    a = [impl.finite_span_slope(2 * np.pi, A) for A in (2, 4, 8)]
    assert np.all(np.diff(a) > 0) and a[1] == pytest.approx(2 * np.pi / (1 + 2 * np.pi / (4 * np.pi)))


def test_low_re_drag_estimate_scaling_and_symmetry():
    cd1 = impl.low_re_cd(2.5e4, 0.0, 0.0)
    cd4 = impl.low_re_cd(1.0e5, 0.0, 0.0)
    assert cd1 / cd4 == pytest.approx(2.0)                                 # laminar friction ~ Re^-1/2
    assert impl.low_re_cd(3e4, 0.4, 5.0) == pytest.approx(impl.low_re_cd(3e4, -0.4, -5.0))
    assert impl.low_re_cd(3e4, 0.4, 5.0, suction_loss=1.0) > impl.low_re_cd(3e4, 0.4, 5.0, suction_loss=0.0)
    assert impl.low_re_cd(3e4, 0.4, 5.0, AR_eff=4.0) > impl.low_re_cd(3e4, 0.4, 5.0)


def test_thwaites_separation_is_reynolds_independent_but_transition_is_not():
    lo = impl.laminar_upper_surface("0012", 2.0, 2.6e4)
    hi = impl.laminar_upper_surface("0012", 2.0, 2.0e6)
    assert lo["x_sep"] == pytest.approx(hi["x_sep"])
    assert hi["transition_first"] and not lo["transition_first"]
    assert hi["Re_theta_sep"] / lo["Re_theta_sep"] == pytest.approx(np.sqrt(2.0e6 / 2.6e4), rel=1e-6)
    seps = [impl.laminar_upper_surface("0012", a, 2.6e4)["x_sep"] for a in (0.0, 2.0, 4.0, 6.0)]
    assert np.all(np.diff(seps) < 0)                                       # separation creeps forward with alpha


def test_prediction_table_is_consistent():
    p = impl.predict_bench([2.0, 6.0], 4.0)
    S = impl.CHORD * impl.SPAN
    npt.assert_allclose(p["lift_g_lo"], p["CL_lo"] * p["q"] * S / 9.81 * 1000)
    assert np.all(p["CL_lo"] < p["CL_hi"]) and np.all(p["CD_lo"] <= p["CD_hi"])
    assert np.all(p["CL_hi"] < p["cl_thin"])
    assert p["Re"] == pytest.approx(impl.reynolds_number(4.0, impl.CHORD))


def test_ladson_comparison_returns_unity_on_its_own_data():
    d = impl.load_ladson("I", 2.0e6)
    ok = np.isfinite(d["cd"]) & (d["alpha_deg"] >= 0)
    c = impl.compare_with_ladson(d["alpha_deg"][ok], d["cl"][ok], d["cd"][ok])
    npt.assert_allclose(c["CL_ratio"][1:], 1.0)
    npt.assert_allclose(c["CD_ratio"], 1.0)


# ---------------------------------------------------------------- SYNTHETIC data

def test_synthetic_dataset_is_labelled_seeded_and_matches_the_shipped_files(tmp_path):
    paths = impl.write_synthetic_csv(tmp_path, seed=230)
    assert len(paths) == 5
    shipped = impl.find_repo_root() / "projects" / "p10-airfoil-bench" / "data" / "synthetic"
    for p in paths:
        assert p.name.startswith("SYNTHETIC_seed230_")
        first = p.read_text().splitlines()[0]
        assert first.startswith("# SYNTHETIC") and "not a measurement" in first
        assert (shipped / p.name).read_text() == p.read_text()
    other = impl.write_synthetic_csv(tmp_path / "s7", seed=7)
    assert other[3].read_text() != paths[3].read_text()


def test_reduction_recovers_the_generator_truth(tmp_path):
    impl.write_synthetic_csv(tmp_path, seed=230)
    ds = impl.load_bench_dataset(tmp_path, seed=230)
    r = impl.reduce_bench(ds)
    truth = impl.make_synthetic_bench_data(230)["truth"]
    assert r["ratio"] == pytest.approx(truth["ratio"], abs=3 * r["u_ratio"] + 1e-3)
    assert r["flow_angle_deg"] == pytest.approx(-truth["flow_angle_deg"], abs=0.25)
    assert r["slope_per_deg"] == pytest.approx(0.060, rel=0.1)          # generator's stated slope
    assert r["u_slope_total"] > r["u_slope"] > 0                        # q's uncertainty scales every CL
    assert 3.6 < r["V"] < 4.4 and 2.0e4 < r["Re"] < 3.2e4
    CD = np.array([x["CD"] for x in r["rows"]])
    assert np.all((CD > 0.02) & (CD < 0.12))
    assert all(x["u_CL"] > 0 and x["u_CD"] > 0 for x in r["rows"])


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "airfoil_bench.png").stat().st_size > 10_000
