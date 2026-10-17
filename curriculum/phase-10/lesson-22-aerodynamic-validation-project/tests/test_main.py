"""Public tests for P10.22 Aerodynamic Validation Project. They assert properties of each layer of a
V&V study, not one magic number. Code verification: a symmetric wing's span load is symmetric, the
half-wing with its symmetry image equals the full wing, lift is linear in sin(alpha), a very long
unswept wing tends to 2 pi, and the lattice converges monotonically at a positive observed order.
Solution verification: the three-grid observed order, Richardson extrapolation and GCI formulas
recover a manufactured solution exactly. Input and validation: the section slope read from Ladson's
published NACA 0012 data is physically sensible; the section correction behaves as derived; the
uncertainty combination is right; and against the twelve measured planforms the VLM is good on
average and fails worst where the linear model should (the AR 1.34 wing with nonlinear lift).

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`.
"""
import numpy as np
import pytest

from conftest import load_impl

impl = load_impl(__file__)

COARSE = ((1, 4), (2, 8), (4, 16))      # keeps the tests fast; main() uses finer lattices


# ---------------------------------------------------------------- code verification

def test_symmetric_wing_has_a_symmetric_span_load():
    assert impl.symmetry_residual(2.61, 0.25, 45.0) < 1e-12
    assert impl.symmetry_residual(5.16, 1.0, 0.0) < 1e-12


def test_half_wing_with_image_equals_full_wing():
    for args in ((5.16, 1.0, 45.0), (2.61, 0.5, 45.0)):
        assert impl.vlm_cla(*args, 3, 6, half=True) == pytest.approx(impl.vlm_cla(*args, 3, 6, half=False), rel=1e-10)


def test_lift_is_linear_in_sin_alpha():
    lat = impl.trapezoidal_lattice(5.16, 1.0, 30.0, 3, 8)
    c1 = impl.lift_coefficient(lat, impl.solve_vlm(lat, np.radians(2.0)))
    c2 = impl.lift_coefficient(lat, impl.solve_vlm(lat, np.radians(7.0)))
    assert c2 / c1 == pytest.approx(np.sin(np.radians(7.0)) / np.sin(np.radians(2.0)), rel=1e-12)


def test_long_unswept_wing_tends_to_two_pi_and_kappa_two_pi():
    assert impl.vlm_cla(1000.0, 1.0, 0.0, 2, 40) == pytest.approx(2 * np.pi, rel=0.01)
    assert impl.vlm_cla(1000.0, 1.0, 0.0, 2, 40, kappa=0.9) == pytest.approx(0.9 * 2 * np.pi, rel=0.01)


def test_lattice_converges_monotonically_with_positive_order():
    vals = [impl.vlm_cla(5.16, 1.0, 0.0, nc, ns) for nc, ns in ((1, 3), (2, 6), (4, 12), (8, 24))]
    assert np.all(np.diff(vals) < 0)
    p = impl.observed_order_3(vals[3], vals[2], vals[1])
    assert 0.5 < p < 1.5                              # the tip singularity limits a uniform lattice to ~first order


# ---------------------------------------------------------------- solution verification formulas

def test_three_grid_formulas_recover_a_manufactured_solution():
    f_exact, C, p = 3.0, 0.4, 1.7
    f = [f_exact + C * h ** p for h in (0.1, 0.2, 0.4)]          # fine, medium, coarse
    assert impl.observed_order_3(*f) == pytest.approx(p, rel=1e-10)
    assert impl.richardson(f[0], f[1], p) == pytest.approx(f_exact, rel=1e-12)
    g = impl.gci_fine(f[0], f[1], p)
    assert g == pytest.approx(1.25 * abs(f[1] - f[0]) / f[0] / (2 ** p - 1))
    assert abs(f[0] - f_exact) / f_exact < g                     # the band covers the true error


def test_grid_study_reports_a_small_positive_gci():
    gs = impl.grid_study(2.61, 1.0, 0.0, grids=COARSE)
    assert 0 < gs["gci"] < 0.1 and gs["u_num"] > 0
    f1, f2, f3 = gs["f"]
    assert abs(gs["extrapolated"] - f1) < abs(f1 - f3)


# ---------------------------------------------------------------- inputs and data

def test_section_slope_from_ladson_is_sensible():
    a0, se = impl.section_slope_from_ladson(2.0e6)
    a6, se6 = impl.section_slope_from_ladson(5.97e6)
    assert 5.5 < a0 < 2 * np.pi * (1 + 0.77 * 0.12)              # below inviscid thin-section value
    assert a6 > a0                                               # thinner boundary layer at higher Re
    assert 0 < se < 0.1 and 0 < se6 < 0.1


def test_finite_wing_data_are_read_with_their_reading_uncertainty():
    rows = impl.load_finite_wing_data()
    assert len(rows) == 15
    assert {round(r["u_read_per_deg"], 4) for r in rows} == {0.0005, 0.0003}
    assert all(0.02 < r["cla_per_deg"] < 0.1 for r in rows)
    assert sum(r["report"] == "NACA TN-1669" for r in rows) == 9


def test_section_correction_lowers_the_slope_less_than_proportionally():
    c1 = impl.vlm_cla(5.16, 1.0, 0.0, 2, 8, kappa=1.0)
    c9 = impl.vlm_cla(5.16, 1.0, 0.0, 2, 8, kappa=0.9)
    assert 0.9 * c1 < c9 < c1          # a finite wing is less sensitive to its section slope than a 2-D one


def test_validation_uncertainty_combination():
    assert impl.validation_uncertainty(0.05, 0.0, 0.0, 0.0, scatter_rel=0.0) == 0.0
    u = impl.validation_uncertainty(0.05, 0.0005, 0.0003, 0.0004)
    expected = np.sqrt(0.0003 ** 2 + 0.0004 ** 2 + (0.0005 / np.sqrt(3)) ** 2 + (0.03 * 0.05 / np.sqrt(3)) ** 2)
    assert u == pytest.approx(expected)


def test_validation_metrics_on_a_constructed_table():
    tab = [{"rel_err": 0.02, "within": True, "planform": "a"}, {"rel_err": -0.10, "within": False, "planform": "b"},
           {"rel_err": 0.04, "within": True, "planform": "c"}, {"rel_err": 0.0, "within": True, "planform": "d"}]
    m = impl.validation_metrics(tab)
    assert m["mean_abs_rel"] == pytest.approx(0.04) and m["max_abs_rel"] == pytest.approx(0.10)
    assert m["worst"] == "b" and m["fraction_within"] == pytest.approx(0.75) and m["mean_rel"] == pytest.approx(-0.01)


def test_validation_against_the_twelve_planforms():
    a0, _ = impl.section_slope_from_ladson(2.0e6)
    a6, _ = impl.section_slope_from_ladson(5.97e6)
    tab = impl.validation_table(a0 / (2 * np.pi), abs(a6 - a0) / (2 * np.pi), grids=COARSE)
    assert len(tab) == 15
    m = impl.validation_metrics(tab)
    assert m["mean_abs_rel"] < 0.08 and m["fraction_within"] >= 0.6
    worst = max(tab, key=lambda t: abs(t["rel_err"]))
    assert worst["AR"] == pytest.approx(1.34) and worst["sweep"] == 0.0 and worst["rel_err"] > 0
    assert all(t["u_val"] > 0 for t in tab)


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path, quick=True)
    assert (tmp_path / "vlm_validation.png").stat().st_size > 10_000
