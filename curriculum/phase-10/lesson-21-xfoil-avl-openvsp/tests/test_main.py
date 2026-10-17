"""Public tests for P10.21 Using XFOIL, AVL & OpenVSP. Nothing here imports or runs XFOIL, AVL or
OpenVSP. The polar files used by these tests are SYNTHETIC fixtures that the tests write into
tmp_path themselves (their headers say so); they test the reader, not any tool. The tests check that
the reader refuses unlabelled files and files missing a tool's settings; that the XFOIL-layout
parser extracts the header settings; that the metrics recover a known line exactly; that comparing a
polar with itself gives zero and a shifted copy gives the shift; and that the learner's own panel and
boundary-layer codes behave correctly (symmetry, thin-airfoil limit, second-order panel convergence)
and sit where they should against the published Ladson NACA 0012 data (NASA TM-4074).

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

FIXTURE_NOTE = "SYNTHETIC test fixture written by test_main.py; invented numbers, not tool output"

XFOIL_LAYOUT_FIXTURE = """
       XFOIL         Version 6.99

 Calculated polar for: NACA 0012 SYNTHETIC TEST FIXTURE (invented numbers, not XFOIL output)

 1 1 Reynolds number fixed          Mach number fixed

 xtrf =   1.000 (top)        1.000 (bottom)
 Mach =   0.150     Re =     6.000 e 6     Ncrit =   9.000

  alpha    CL        CD       CDp       CM     Top_Xtr  Bot_Xtr
 ------- -------- --------- --------- -------- -------- --------
  -2.000  -0.2200   0.00600   0.00100   0.0010   0.6000   0.3000
   0.000   0.0000   0.00550   0.00090   0.0000   0.4500   0.4500
   2.000   0.2200   0.00600   0.00100  -0.0010   0.3000   0.6000
   4.000   0.4400   0.00700   0.00130  -0.0020   0.1500   0.7500
"""


def _write_fixture(path, meta, a=(-2, 0, 2, 4, 6), slope=0.11, a0=-1.0):
    a = np.asarray(a, float)
    cl = slope * (a - a0)
    cd = 0.006 + 0.004 * cl ** 2
    return impl.write_tool_polar(path, meta, a, cl, cd, -0.05 + 0 * a)


def _good_meta():
    return {"tool": "XFOIL", "version": "SYNTHETIC-0.0", "note": FIXTURE_NOTE,
            "settings": {"re": "6e6", "mach": 0.15, "ncrit": 9, "panels": 160, "iter": 100}}


# ---------------------------------------------------------------- the reader refuses unlabelled files

def test_reader_refuses_a_file_without_a_header(tmp_path):
    p = tmp_path / "unlabelled.csv"
    p.write_text("alpha_deg,cl,cd,cm\n0,0.1,0.006,0\n2,0.3,0.007,0\n")
    with pytest.raises(ValueError):
        impl.read_tool_polar(p)


def test_reader_refuses_missing_version_or_settings(tmp_path):
    m = _good_meta(); del m["version"]
    with pytest.raises(ValueError):
        impl.read_tool_polar(_write_fixture(tmp_path / "nover.csv", m))
    m = _good_meta(); del m["settings"]
    with pytest.raises(ValueError):
        impl.read_tool_polar(_write_fixture(tmp_path / "noset.csv", m))


def test_reader_refuses_a_tool_polar_missing_a_required_setting(tmp_path):
    m = _good_meta(); del m["settings"]["ncrit"]
    with pytest.raises(ValueError):
        impl.read_tool_polar(_write_fixture(tmp_path / "noncrit.csv", m))
    m = {"tool": "AVL", "version": "SYNTHETIC-0.0", "note": FIXTURE_NOTE, "settings": {"mach": 0.0, "nchord": 8}}
    with pytest.raises(ValueError):
        impl.read_tool_polar(_write_fixture(tmp_path / "avl.csv", m))


def test_reader_accepts_a_labelled_file_and_round_trips(tmp_path):
    p = _write_fixture(tmp_path / "SYNTHETIC_fixture_xfoil.csv", _good_meta())
    r = impl.read_tool_polar(p)
    assert r["meta"]["tool"] == "XFOIL" and r["meta"]["version"] == "SYNTHETIC-0.0"
    assert r["meta"]["settings"]["ncrit"] == "9" and r["meta"]["settings"]["panels"] == "160"
    assert "SYNTHETIC" in r["meta"]["note"]
    npt.assert_allclose(r["cl"], 0.11 * (r["alpha_deg"] + 1.0), rtol=1e-5)
    assert r["alpha_deg"].shape == (5,)


def test_xfoil_layout_parser_reads_the_header_settings():
    p = impl.parse_xfoil_polar(XFOIL_LAYOUT_FIXTURE)
    m = p["meta"]
    assert m["version"] == "6.99" and m["re"] == pytest.approx(6.0e6)
    assert m["mach"] == pytest.approx(0.15) and m["ncrit"] == pytest.approx(9.0)
    assert m["xtr_top"] == pytest.approx(1.0) and "NACA 0012" in m["airfoil"]
    npt.assert_allclose(p["alpha_deg"], [-2, 0, 2, 4])
    npt.assert_allclose(p["cd"], [0.006, 0.0055, 0.006, 0.007])
    with pytest.raises(ValueError):
        impl.parse_xfoil_polar("alpha CL CD\n0 0 0\n")


def test_xfoil_conversion_adds_the_settings_xfoil_does_not_write(tmp_path):
    path = impl.xfoil_polar_to_labelled_csv(XFOIL_LAYOUT_FIXTURE, tmp_path / "SYNTHETIC_fixture_converted.csv",
                                            panels=160, iter_limit=100, produced_by="pytest fixture")
    r = impl.read_tool_polar(path)
    assert r["meta"]["settings"]["panels"] == "160" and r["meta"]["settings"]["iter"] == "100"
    assert float(r["meta"]["settings"]["re"]) == pytest.approx(6e6)
    npt.assert_allclose(r["cl"], [-0.22, 0.0, 0.22, 0.44])


# ---------------------------------------------------------------- metrics

def test_metrics_recover_a_known_line_exactly():
    a = np.array([-4, -2, 0, 2, 4, 6, 8, 12, 14.0])
    cl = 0.105 * (a + 2.0)
    cl[-2:] = [1.0, 0.9]                                   # stall points outside the linear range
    cd = 0.007 + 0.01 * (a / 10) ** 2
    m = impl.polar_metrics(a, cl, cd, (-4, 8))
    assert m["cl_alpha_deg"] == pytest.approx(0.105)
    assert m["cl_alpha_rad"] == pytest.approx(0.105 * 180 / np.pi)
    assert m["alpha_L0_deg"] == pytest.approx(-2.0)
    assert m["rms_line"] < 1e-12 and m["cd_min"] == pytest.approx(0.007) and m["alpha_cd_min"] == 0.0
    assert m["cl_max"] == pytest.approx(1.05)


def test_comparing_a_polar_with_itself_and_with_a_shifted_copy():
    a = np.arange(-4, 9, 1.0)
    p = {"alpha_deg": a, "cl": 0.1 * a, "cd": 0.006 + 0.0 * a}
    c = impl.compare_polars(p, p)
    assert c["rms_dcl"] == 0.0 and c["d_alpha_L0_deg"] == pytest.approx(0.0, abs=1e-12)
    q = {"alpha_deg": a + 0.3, "cl": 0.1 * a, "cd": 0.0066 + 0.0 * a}   # same curve shifted by 0.3 deg
    c = impl.compare_polars(p, q)
    assert c["d_alpha_L0_deg"] == pytest.approx(0.3)
    assert c["rms_dcl"] == pytest.approx(0.03)
    assert c["rel_d_cl_alpha"] == pytest.approx(0.0, abs=1e-12)
    assert c["rel_d_cd_min"] == pytest.approx(0.1)


# ---------------------------------------------------------------- authentic data and the own codes

def test_ladson_data_is_found_and_read_as_published():
    d = impl.load_ladson("I", 5.97e6)
    assert len(d["alpha_deg"]) == 18 and "TM-4074" in d["source"]
    m = impl.polar_metrics(d["alpha_deg"], d["cl"], d["cd"], (-4.5, 8.5))
    assert abs(m["alpha_L0_deg"]) < 0.1                   # a symmetric section
    assert 0.09 < m["cl_alpha_deg"] < 0.12                # below thin-airfoil 2 pi = 0.1097/deg plus thickness
    assert np.isnan(d["cd"][-1])                          # post-stall rows print no drag
    assert impl.find_phase10_data().name == "data"


def test_panel_code_symmetry_thin_limit_and_convergence():
    X, Y = impl.naca4("0012", 60)
    assert abs(impl.vortex_panel(X, Y, 0.0)["cl"]) < 1e-12
    X, Y = impl.naca4("0001", 100)
    cl = impl.vortex_panel(X, Y, np.radians(2.0))["cl"]
    assert cl == pytest.approx(2 * np.pi * np.radians(2.0), rel=0.02)
    cls = [impl.vortex_panel(*impl.naca4("0012", n), np.radians(4.0))["cl"] for n in (25, 50, 100, 200)]
    e = np.abs(np.diff(cls))
    p = np.log(e[:-1] / e[1:]) / np.log(2)
    assert np.all(np.abs(p - 2) < 0.2)


def test_boundary_layer_estimate_responds_to_reynolds_number_and_tripping():
    X, Y = impl.naca4("0012", 100)
    p = impl.vortex_panel(X, Y, 0.0)
    free_lo = impl.bl_drag(p, 2e6)["cd"]
    free_hi = impl.bl_drag(p, 6e6)["cd"]
    tripped = impl.bl_drag(p, 6e6, (0.05, 0.05))["cd"]
    assert free_hi < free_lo
    assert tripped > free_hi
    assert impl.bl_drag(p, 6e6)["upper"]["reason"] == "Michel"


def test_own_codes_against_ladson_sit_where_the_physics_says():
    v = impl.validate_against_ladson(5.97e6, n=80)
    c = v["compare"]
    # inviscid panel lift is too steep: the real boundary layer decambers the section (about 5-20 %)
    assert 0.05 < c["rel_d_cl_alpha"] < 0.20
    assert abs(c["d_alpha_L0_deg"]) < 0.1
    # the one-way drag estimate is within 25 % of the measured minimum drag at this Re
    assert abs(c["rel_d_cd_min"]) < 0.25
    assert np.all((v["cd_ratio"] > 0.7) & (v["cd_ratio"] < 1.4))


def test_own_polar_carries_its_own_label_and_round_trips(tmp_path):
    mine = impl.own_polar("0012", [0.0, 2.0, 4.0], 6e6, n=60)
    assert {"tool", "version", "settings"} <= set(mine["meta"])
    path = impl.write_tool_polar(tmp_path / "own_naca0012_Re6e6.csv", mine["meta"], mine["alpha_deg"], mine["cl"], mine["cd"])
    r = impl.read_tool_polar(path)
    npt.assert_allclose(r["cl"], mine["cl"], rtol=1e-5)


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path, quick=True)
    assert (tmp_path / "tool_validation.png").stat().st_size > 10_000
