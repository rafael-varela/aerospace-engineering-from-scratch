"""Public tests for P07.03 Properties of Pure Substances. They assert what a property table must satisfy
and what a lookup must do, not magic numbers: interpolation reproduces the table nodes; the logarithmic
p_sat interpolation beats the linear one; the saturation temperature inverts the saturation pressure;
phase identification and quality are consistent; Clausius-Clapeyron dp/dT = h_fg/(T v_fg) holds on the
table; h = u + p v to table round-off; Z -> 1 as p -> 0; superheated-steam enthalpy differences at low
pressure match the NASA ideal-gas H2O fit; and the rigid-tank first law and fill limit are consistent.

Data: NIST SRD 69 values as tabulated in Yan, Introduction to Engineering Thermodynamics (LibreTexts).
Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. The only file written goes to tmp_path.
"""
import numpy as np
import pytest

from conftest import load_impl

impl = load_impl(__file__)

FLUIDS = ["water", "ammonia"]


# ---------------------------------------------------------------- saturation lookup

@pytest.mark.parametrize("fl", FLUIDS)
def test_saturation_lookup_reproduces_the_table_nodes(fl):
    t = impl.sat_table(fl)
    for i in (0, 5, len(t["T"]) // 2, len(t["T"]) - 2):
        st = impl.saturation_state(fl, t["T"][i])
        for k in ("vf", "vg", "uf", "ug", "hf", "hg", "sf", "sg"):
            assert st[k] == pytest.approx(t[k][i], rel=1e-12)
        for meth in ("linear", "log"):
            assert impl.saturation_pressure(fl, t["T"][i], meth) == pytest.approx(t["p"][i], rel=1e-12)
        assert st["hfg"] == pytest.approx(t["hg"][i] - t["hf"][i])
    with pytest.raises(ValueError):
        impl.saturation_pressure(fl, t["T"][-1] + 1.0)                       # no extrapolation past T_c


@pytest.mark.parametrize("fl", FLUIDS)
def test_log_interpolation_of_psat_beats_linear(fl):
    _, e_lin = impl.leave_one_out_psat(fl, "linear")
    _, e_log = impl.leave_one_out_psat(fl, "log")
    assert np.max(np.abs(e_log)) < 0.2 * np.max(np.abs(e_lin))
    assert np.median(np.abs(e_log)) < 1e-3
    assert np.all(e_lin > -1e-3)                     # p_sat(T) is convex: the chord lies above it


@pytest.mark.parametrize("fl", FLUIDS)
def test_saturation_temperature_inverts_saturation_pressure(fl):
    t = impl.sat_table(fl)
    for T in np.linspace(t["T"][0] + 0.3, t["T"][-1] - 0.3, 7):
        assert impl.saturation_temperature(fl, impl.saturation_pressure(fl, T)) == pytest.approx(T, abs=1e-7)
    assert impl.saturation_temperature("water", 101.42) == pytest.approx(100.0, abs=0.01)


# ---------------------------------------------------------------- phase, quality, mixtures

def test_phase_identification_from_t_and_p():
    ps = impl.saturation_pressure("water", 120.0)
    assert impl.phase_from_Tp("water", 120.0, 2 * ps) == "compressed liquid"
    assert impl.phase_from_Tp("water", 120.0, 0.5 * ps) == "superheated vapor"
    assert impl.phase_from_Tp("water", 120.0, ps) == "saturated"
    assert impl.phase_from_Tp("water", 400.0, 30000.0) == "supercritical"
    # ammonia at room temperature and 1 atm is a vapour; at 1 MPa it is a liquid
    assert impl.phase_from_Tp("ammonia", 20.0, 101.325) == "superheated vapor"
    assert impl.phase_from_Tp("ammonia", 20.0, 1000.0) == "compressed liquid"


@pytest.mark.parametrize("fl, T", [("water", 150.0), ("water", 62.5), ("ammonia", -12.0), ("ammonia", 47.0)])
def test_quality_and_mixture_properties_are_consistent(fl, T):
    st = impl.saturation_state(fl, T)
    for x in (0.0, 0.3, 1.0):
        v = impl.mixture_property(fl, T, x, "v")
        assert impl.phase_from_Tv(fl, T, v)[0] == "saturated mixture"
        assert impl.phase_from_Tv(fl, T, v)[1] == pytest.approx(x, abs=1e-12)
        for prop in "uhs":
            assert impl.quality(fl, T, prop, impl.mixture_property(fl, T, x, prop)) == pytest.approx(x, abs=1e-12)
    assert impl.phase_from_Tv(fl, T, 0.5 * st["vf"])[0] == "compressed liquid"
    assert impl.phase_from_Tv(fl, T, 2.0 * st["vg"])[0] == "superheated vapor"
    with pytest.raises(ValueError):
        impl.mixture_property(fl, T, 1.2, "h")


# ---------------------------------------------------------------- superheated vapour and (T, v) states

@pytest.mark.parametrize("fl", FLUIDS)
def test_superheated_lookup_reproduces_nodes_and_brackets_between_pressures(fl):
    sup = impl.FLUIDS[fl]["sup"]
    ps = sorted(sup)
    for p in ps:
        row = sup[p][len(sup[p]) // 2]
        st = impl.superheated_state(fl, p, row[0])
        assert [st["v"], st["u"], st["h"], st["s"]] == pytest.approx(list(row[1:]), rel=1e-12)
    # between two pressures, v, h and s lie between their values on the neighbouring isobars
    p0, p1 = ps[1], ps[2]
    T = max(sup[p0][0, 0], sup[p1][0, 0]) + 10.0
    a, b, m = impl.superheated_state(fl, p0, T), impl.superheated_state(fl, p1, T), impl.superheated_state(fl, np.sqrt(p0 * p1), T)
    for k in ("v", "h", "s"):
        assert min(a[k], b[k]) < m[k] < max(a[k], b[k])
    with pytest.raises(ValueError):
        impl.superheated_state(fl, 1e6, T)


def test_state_from_t_and_v_round_trips_in_both_regions():
    sh = impl.superheated_state("water", 1000, 350.0)
    st = impl.state_from_Tv("water", 350.0, sh["v"])
    assert st["phase"] == "superheated vapor" and st["p"] == pytest.approx(1000.0, rel=1e-8)
    assert st["h"] == pytest.approx(sh["h"], rel=1e-10)
    sc = impl.superheated_state("water", 2500.0, 600.0)                     # above T_c, between two isobars
    st = impl.state_from_Tv("water", 600.0, sc["v"])
    assert st["phase"] == "supercritical" and st["p"] == pytest.approx(2500.0, rel=1e-8)
    v = impl.mixture_property("ammonia", 10.0, 0.25, "v")
    st = impl.state_from_Tv("ammonia", 10.0, v)
    assert st["phase"] == "saturated mixture" and st["x"] == pytest.approx(0.25)
    assert st["h"] == pytest.approx(impl.mixture_property("ammonia", 10.0, 0.25, "h"))


# ---------------------------------------------------------------- independent checks of the table

@pytest.mark.parametrize("fl", FLUIDS)
def test_clausius_clapeyron_holds_on_the_table(fl):
    T, hcc, htab = impl.clausius_clapeyron_hfg(fl)
    err = hcc / htab - 1
    # where p_sat is printed to at least four significant figures (p >= 0.01 MPa) and away from the critical point
    t = impl.sat_table(fl)
    p = t["p"][1:-1]
    good = (p >= 10.0) & (T < impl.critical_point(fl)[0] - 30.0)
    assert good.sum() > 20
    assert np.max(np.abs(err[good])) < 2e-3
    assert np.median(np.abs(err)) < 5e-4


@pytest.mark.parametrize("fl", FLUIDS)
def test_h_equals_u_plus_pv_to_table_round_off(fl):
    r, b = impl.enthalpy_identity_residual(fl)
    assert len(r) > 100
    assert np.all(np.abs(r) <= b)
    assert np.max(b) < 1.1                           # and the round-off bound itself is about 1 kJ/kg at most


@pytest.mark.parametrize("fl", FLUIDS)
def test_compressibility_factor_tends_to_one_at_low_pressure(fl):
    sup = impl.FLUIDS[fl]["sup"]
    ps = sorted(sup)
    T = max(sup[p][0, 0] for p in ps) + 20.0                     # superheated on every tabulated isobar
    Z = np.array([impl.compressibility_factor(fl, p, T, impl.superheated_state(fl, p, T)["v"]) for p in ps])
    assert np.all(Z < 1.0) and np.all(np.diff(Z) < 0)             # attraction lowers Z, more at higher p
    assert abs(1 - Z[0]) < 0.01
    # 1 - Z is roughly proportional to p at low pressure (second virial coefficient)
    r = (1 - Z[1]) / (1 - Z[0])
    assert 0.5 * ps[1] / ps[0] < r < 1.5 * ps[1] / ps[0]


def test_superheated_steam_matches_the_nasa_ideal_gas_fit_at_low_pressure():
    for T1, T2 in ((100.0, 300.0), (200.0, 700.0), (300.0, 1000.0)):
        tab, ig = impl.superheated_dh_vs_ideal_gas(10, T1, T2)
        assert tab == pytest.approx(ig, rel=5e-3)
    tab, ig = impl.superheated_dh_vs_ideal_gas(4000, 300.0, 600.0)
    assert tab / ig - 1 > 0.05                                    # at 4 MPa the vapour is far from ideal


# ---------------------------------------------------------------- the rigid tank

def test_rigid_tank_heating_first_law_and_fill_limit():
    V, Tf, Tm = 1.0, 15.0, 40.0
    st = impl.saturation_state("ammonia", Tf)
    m = 0.9 * V / st["vf"] + 0.1 * V / st["vg"]
    r = impl.rigid_tank_heating("ammonia", V, m, Tf, Tm)
    assert r["liquid_volume_fraction1"] == pytest.approx(0.9, rel=1e-10)
    assert r["p2"] == pytest.approx(impl.saturation_pressure("ammonia", Tm))
    assert r["liquid_volume_fraction2"] > r["liquid_volume_fraction1"]          # the liquid swells
    u1 = impl.mixture_property("ammonia", Tf, r["x1"], "u")
    u2 = impl.mixture_property("ammonia", Tm, r["x2"], "u")
    assert r["Q"] == pytest.approx(m * (u2 - u1), rel=1e-10) and r["Q"] > 0
    # at the fill limit the tank is exactly liquid-full at T_max
    f = impl.max_liquid_fill("ammonia", Tf, Tm)
    assert 0.9 < f < 1.0
    m_lim = f * V / st["vf"] + (1 - f) * V / st["vg"]
    assert V / m_lim == pytest.approx(impl.saturation_state("ammonia", Tm)["vf"], rel=1e-10)
    assert impl.max_liquid_fill("ammonia", Tf, 50.0) < f                     # a hotter limit needs more ullage


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "property_tables.png").exists()
