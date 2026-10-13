"""Public tests for P00.08 Reading Data. They assert properties of the parser, not magic numbers.

All data are synthetic and generated in code; files are written only under ``tmp_path``.
The HDF5 tests are skipped with a clear reason when ``h5py`` is not installed.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-00/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

SMALL_LOG = """# a comment line
# another: units are in the header

time [s],altitude [ft],airspeed [kt],pitch [deg]
0.0,1000,100,2.0
0.5,1010,,2.1
1.0,nan,102,2.2
1.5,1030,103,---
2.0,1040,104,2.4
"""


# ---------------------------------------------------------------- header and units

def test_parse_header_extracts_names_and_units():
    pairs = impl.parse_header("time [s], Altitude [ ft ],airspeed [m/s],pitch")
    assert pairs == [("time", "s"), ("altitude", "ft"), ("airspeed", "m/s"), ("pitch", "")]


def test_parse_header_rejects_empty_name():
    with pytest.raises(ValueError):
        impl.parse_header("time [s],[m],pitch")


@pytest.mark.parametrize("unit, factor, target", [
    ("ft", 0.3048, "m"), ("kt", 1852.0 / 3600.0, "m/s"), ("km/h", 1.0 / 3.6, "m/s"),
    ("ms", 1e-3, "s"), ("rad", 180.0 / np.pi, "deg"), ("m", 1.0, "m"),
])
def test_convert_to_si_known_units(unit, factor, target):
    values = np.array([0.0, 1.0, 2.5])
    out, u = impl.convert_to_si(values, unit)
    assert u == target
    npt.assert_allclose(out, values * factor, rtol=1e-15)


def test_convert_to_si_passes_through_empty_unit_and_rejects_unknown():
    out, u = impl.convert_to_si(np.array([3.0]), "")
    assert u == "" and out[0] == 3.0
    with pytest.raises(ValueError):
        impl.convert_to_si(np.array([1.0]), "furlongs")


def test_units_round_trip_feet_to_metres():
    """The same flight written in ft/kt and in m/(m/s) must parse to the same SI numbers."""
    si_text = impl.generate_flight_log(n=50, seed=3, missing_fraction=0.0)
    imp_text = impl.generate_flight_log(n=50, seed=3, missing_fraction=0.0,
                                        altitude_unit="ft", airspeed_unit="kt")
    a = impl.parse_flight_log(si_text)
    b = impl.parse_flight_log(imp_text)
    assert a.units == b.units == {"time": "s", "altitude": "m", "airspeed": "m/s", "pitch": "deg"}
    assert b.source_units["altitude"] == "ft" and b.source_units["airspeed"] == "kt"
    # The files are written to 3 decimals, so agreement is limited by that rounding.
    npt.assert_allclose(a.altitude, b.altitude, atol=0.3048 * 1e-3 + 1e-3)
    npt.assert_allclose(a.airspeed, b.airspeed, atol=0.515 * 1e-3 + 1e-3)
    npt.assert_allclose(a.pitch, b.pitch, atol=1e-9)


# ---------------------------------------------------------------- parsing

def test_parse_small_log_handles_comments_missing_values_and_units():
    log = impl.parse_flight_log(SMALL_LOG)
    assert len(log) == 5
    npt.assert_allclose(log.time, [0.0, 0.5, 1.0, 1.5, 2.0])
    assert log.altitude[0] == pytest.approx(1000 * 0.3048)
    assert np.isnan(log.altitude[2]) and np.isnan(log.airspeed[1]) and np.isnan(log.pitch[3])
    assert np.isnan(log.altitude).sum() == 1 and np.isnan(log.airspeed).sum() == 1
    assert log.airspeed[0] == pytest.approx(100 * 1852.0 / 3600.0)
    assert log.pitch[0] == pytest.approx(2.0)


def test_parse_accepts_column_reordering_and_extra_columns():
    text = "pitch [deg],time [s],battery [V],airspeed [m/s],altitude [m]\n1,0,12.1,30,100\n2,1,12.0,31,110\n"
    log = impl.parse_flight_log(text)
    npt.assert_allclose(log.time, [0.0, 1.0])
    npt.assert_allclose(log.altitude, [100.0, 110.0])
    npt.assert_allclose(log.pitch, [1.0, 2.0])


@pytest.mark.parametrize("text, why", [
    ("time [s],altitude [m],airspeed [m/s],pitch [deg]\n0,1,2\n", "ragged row"),
    ("time [s],altitude [m],airspeed [m/s]\n0,1,2\n", "missing channel"),
    ("time [s],altitude [m],airspeed [m/s],pitch [deg]\n0,1,2,3\n0,1,2,3\n", "repeated time stamp"),
    ("time [s],altitude [m],airspeed [m/s],pitch [deg]\n1,1,2,3\n0,1,2,3\n", "time goes backwards"),
    ("time [s],altitude [m],airspeed [m/s],pitch [deg]\n0,1,2,3\n,1,2,3\n", "missing time stamp"),
    ("time [s],altitude [parsecs],airspeed [m/s],pitch [deg]\n0,1,2,3\n", "unknown unit"),
    ("time [s],altitude [m],airspeed [m/s],pitch [deg]\n0,abc,2,3\n", "non-numeric cell"),
    ("# only comments\n", "no header"),
    ("time [s],altitude [m],airspeed [m/s],pitch [deg]\n", "no data rows"),
])
def test_parse_rejects_malformed_logs(text, why):
    with pytest.raises(ValueError):
        impl.parse_flight_log(text)


def test_read_flight_log_from_disk(tmp_path):
    path = tmp_path / "log.csv"
    path.write_text(impl.generate_flight_log(n=120, seed=1), encoding="utf-8")
    log = impl.read_flight_log(path)
    assert len(log) == 120
    assert np.all(np.diff(log.time) > 0.0)
    assert np.nanmax(log.altitude) == pytest.approx(300.0, abs=10.0)  # noise sigma is 1.5 m


def test_generated_log_is_deterministic_and_has_requested_gaps():
    a = impl.generate_flight_log(n=200, seed=7, missing_fraction=0.05)
    b = impl.generate_flight_log(n=200, seed=7, missing_fraction=0.05)
    assert a == b
    log = impl.parse_flight_log(a)
    assert impl.summarize(log)["missing_fraction"] == pytest.approx(0.05, abs=0.5 / 600)


# ---------------------------------------------------------------- gap filling

def test_fill_missing_recovers_a_linear_signal_exactly():
    t = np.array([0.0, 1.0, 2.5, 3.0, 4.0, 6.0])
    y_true = 2.0 * t + 1.0
    y = y_true.copy()
    y[[1, 2, 4]] = np.nan
    filled = impl.fill_missing(t, y)
    npt.assert_allclose(filled, y_true, rtol=1e-14)
    assert np.isnan(y).sum() == 3  # input untouched


def test_fill_missing_holds_edges_and_rejects_all_nan():
    t = np.arange(5.0)
    y = np.array([np.nan, 1.0, np.nan, 3.0, np.nan])
    filled = impl.fill_missing(t, y)
    npt.assert_allclose(filled, [1.0, 1.0, 2.0, 3.0, 3.0])
    assert not np.any(np.isnan(filled))
    with pytest.raises(ValueError):
        impl.fill_missing(t, np.full(5, np.nan))


# ---------------------------------------------------------------- resampling

def test_resample_uniform_grid_properties():
    rng = np.random.default_rng(0)
    t = np.sort(rng.uniform(0.0, 10.0, 80))
    t[0], t[-1] = 0.0, 10.0
    y = np.sin(t)
    t_u, y_u = impl.resample_uniform(t, y, dt=0.25)
    assert t_u[0] == t[0]
    assert t_u[-1] <= t[-1] + 1e-12                       # never extrapolates
    npt.assert_allclose(np.diff(t_u), 0.25, rtol=1e-12)   # uniform spacing
    assert len(t_u) == 41
    assert y_u.shape == t_u.shape


def test_resample_uniform_is_exact_for_linear_and_second_order_otherwise():
    t = np.array([0.0, 0.3, 0.7, 1.2, 2.0, 2.4, 3.0])
    t_u, y_u = impl.resample_uniform(t, 3.0 * t - 1.0, dt=0.1)
    npt.assert_allclose(y_u, 3.0 * t_u - 1.0, rtol=1e-13, atol=1e-13)
    # Linear interpolation error on a smooth signal scales with the SOURCE spacing squared.
    errs = []
    for n in (11, 21, 41):
        ts = np.linspace(0.0, 2.0, n)
        tu, yu = impl.resample_uniform(ts, ts ** 2, dt=0.01)
        errs.append(np.max(np.abs(yu - tu ** 2)))
    rates = np.log2(np.array(errs[:-1]) / np.array(errs[1:]))
    assert np.all(rates > 1.8), f"expected second-order error decay, got {rates}"


def test_resample_uniform_rejects_bad_inputs():
    t = np.arange(5.0)
    with pytest.raises(ValueError):
        impl.resample_uniform(t, t, dt=0.0)
    with pytest.raises(ValueError):
        impl.resample_uniform(np.array([0.0, 2.0, 1.0]), np.zeros(3), dt=0.5)
    with pytest.raises(ValueError):
        impl.resample_uniform(t, np.array([0.0, np.nan, 2.0, 3.0, 4.0]), dt=0.5)


def test_resample_log_preserves_physics_of_the_profile():
    log = impl.parse_flight_log(impl.generate_flight_log(n=300, seed=2, missing_fraction=0.05))
    uniform = impl.resample_log(log, dt=0.5)
    assert not any(np.isnan(getattr(uniform, c)).any() for c in ("altitude", "airspeed", "pitch"))
    npt.assert_allclose(np.diff(uniform.time), 0.5, rtol=1e-12)
    # Peak altitude and mean airspeed are not invented by the resampling.
    assert np.max(uniform.altitude) <= np.nanmax(log.altitude) + 1e-9
    assert np.mean(uniform.airspeed) == pytest.approx(np.nanmean(log.airspeed), rel=0.02)
    assert uniform.units == log.units


def test_summarize_keys_and_values():
    log = impl.parse_flight_log(SMALL_LOG)
    s = impl.summarize(log)
    assert set(s) >= {"duration_s", "max_altitude_m", "max_airspeed_m_s", "mean_dt_s", "missing_fraction"}
    assert s["duration_s"] == pytest.approx(2.0)
    assert s["mean_dt_s"] == pytest.approx(0.5)
    assert s["max_altitude_m"] == pytest.approx(1040 * 0.3048)
    assert s["missing_fraction"] == pytest.approx(3 / 15)


# ---------------------------------------------------------------- HDF5 (optional dependency)

def test_hdf5_round_trip(tmp_path):
    pytest.importorskip("h5py", reason="h5py not installed; HDF5 part of the lesson is optional")
    log = impl.resample_log(impl.parse_flight_log(impl.generate_flight_log(n=100, seed=5)), dt=0.5)
    path = tmp_path / "log.h5"
    assert impl.write_hdf5(path, log) is True
    back = impl.read_hdf5(path)
    for name in ("time", "altitude", "airspeed", "pitch"):
        npt.assert_array_equal(getattr(back, name), getattr(log, name))
    assert back.units == log.units


def test_write_hdf5_without_h5py_is_a_clear_skip(tmp_path, capsys, monkeypatch):
    """When h5py is absent the writer must say so and return False, writing nothing."""
    monkeypatch.setattr(impl, "HAVE_H5PY", False)
    log = impl.parse_flight_log(SMALL_LOG)
    path = tmp_path / "log.h5"
    assert impl.write_hdf5(path, log) is False
    assert not path.exists()
    assert "h5py" in capsys.readouterr().out.lower()
    with pytest.raises(ImportError):
        impl.read_hdf5(path)
