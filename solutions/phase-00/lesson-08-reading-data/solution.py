"""P00.08 Reading Data: CSV, HDF5, Telemetry Logs: reference solution.

A parser for a small flight log. The log is a CSV file with comment lines, a header that
carries units in square brackets (``altitude [ft]``), irregular time stamps and missing
values. The parser converts every column to SI, fills gaps, resamples to a uniform time
base, and (when ``h5py`` is installed) round-trips the result through an HDF5 file.

The log is synthetic: ``generate_flight_log`` writes it from a simple climb-cruise-descent
profile with noise, so nothing here is a measurement of a real aircraft.

Run:  python solution.py          (writes flight_log.csv, flight_log.png, optionally flight_log.h5)
Test: AEFS_IMPL=solution pytest curriculum/phase-00/lesson-08-reading-data
"""
from __future__ import annotations

import io
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

try:  # HDF5 is optional: the lesson works without it and says so.
    import h5py  # type: ignore

    HAVE_H5PY = True
except ImportError:  # pragma: no cover - depends on the environment
    h5py = None
    HAVE_H5PY = False

CHANNELS = ("time", "altitude", "airspeed", "pitch")

# Multiply a value in the key unit by the factor to get the SI (or lesson-standard) unit.
# Angles are kept in degrees because that is how flight data is reported; "deg" is the
# lesson-standard angle unit, not an SI unit.
UNIT_FACTORS: dict[str, tuple[str, float]] = {
    "s": ("s", 1.0),
    "ms": ("s", 1e-3),
    "min": ("s", 60.0),
    "m": ("m", 1.0),
    "km": ("m", 1000.0),
    "ft": ("m", 0.3048),
    "m/s": ("m/s", 1.0),
    "km/h": ("m/s", 1000.0 / 3600.0),
    "kt": ("m/s", 1852.0 / 3600.0),
    "kn": ("m/s", 1852.0 / 3600.0),
    "deg": ("deg", 1.0),
    "rad": ("deg", 180.0 / np.pi),
}

MISSING_TOKENS = {"", "nan", "NaN", "NAN", "null", "NULL", "None", "-", "--", "---", "N/A", "n/a"}


@dataclass
class FlightLog:
    """Four channels on a common time base. All arrays have the same length."""

    time: np.ndarray            # s
    altitude: np.ndarray        # m
    airspeed: np.ndarray        # m/s
    pitch: np.ndarray           # deg
    units: dict[str, str] = field(default_factory=lambda: {"time": "s", "altitude": "m",
                                                           "airspeed": "m/s", "pitch": "deg"})
    source_units: dict[str, str] = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.time)


def parse_header(line: str) -> list[tuple[str, str]]:
    """Split a CSV header line into ``(name, unit)`` pairs.

    ``"time [s],altitude [ft],pitch"`` -> ``[("time", "s"), ("altitude", "ft"), ("pitch", "")]``.
    Names are lower-cased and stripped; a column without brackets gets unit ``""``.
    """
    pairs: list[tuple[str, str]] = []
    for raw in line.strip().split(","):
        m = re.fullmatch(r"\s*([^\[\]]+?)\s*(?:\[\s*([^\[\]]*?)\s*\])?\s*", raw)
        if m is None or not m.group(1).strip():
            raise ValueError(f"malformed header field {raw!r}")
        pairs.append((m.group(1).strip().lower(), (m.group(2) or "").strip()))
    return pairs


def convert_to_si(values: np.ndarray, unit: str) -> tuple[np.ndarray, str]:
    """Convert ``values`` from ``unit`` to the lesson-standard unit; return ``(values, unit)``.

    Known units are the keys of ``UNIT_FACTORS``. An empty unit string is passed through
    unchanged. Any other unit raises ``ValueError``: silently assuming a unit is how the
    Mars Climate Orbiter was lost.
    """
    if unit == "":
        return np.asarray(values, dtype=float), ""
    if unit not in UNIT_FACTORS:
        raise ValueError(f"unknown unit {unit!r}; add it to UNIT_FACTORS if it is legitimate")
    target, factor = UNIT_FACTORS[unit]
    return np.asarray(values, dtype=float) * factor, target


def _parse_value(token: str) -> float:
    token = token.strip()
    if token in MISSING_TOKENS:
        return np.nan
    return float(token)


def parse_flight_log(text: str) -> FlightLog:
    """Parse the text of a flight-log CSV into a ``FlightLog`` in SI units.

    Rules: lines starting with ``#`` and blank lines are ignored; the first remaining line
    is the header; a missing value (empty cell, ``nan``, ``---`` ...) becomes ``NaN``;
    every row must have as many fields as the header; the four channels in ``CHANNELS``
    must all be present (any order, extra columns are ignored); time must be strictly
    increasing after conversion to seconds.
    """
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.lstrip().startswith("#")]
    if not lines:
        raise ValueError("empty log: no header line found")
    header = parse_header(lines[0])
    names = [n for n, _ in header]
    missing = [c for c in CHANNELS if c not in names]
    if missing:
        raise ValueError(f"log is missing channels {missing}; header has {names}")

    rows: list[list[float]] = []
    for lineno, ln in enumerate(lines[1:], start=2):
        fields = ln.split(",")
        if len(fields) != len(header):
            raise ValueError(f"line {lineno}: expected {len(header)} fields, got {len(fields)}")
        try:
            rows.append([_parse_value(f) for f in fields])
        except ValueError as exc:
            raise ValueError(f"line {lineno}: {exc}") from None
    if not rows:
        raise ValueError("log has a header but no data rows")
    data = np.array(rows, dtype=float)

    columns: dict[str, np.ndarray] = {}
    units: dict[str, str] = {}
    source_units: dict[str, str] = {}
    for name in CHANNELS:
        j = names.index(name)
        unit_in = header[j][1]
        values, unit_out = convert_to_si(data[:, j], unit_in)
        columns[name] = values
        units[name] = unit_out
        source_units[name] = unit_in

    t = columns["time"]
    if np.any(np.isnan(t)):
        raise ValueError("time column has missing values; a sample without a time stamp is unusable")
    if np.any(np.diff(t) <= 0.0):
        bad = int(np.argmax(np.diff(t) <= 0.0)) + 2
        raise ValueError(f"time must be strictly increasing; first violation near data row {bad}")
    return FlightLog(t, columns["altitude"], columns["airspeed"], columns["pitch"], units, source_units)


def read_flight_log(path: str | Path) -> FlightLog:
    """Read a flight-log CSV from disk (see ``parse_flight_log`` for the format)."""
    return parse_flight_log(Path(path).read_text(encoding="utf-8"))


def fill_missing(t: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Return a copy of ``y`` with NaNs replaced by linear interpolation in ``t``.

    Interior gaps are interpolated between the nearest valid neighbours; gaps at either
    end are filled by holding the nearest valid value (``np.interp`` does this). Raises
    ``ValueError`` if ``y`` has no valid value at all.
    """
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = ~np.isnan(y)
    if not np.any(ok):
        raise ValueError("cannot fill a channel that has no valid samples")
    if np.all(ok):
        return y.copy()
    return np.interp(t, t[ok], y[ok])


def resample_uniform(t: np.ndarray, y: np.ndarray, dt: float) -> tuple[np.ndarray, np.ndarray]:
    """Resample ``y(t)`` by linear interpolation onto a uniform grid of spacing ``dt``.

    The grid starts at ``t[0]`` and ends at or before ``t[-1]`` (never extrapolates).
    ``t`` must be strictly increasing and ``y`` must have no NaN (call ``fill_missing``
    first). ``dt`` must be positive. Returns ``(t_uniform, y_uniform)``.
    """
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    if dt <= 0.0:
        raise ValueError("dt must be positive")
    if t.ndim != 1 or t.shape != y.shape:
        raise ValueError("t and y must be 1-D arrays of the same length")
    if np.any(np.diff(t) <= 0.0):
        raise ValueError("t must be strictly increasing")
    if np.any(np.isnan(y)):
        raise ValueError("y contains NaN; fill missing values before resampling")
    n = int(np.floor((t[-1] - t[0]) / dt + 1e-9)) + 1
    t_u = t[0] + dt * np.arange(n)
    return t_u, np.interp(t_u, t, y)


def resample_log(log: FlightLog, dt: float) -> FlightLog:
    """Fill gaps in every channel and resample the whole log to a uniform time base."""
    t_u = None
    out: dict[str, np.ndarray] = {}
    for name in CHANNELS[1:]:
        y = fill_missing(log.time, getattr(log, name))
        t_u, out[name] = resample_uniform(log.time, y, dt)
    assert t_u is not None
    return FlightLog(t_u, out["altitude"], out["airspeed"], out["pitch"], dict(log.units), dict(log.source_units))


def summarize(log: FlightLog) -> dict[str, float]:
    """Return headline numbers: duration (s), max altitude (m), max airspeed (m/s),
    mean sample interval (s), and the fraction of missing samples over the three channels."""
    n = len(log)
    missing = sum(int(np.isnan(getattr(log, c)).sum()) for c in CHANNELS[1:])
    return {
        "duration_s": float(log.time[-1] - log.time[0]),
        "max_altitude_m": float(np.nanmax(log.altitude)),
        "max_airspeed_m_s": float(np.nanmax(log.airspeed)),
        "mean_dt_s": float(np.mean(np.diff(log.time))) if n > 1 else float("nan"),
        "missing_fraction": missing / (3 * n),
    }


def generate_flight_log(
    n: int = 600,
    dt_nominal: float = 0.5,
    seed: int = 0,
    missing_fraction: float = 0.03,
    altitude_unit: str = "m",
    airspeed_unit: str = "m/s",
    time_jitter: float = 0.1,
) -> str:
    """Return the text of a synthetic flight log: a climb, a cruise and a descent.

    The time stamps are jittered by up to ``time_jitter`` of ``dt_nominal`` (telemetry is
    rarely uniform), about ``missing_fraction`` of the altitude/airspeed/pitch cells are
    blank or ``nan``, and the altitude and airspeed columns are written in the requested
    units with the unit in the header. Deterministic for a given ``seed``.
    """
    if n < 2:
        raise ValueError("need at least two samples")
    rng = np.random.default_rng(seed)
    t = dt_nominal * np.arange(n) + rng.uniform(-time_jitter, time_jitter, n) * dt_nominal
    t[0] = 0.0
    t = np.maximum.accumulate(t + 1e-6 * np.arange(n))  # keep strictly increasing
    duration = t[-1]
    # Profile: climb to 300 m at 3 m/s for the first third, cruise, descend.
    t1, t2 = duration / 3.0, 2.0 * duration / 3.0
    climb_rate = 300.0 / t1
    alt = np.where(t < t1, climb_rate * t,
                   np.where(t < t2, 300.0, np.maximum(300.0 - climb_rate * (t - t2), 0.0)))
    speed = np.where(t < t1, 28.0, np.where(t < t2, 35.0, 30.0))
    pitch = np.where(t < t1, 8.0, np.where(t < t2, 2.0, -4.0))
    alt = alt + rng.normal(0.0, 1.5, n)
    speed = speed + rng.normal(0.0, 0.6, n)
    pitch = pitch + rng.normal(0.0, 0.4, n)

    alt_out = alt / UNIT_FACTORS[altitude_unit][1]
    spd_out = speed / UNIT_FACTORS[airspeed_unit][1]

    cells = np.array([alt_out, spd_out, pitch], dtype=object)
    k = int(round(missing_fraction * 3 * n))
    if k:
        flat = rng.choice(3 * n, size=k, replace=False)
        for idx in flat:
            cells[idx // n, idx % n] = "" if idx % 2 == 0 else "nan"

    buf = io.StringIO()
    buf.write("# synthetic flight log, generated by P00.08 generate_flight_log\n")
    buf.write(f"# seed={seed} n={n} dt_nominal={dt_nominal} s\n")
    buf.write(f"time [s],altitude [{altitude_unit}],airspeed [{airspeed_unit}],pitch [deg]\n")
    for i in range(n):
        row = [f"{t[i]:.3f}"]
        for c in range(3):
            v = cells[c, i]
            row.append(v if isinstance(v, str) else f"{float(v):.3f}")
        buf.write(",".join(row) + "\n")
    return buf.getvalue()


def write_hdf5(path: str | Path, log: FlightLog) -> bool:
    """Write the log to an HDF5 file with one dataset per channel and ``units`` attributes.

    Returns ``True`` on success. If ``h5py`` is not importable, prints a clear message and
    returns ``False`` without writing anything.
    """
    if not HAVE_H5PY:
        print("h5py is not installed: skipping HDF5 output (pip install h5py to enable it)")
        return False
    with h5py.File(Path(path), "w") as f:
        f.attrs["source"] = "P00.08 FlightLog"
        for name in CHANNELS:
            ds = f.create_dataset(name, data=getattr(log, name))
            ds.attrs["units"] = log.units[name]
    return True


def read_hdf5(path: str | Path) -> FlightLog:
    """Read a log written by ``write_hdf5``. Raises ``ImportError`` if h5py is unavailable."""
    if not HAVE_H5PY:
        raise ImportError("h5py is not installed; cannot read HDF5")
    with h5py.File(Path(path), "r") as f:
        arrays = {name: np.asarray(f[name][...], dtype=float) for name in CHANNELS}
        units = {name: str(f[name].attrs["units"]) for name in CHANNELS}
    return FlightLog(arrays["time"], arrays["altitude"], arrays["airspeed"], arrays["pitch"], units)


def main(out_dir: str | Path = ".") -> None:
    """Generate, write, read, clean and plot a synthetic flight log."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    csv_path = out_dir / "flight_log.csv"
    csv_path.write_text(generate_flight_log(altitude_unit="ft", airspeed_unit="kt"), encoding="utf-8")
    log = read_flight_log(csv_path)
    print(f"read {len(log)} samples from {csv_path}")
    print("source units:", log.source_units, "-> stored as", log.units)
    for key, val in summarize(log).items():
        print(f"  {key:>18s}: {val:10.3f}")
    uniform = resample_log(log, dt=0.5)
    print(f"resampled to {len(uniform)} samples at dt = 0.5 s "
          f"(spacing spread {np.ptp(np.diff(uniform.time)):.2e} s)")

    if write_hdf5(out_dir / "flight_log.h5", uniform):
        back = read_hdf5(out_dir / "flight_log.h5")
        print(f"HDF5 round trip: max |altitude difference| = "
              f"{np.max(np.abs(back.altitude - uniform.altitude)):.2e} m")

    fig, axes = plt.subplots(3, 1, figsize=(7, 7), sharex=True)
    for ax, name, unit in zip(axes, CHANNELS[1:], ("m", "m/s", "deg")):
        ax.plot(log.time, getattr(log, name), ".", ms=3, alpha=0.5, label="raw (gaps as NaN)")
        ax.plot(uniform.time, getattr(uniform, name), "-", lw=1, label="filled, resampled")
        ax.set_ylabel(f"{name} [{unit}]")
        ax.grid(True)
    axes[0].legend()
    axes[-1].set_xlabel("time [s]")
    fig.tight_layout()
    fig.savefig(out_dir / "flight_log.png", dpi=150)
    print(f"saved {out_dir / 'flight_log.png'}")


if __name__ == "__main__":
    main()
