"""P00.08 Reading Data: CSV, HDF5, Telemetry Logs: starter code.

A parser for a small flight log. The log is a CSV file with comment lines, a header that
carries units in square brackets (``altitude [ft]``), irregular time stamps and missing
values. Your parser converts every column to SI, fills gaps, resamples to a uniform time
base, and (when ``h5py`` is installed) round-trips the result through an HDF5 file.

The log is synthetic: ``generate_flight_log`` (given, complete) writes it from a simple
climb-cruise-descent profile with noise, so nothing here is a measurement of a real aircraft.

Complete every function marked NotImplementedError.

Run:  python main.py          (writes flight_log.csv, flight_log.png, optionally flight_log.h5)
Test: pytest curriculum/phase-00/lesson-08-reading-data
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
    Raise ``ValueError`` for an empty field name.
    """
    raise NotImplementedError


def convert_to_si(values: np.ndarray, unit: str) -> tuple[np.ndarray, str]:
    """Convert ``values`` from ``unit`` to the lesson-standard unit; return ``(values, unit)``.

    Known units are the keys of ``UNIT_FACTORS``. An empty unit string is passed through
    unchanged. Any other unit raises ``ValueError``: silently assuming a unit is how the
    Mars Climate Orbiter was lost.
    """
    raise NotImplementedError


def parse_flight_log(text: str) -> FlightLog:
    """Parse the text of a flight-log CSV into a ``FlightLog`` in SI units.

    Rules: lines starting with ``#`` and blank lines are ignored; the first remaining line
    is the header; a missing value (any token in ``MISSING_TOKENS``) becomes ``NaN``;
    every row must have as many fields as the header (``ValueError`` otherwise); the four
    channels in ``CHANNELS`` must all be present (any order, extra columns are ignored);
    time must have no NaN and be strictly increasing after conversion to seconds.
    Record the header units in ``source_units`` and the converted units in ``units``.
    """
    raise NotImplementedError


def read_flight_log(path: str | Path) -> FlightLog:
    """Read a flight-log CSV from disk (see ``parse_flight_log`` for the format)."""
    raise NotImplementedError


def fill_missing(t: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Return a copy of ``y`` with NaNs replaced by linear interpolation in ``t``.

    Interior gaps are interpolated between the nearest valid neighbours; gaps at either
    end are filled by holding the nearest valid value (``np.interp`` does this). Raise
    ``ValueError`` if ``y`` has no valid value at all.
    """
    raise NotImplementedError


def resample_uniform(t: np.ndarray, y: np.ndarray, dt: float) -> tuple[np.ndarray, np.ndarray]:
    """Resample ``y(t)`` by linear interpolation onto a uniform grid of spacing ``dt``.

    The grid starts at ``t[0]`` and ends at or before ``t[-1]`` (never extrapolates).
    ``t`` must be strictly increasing and ``y`` must have no NaN (call ``fill_missing``
    first); ``dt`` must be positive. Raise ``ValueError`` otherwise.
    Returns ``(t_uniform, y_uniform)``.
    """
    raise NotImplementedError


def resample_log(log: FlightLog, dt: float) -> FlightLog:
    """Fill gaps in every channel and resample the whole log to a uniform time base."""
    raise NotImplementedError


def summarize(log: FlightLog) -> dict[str, float]:
    """Return headline numbers with these keys: ``duration_s``, ``max_altitude_m``,
    ``max_airspeed_m_s``, ``mean_dt_s`` and ``missing_fraction`` (missing cells over the
    three data channels divided by ``3 * len(log)``). Ignore NaN when taking maxima."""
    raise NotImplementedError


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

    (Given, complete.) The time stamps are jittered by up to ``time_jitter`` of
    ``dt_nominal`` (telemetry is rarely uniform), about ``missing_fraction`` of the
    altitude/airspeed/pitch cells are blank or ``nan``, and the altitude and airspeed
    columns are written in the requested units with the unit in the header.
    Deterministic for a given ``seed``.
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
    """Write the log to an HDF5 file with one dataset per channel and a ``units`` attribute
    on each dataset. Return ``True`` on success. If ``h5py`` is not importable, print a clear
    message and return ``False`` without writing anything."""
    raise NotImplementedError


def read_hdf5(path: str | Path) -> FlightLog:
    """Read a log written by ``write_hdf5``. Raise ``ImportError`` if h5py is unavailable."""
    raise NotImplementedError


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
