"""P00.06 SciPy Tour: reference solution.

One synthetic engineering dataset, a model-rocket motor thrust curve, worked through
four SciPy submodules:

* ``scipy.integrate``   total impulse from the thrust curve
* ``scipy.optimize``    fit the thrust model and find the peak
* ``scipy.interpolate`` look up the atmosphere between table rows
* ``scipy.linalg``      resolve a three-cell wind-tunnel balance

The thrust samples are *synthetic*: generated in code from a known model with added
noise, so every answer has an analytic reference. Nothing here is a measurement.

Run:  python solution.py            (prints the report, writes figures/thrust_fit.svg)
Test: AEFS_IMPL=solution pytest curriculum/phase-00/lesson-06-scipy-tour
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy import integrate, interpolate, linalg, optimize  # noqa: E402

# --------------------------------------------------------------------------- thrust model

def thrust_model(t: np.ndarray, a: float, tau: float) -> np.ndarray:
    """Thrust in N at time ``t`` (s): ``F(t) = a * t * exp(-t / tau)``.

    ``a`` has units N/s and ``tau`` seconds. The curve rises linearly, peaks at
    ``t = tau`` with ``F = a * tau / e``, and decays. It is zero for ``t < 0``.
    """
    t = np.asarray(t, dtype=float)
    return np.where(t >= 0.0, a * t * np.exp(-t / tau), 0.0)


def analytic_impulse(a: float, tau: float, t_end: float) -> float:
    """Exact integral of :func:`thrust_model` from 0 to ``t_end`` (N·s).

    By parts: ``∫0^T a t e^{-t/τ} dt = a τ² [1 − e^{−T/τ}(1 + T/τ)]``.
    """
    x = t_end / tau
    return a * tau**2 * (1.0 - np.exp(-x) * (1.0 + x))


def synthesize_thrust_data(
    a: float = 60.0, tau: float = 0.3, t_end: float = 2.0, n: int = 81, noise: float = 0.3, seed: int = 0
) -> tuple[np.ndarray, np.ndarray]:
    """Synthetic thrust samples ``(t, F)``: the model plus Gaussian noise of std ``noise`` N.

    ``noise=0`` gives the exact model, which the tests use for analytic comparisons.
    """
    rng = np.random.default_rng(seed)
    t = np.linspace(0.0, t_end, n)
    f = thrust_model(t, a, tau) + rng.normal(0.0, noise, size=t.shape) * (noise > 0)
    return t, f


# --------------------------------------------------------------------------- integrate

def total_impulse_samples(t: np.ndarray, f: np.ndarray) -> float:
    """Total impulse (N·s) of sampled thrust by the trapezoidal rule."""
    t = np.asarray(t, dtype=float)
    f = np.asarray(f, dtype=float)
    if t.shape != f.shape or t.ndim != 1:
        raise ValueError("t and f must be 1-D arrays of the same length")
    return float(integrate.trapezoid(f, t))


def total_impulse_model(a: float, tau: float, t_end: float) -> float:
    """Total impulse (N·s) of the continuous model by adaptive quadrature (``quad``)."""
    value, _err = integrate.quad(lambda t: thrust_model(t, a, tau), 0.0, t_end)
    return float(value)


def trapezoid_convergence(a: float, tau: float, t_end: float, ns: list[int]) -> np.ndarray:
    """Absolute error of the trapezoidal impulse against the analytic value for each ``n``."""
    exact = analytic_impulse(a, tau, t_end)
    errs = []
    for n in ns:
        t = np.linspace(0.0, t_end, n)
        errs.append(abs(total_impulse_samples(t, thrust_model(t, a, tau)) - exact))
    return np.array(errs)


def motor_class(impulse_ns: float) -> str:
    """Model-rocket motor letter class for a total impulse in N·s.

    Classes double: 'A' is 1.26–2.50 N·s, 'B' 2.51–5.00, 'C' 5.01–10.0, 'D' 10.01–20.0,
    and so on up to 'O' (20 480–40 960). Below 1.26 N·s the class is '1/2A' or smaller;
    this function returns 'micro' for those.
    """
    if impulse_ns <= 0.0:
        raise ValueError("impulse must be positive")
    if impulse_ns <= 1.25:
        return "micro"
    upper = 2.5
    for letter in "ABCDEFGHIJKLMNO":
        if impulse_ns <= upper:
            return letter
        upper *= 2.0
    raise ValueError("impulse exceeds the hobby motor classes")


# --------------------------------------------------------------------------- optimize

def fit_thrust_model(t: np.ndarray, f: np.ndarray) -> tuple[float, float]:
    """Least-squares fit of ``(a, tau)`` in :func:`thrust_model` to samples (``curve_fit``).

    The initial guess comes from the data: ``tau`` ≈ time of the largest sample and
    ``a`` ≈ peak * e / tau.
    """
    t = np.asarray(t, dtype=float)
    f = np.asarray(f, dtype=float)
    i_peak = int(np.argmax(f))
    tau0 = max(t[i_peak], 1e-3)
    a0 = max(f[i_peak], 1e-3) * np.e / tau0
    popt, _pcov = optimize.curve_fit(thrust_model, t, f, p0=(a0, tau0))
    return float(popt[0]), float(popt[1])


def peak_thrust(a: float, tau: float, t_end: float) -> tuple[float, float]:
    """Time and value of the thrust maximum on ``[0, t_end]`` (``minimize_scalar``).

    Analytic answer: ``t = tau``, ``F = a * tau / e`` (when ``tau <= t_end``).
    """
    res = optimize.minimize_scalar(
        lambda t: -thrust_model(t, a, tau), bounds=(0.0, t_end), method="bounded",
        options={"xatol": 1e-10},
    )
    return float(res.x), float(-res.fun)


def burnout_time(a: float, tau: float, threshold: float, t_end: float) -> float:
    """First time after the peak at which thrust falls to ``threshold`` N (``brentq``).

    Raises ``ValueError`` if the thrust never drops to ``threshold`` before ``t_end``
    or if ``threshold`` is above the peak.
    """
    t_peak, f_peak = tau, a * tau / np.e
    if not 0.0 < threshold < f_peak:
        raise ValueError("threshold must lie between 0 and the peak thrust")
    g = lambda t: thrust_model(t, a, tau) - threshold  # noqa: E731
    if g(t_end) > 0.0:
        raise ValueError("thrust has not decayed to the threshold by t_end")
    return float(optimize.brentq(g, t_peak, t_end, xtol=1e-12))


# --------------------------------------------------------------------------- interpolate

#: A coarse ISA table, altitude (m) vs pressure (Pa) and temperature (K), every 2 km.
ISA_TABLE_H = np.arange(0.0, 20_001.0, 2_000.0)
_G0, _R, _T0, _P0, _L = 9.80665, 287.05287, 288.15, 101_325.0, -0.0065


def isa_pressure_exact(h: np.ndarray) -> np.ndarray:
    """Closed-form ISA pressure in Pa for 0–20 km (troposphere + isothermal layer)."""
    h = np.asarray(h, dtype=float)
    h11 = 11_000.0
    t11 = _T0 + _L * h11
    p11 = _P0 * (t11 / _T0) ** (-_G0 / (_L * _R))
    p_tropo = _P0 * ((_T0 + _L * np.minimum(h, h11)) / _T0) ** (-_G0 / (_L * _R))
    p_strato = p11 * np.exp(-_G0 * (h - h11) / (_R * t11))
    return np.where(h <= h11, p_tropo, p_strato)


ISA_TABLE_P = isa_pressure_exact(ISA_TABLE_H)


def interpolate_pressure(h: np.ndarray, kind: str = "log-linear") -> np.ndarray:
    """Pressure at altitude ``h`` (m) from the coarse table.

    ``kind="linear"`` interpolates pressure directly; ``kind="log-linear"``
    interpolates ``ln p`` linearly, which is exact for an isothermal layer and far
    better in the troposphere. ``kind="cubic"`` uses a cubic spline on ``ln p``.
    Raises ``ValueError`` outside the table range (no silent extrapolation).
    """
    h = np.asarray(h, dtype=float)
    if np.any(h < ISA_TABLE_H[0]) or np.any(h > ISA_TABLE_H[-1]):
        raise ValueError("altitude outside the table range 0–20 km; refusing to extrapolate")
    if kind == "linear":
        return interpolate.interp1d(ISA_TABLE_H, ISA_TABLE_P)(h)
    if kind == "log-linear":
        return np.exp(interpolate.interp1d(ISA_TABLE_H, np.log(ISA_TABLE_P))(h))
    if kind == "cubic":
        return np.exp(interpolate.CubicSpline(ISA_TABLE_H, np.log(ISA_TABLE_P))(h))
    raise ValueError(f"unknown kind {kind!r}")


# --------------------------------------------------------------------------- linalg

#: Geometry of a three-cell wind-tunnel balance. The model is held by a sting; a
#: forward cell and an aft cell measure vertical force at x = +0.15 m and x = -0.15 m
#: from the reference point, and a drag cell measures axial force. Readings are
#: ``r = B @ [L, D, M]`` with ``M`` the pitching moment about the reference point.
BALANCE_ARM = 0.15  # m


def balance_matrix(arm: float = BALANCE_ARM) -> np.ndarray:
    """Matrix ``B`` mapping ``[L, D, M]`` to the three cell readings ``[r_fwd, r_aft, r_drag]``.

    Vertical equilibrium: ``r_fwd + r_aft = L``. Moment about the reference point:
    ``arm * r_fwd − arm * r_aft = M``. The drag cell reads ``D`` plus 2 % of lift from
    a slight misalignment (a realistic interaction term).
    """
    return np.array([
        [0.5, 0.0, 0.5 / arm],
        [0.5, 0.0, -0.5 / arm],
        [0.02, 1.0, 0.0],
    ])


def resolve_balance(readings: np.ndarray, arm: float = BALANCE_ARM) -> np.ndarray:
    """Recover ``[L, D, M]`` (N, N, N·m) from the three cell readings by ``linalg.solve``."""
    r = np.asarray(readings, dtype=float)
    if r.shape != (3,):
        raise ValueError("expected three readings [r_fwd, r_aft, r_drag]")
    return linalg.solve(balance_matrix(arm), r)


def balance_condition_number(arm: float = BALANCE_ARM) -> float:
    """2-norm condition number of the balance matrix: how much reading noise is amplified."""
    return float(np.linalg.cond(balance_matrix(arm)))


# --------------------------------------------------------------------------- report

def main(out_dir: str | Path = "figures") -> None:
    """Run the four-submodule tour on the synthetic thrust data and print a report."""
    a, tau, t_end = 60.0, 0.3, 2.0
    t, f = synthesize_thrust_data(a, tau, t_end, noise=0.3, seed=0)

    print("integrate:")
    i_trap = total_impulse_samples(t, f)
    i_quad = total_impulse_model(a, tau, t_end)
    i_exact = analytic_impulse(a, tau, t_end)
    print(f"  trapezoid on noisy samples : {i_trap:8.4f} N*s  -> class {motor_class(i_trap)}")
    print(f"  quad on the model          : {i_quad:8.4f} N*s")
    print(f"  analytic                   : {i_exact:8.4f} N*s")
    ns = [11, 21, 41, 81, 161]
    errs = trapezoid_convergence(a, tau, t_end, ns)
    rates = np.log2(errs[:-1] / errs[1:])
    print(f"  trapezoid error vs n {ns}: {np.array2string(errs, precision=2)}")
    print(f"  observed convergence rates : {np.round(rates, 2)}  (trapezoid: 2)")

    print("optimize:")
    a_fit, tau_fit = fit_thrust_model(t, f)
    print(f"  fitted a = {a_fit:.3f} N/s (true {a}), tau = {tau_fit:.4f} s (true {tau})")
    t_pk, f_pk = peak_thrust(a_fit, tau_fit, t_end)
    print(f"  peak {f_pk:.3f} N at t = {t_pk:.4f} s   (analytic a*tau/e = {a * tau / np.e:.3f} N at {tau} s)")
    print(f"  thrust falls below 1 N at t = {burnout_time(a_fit, tau_fit, 1.0, t_end):.4f} s")

    print("interpolate:")
    h = np.array([1_000.0, 5_500.0, 13_000.0])
    for kind in ("linear", "log-linear", "cubic"):
        p = interpolate_pressure(h, kind)
        rel = (p - isa_pressure_exact(h)) / isa_pressure_exact(h)
        print(f"  {kind:<10s} max |rel error| at {h / 1000} km: {np.max(np.abs(rel)) * 100:.3f} %")

    print("linalg:")
    truth = np.array([120.0, 8.0, 3.0])  # L N, D N, M N*m
    readings = balance_matrix() @ truth
    recovered = resolve_balance(readings)
    print(f"  readings {np.round(readings, 4)} -> [L, D, M] = {np.round(recovered, 6)}")
    print(f"  condition number of the balance: {balance_condition_number():.2f}")

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(t, f, "o", ms=3, color="#1f5fa8", label="synthetic samples")
    tt = np.linspace(0.0, t_end, 400)
    ax.plot(tt, thrust_model(tt, a_fit, tau_fit), color="#2a9d8f", label="curve_fit")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Thrust [N]")
    ax.legend()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / "thrust_fit.svg")
    plt.close(fig)
    print(f"wrote {out / 'thrust_fit.svg'}")


if __name__ == "__main__":
    main()
