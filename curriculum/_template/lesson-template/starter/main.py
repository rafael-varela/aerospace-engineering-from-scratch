"""Pxx.yy <Lesson Title>: starter code.

Complete every function marked NotImplementedError. When finished, running this file
produces the lesson's plots and prints the verification table.

Run:  python main.py
Test: pytest tests/
"""
from __future__ import annotations

import numpy as np


def model(x: np.ndarray, *, param: float) -> np.ndarray:
    """Return the model prediction at points ``x``.

    Parameters
    ----------
    x : array of independent variable values, SI units.
    param : the model parameter, SI units.
    """
    raise NotImplementedError


def analytic_reference(x: np.ndarray, *, param: float) -> np.ndarray:
    """Closed-form solution for the verification case (from part 2 of the lesson)."""
    raise NotImplementedError


def convergence_study(ns: list[int]) -> np.ndarray:
    """Return the error against the analytic reference for each resolution in ``ns``."""
    raise NotImplementedError


def main() -> None:
    xs = np.linspace(0.0, 1.0, 101)
    y = model(xs, param=1.0)
    y_ref = analytic_reference(xs, param=1.0)
    print(f"max abs error vs analytic: {np.max(np.abs(y - y_ref)):.3e}")
    errs = convergence_study([10, 20, 40, 80])
    rates = np.log2(errs[:-1] / errs[1:])
    print("observed convergence rates:", np.round(rates, 2))


if __name__ == "__main__":
    main()
