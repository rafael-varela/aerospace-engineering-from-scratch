"""Public tests for Pxx.yy. They assert physics, not magic numbers.

Hidden validation cases for assessment live in assessments/phase-xx/parameterized-cases/.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "starter"))
import main  # noqa: E402

pytestmark = pytest.mark.skip(reason="template lesson; remove this line in a real lesson")


def test_matches_analytic_case():
    xs = np.linspace(0.0, 1.0, 201)
    assert np.allclose(main.model(xs, param=1.0), main.analytic_reference(xs, param=1.0), rtol=1e-6, atol=1e-9)


def test_limiting_case():
    # e.g. the response vanishes when the forcing parameter is zero
    xs = np.linspace(0.0, 1.0, 11)
    assert np.allclose(main.model(xs, param=0.0), 0.0, atol=1e-12)


def test_convergence_rate():
    errs = main.convergence_study([10, 20, 40, 80])
    rates = np.log2(errs[:-1] / errs[1:])
    assert np.all(rates > 1.8), f"expected second-order convergence, got rates {rates}"
