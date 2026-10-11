"""Public tests for Pxx.yy. They assert physics, not magic numbers.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-xx/parameterized-cases/.
"""
import numpy as np
import pytest

from conftest import load_impl

pytestmark = pytest.mark.skip(reason="template lesson; remove this line in a real lesson")
impl = None  # replaced by load_impl(__file__) in a real lesson


def test_matches_analytic_case():
    xs = np.linspace(0.0, 1.0, 201)
    assert np.allclose(impl.model(xs, param=1.0), impl.analytic_reference(xs, param=1.0), rtol=1e-6, atol=1e-9)


def test_limiting_case():
    xs = np.linspace(0.0, 1.0, 11)
    assert np.allclose(impl.model(xs, param=0.0), 0.0, atol=1e-12)


def test_convergence_rate():
    errs = impl.convergence_study([10, 20, 40, 80])
    rates = np.log2(errs[:-1] / errs[1:])
    assert np.all(rates > 1.8), f"expected second-order convergence, got rates {rates}"
