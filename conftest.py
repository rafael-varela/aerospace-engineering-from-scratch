"""Shared pytest helpers for lesson tests.

A lesson's tests import the learner's code from ``starter/main.py``. In CI, and whenever
``AEFS_IMPL=solution`` is set, the same tests run against the reference implementation in
``solutions/phase-XX/lesson-YY-<slug>/solution.py`` instead, which proves the tests are
satisfiable and the solution is correct.

Usage in ``tests/test_main.py``::

    from conftest import load_impl
    impl = load_impl(__file__)

Run against your own code:      pytest curriculum/phase-00/lesson-03-numpy-engineering-calculator
Run against the reference:      AEFS_IMPL=solution pytest curriculum/phase-00/...
"""
from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))  # so `from conftest import load_impl` works from any test file


def load_impl(test_file: str):
    """Import the implementation under test for the lesson that owns ``test_file``."""
    lesson_dir = Path(test_file).resolve().parents[1]          # .../lesson-YY-slug
    rel = lesson_dir.relative_to(ROOT / "curriculum")           # phase-XX/lesson-YY-slug
    if os.environ.get("AEFS_IMPL") == "solution":
        target = ROOT / "solutions" / rel / "solution.py"
    else:
        target = lesson_dir / "starter" / "main.py"
    if not target.exists():
        raise FileNotFoundError(f"implementation not found: {target}")
    name = "aefs_" + "_".join(rel.parts).replace("-", "_")
    spec = importlib.util.spec_from_file_location(name, target)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module
