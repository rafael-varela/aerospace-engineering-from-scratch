# Parameterized validation cases

Alternate inputs and expected outputs for the phase's computational problems, generated from a seed so answers cannot be copied. Each case is a YAML or JSON file: `{"seed": ..., "inputs": {...}, "expected": {...}, "tolerance": ...}`. A small runner compares the learner's code against them.

Cases here are the "unseen" tests for assessment; the public tests in each lesson's `tests/` are for learning.
