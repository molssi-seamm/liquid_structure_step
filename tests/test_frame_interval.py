#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Frames written more often than the analysis needs are skipped."""

from contextlib import contextmanager
from types import SimpleNamespace

import numpy as np

import liquid_structure_step.liquid_structure as module


class _Times:
    """The time axis of the positions, 0.1 ps apart."""

    def __init__(self, n, dt=0.1):
        self.values = np.arange(n) * dt
        self.shape = self.values.shape
        self.attrs = {"unit": "ps"}

    def __getitem__(self, i):
        return self.values[i]


def _fake_md(n):
    frames = [
        SimpleNamespace(positions=np.full((2, 3), float(i)), dimensions=None, volume=0)
        for i in range(n)
    ]
    md = SimpleNamespace(
        first={"particles/all/position/time": _Times(n)},
        universe=lambda configuration: SimpleNamespace(trajectory=frames),
    )

    @contextmanager
    def open_md(path):
        yield md

    return open_md


def _analyzed(monkeypatch, wanted, n=100):
    monkeypatch.setattr(module, "open_md", _fake_md(n))
    step = SimpleNamespace(paths=["run1"], _wanted_dt=wanted)
    frames = list(module.LiquidStructure._frames(step, None, periodic=False))
    return [int(x[0][0, 0]) for x in frames], step._dt, step._run_frames


def test_frames_every_ps_from_frames_every_tenth(monkeypatch):
    kept, dt, runs = _analyzed(monkeypatch, wanted=1.0)
    assert kept == list(range(0, 100, 10))
    assert dt == 1.0 and runs == [10]


def test_frames_every_tenth_kept_for_residence_times(monkeypatch):
    kept, dt, runs = _analyzed(monkeypatch, wanted=0.1)
    assert kept == list(range(100)) and runs == [100]
    assert abs(dt - 0.1) < 1e-12


def test_an_uneven_ratio_stays_within_the_interval(monkeypatch):
    kept, dt, _ = _analyzed(monkeypatch, wanted=0.25)
    assert kept[:3] == [0, 2, 4] and abs(dt - 0.2) < 1e-12
