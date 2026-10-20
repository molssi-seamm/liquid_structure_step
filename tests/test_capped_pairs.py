# -*- coding: utf-8 -*-

"""capped_pairs finds the same pairs and distances as MDAnalysis' capped_distance,
faster, for rectangular boxes (job 5276's RDFs took ~5 s a frame, 2026-10-20)."""

import numpy as np
import pytest
from MDAnalysis.lib.distances import capped_distance

from liquid_structure_step.analysis import capped_pairs


def _as_set(pairs, distances):
    return {(int(i), int(j)): float(d) for (i, j), d in zip(pairs, distances)}


@pytest.mark.parametrize("box", [[20.0, 22.0, 25.0, 90, 90, 90], None])
def test_same_pairs_as_mdanalysis(box):
    rng = np.random.default_rng(5)
    a = rng.random((300, 3)) * 20.0 - 2.0  # some outside the box, as in a trajectory
    b = rng.random((400, 3)) * 20.0
    box = None if box is None else np.array(box, dtype=float)
    p1, d1 = capped_distance(a, b, 7.5, box=box, return_distances=True)
    p2, d2 = capped_pairs(a, b, 7.5, box=box)
    s1, s2 = _as_set(p1, d1), _as_set(p2, d2)
    assert set(s1) == set(s2)
    assert max(abs(s1[k] - s2[k]) for k in s1) < 1.0e-4


def test_same_atoms_in_both_sets():
    """An RDF of one element with itself: every ordered pair, both ways."""
    rng = np.random.default_rng(6)
    x = rng.random((200, 3)) * 15.0
    box = np.array([15.0, 15.0, 15.0, 90, 90, 90])
    p1, d1 = capped_distance(x, x, 6.0, box=box, return_distances=True)
    p2, d2 = capped_pairs(x, x, 6.0, box=box)
    nonzero = lambda p, d: {k for k, v in _as_set(p, d).items() if v > 0}  # noqa
    assert nonzero(p1, d1) == nonzero(p2, d2)


def test_triclinic_and_long_cutoffs_use_mdanalysis():
    rng = np.random.default_rng(7)
    a, b = rng.random((50, 3)) * 10.0, rng.random((60, 3)) * 10.0
    triclinic = np.array([10.0, 10.0, 10.0, 90, 90, 75])
    p1, d1 = capped_distance(a, b, 4.0, box=triclinic, return_distances=True)
    p2, d2 = capped_pairs(a, b, 4.0, box=triclinic)
    assert _as_set(p1, d1) == pytest.approx(_as_set(p2, d2))
    cubic = np.array([10.0, 10.0, 10.0, 90, 90, 90])
    p2 = capped_pairs(a, b, 6.0, box=cubic, return_distances=False)  # > L/2
    assert len(p2) > 0


def test_empty():
    pairs, d = capped_pairs(np.empty((0, 3)), np.zeros((3, 3)), 5.0)
    assert pairs.shape == (0, 2) and len(d) == 0
