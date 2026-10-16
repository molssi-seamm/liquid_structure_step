# -*- coding: utf-8 -*-

"""Tests of the structural analysis (liquid_structure_step.analysis) on systems
with known answers."""

import math

import numpy as np
import pytest

from liquid_structure_step.analysis import (
    AnalysisError,
    PairingAccumulator,
    RDF,
    ShellAccumulator,
    Species,
    classify_ion_pairs,
    coordination_number,
    hill_formula,
    parse_cutoffs,
    parse_pairs,
    peak_and_minimum,
    residence_time,
    second_minimum,
    split_pair,
)


def cubic_lattice(n=6, a=3.0):
    xyz = np.array(
        [[i, j, k] for i in range(n) for j in range(n) for k in range(n)], float
    )
    return xyz * a, np.array([n * a] * 3 + [90.0] * 3)


# ----------------------------------------------------------------------
# Species and sites
# ----------------------------------------------------------------------
def test_hill_formula():
    assert hill_formula(["O", "H", "H"]) == "H2O"
    assert hill_formula(["C", "O", "O", "O", "C", "C", "H", "H", "H", "H"]) == "C3H4O3"
    assert hill_formula(["F", "B", "F", "F", "F"]) == "BF4"
    assert hill_formula(["Li"]) == "Li"


@pytest.fixture
def nacl_water():
    """Two waters, a sodium and a chloride."""
    symbols = ["O", "H", "H", "Na", "O", "H", "H", "Cl"]
    molecules = [[0, 1, 2], [3], [4, 5, 6], [7]]
    return Species(symbols, molecules)


def test_species_and_sites(nacl_water):
    s = nacl_water
    assert s.names == ["H2O", "Na", "Cl"]
    assert s.count("H2O") == 2
    assert s.site("O(H2O)").tolist() == [0, 4]
    assert s.site("O").tolist() == [0, 4]
    assert s.site("H2O").tolist() == [0, 1, 2, 4, 5, 6]
    assert s.site("Na").tolist() == [3]
    assert s.heavy_atoms("H2O").tolist() == [0, 4]
    assert s.element_pairs() == [
        ("Cl", "Cl"),
        ("Cl", "H"),
        ("Cl", "Na"),
        ("Cl", "O"),
        ("H", "H"),
        ("H", "Na"),
        ("H", "O"),
        ("Na", "Na"),
        ("Na", "O"),
        ("O", "O"),
    ]
    assert s.site_label(s.site("O(H2O)")) == "O(H2O)"
    with pytest.raises(AnalysisError, match="no species"):
        s.site("O(EC)")
    with pytest.raises(AnalysisError, match="No atoms"):
        s.site("Li")
    with pytest.raises(AnalysisError, match="selection needs"):
        s.site("select: name O")


def test_ions(nacl_water):
    assert nacl_water.ions() == (["Na"], ["Cl"])
    charges = [-0.8476, 0.4238, 0.4238, 1.0, -0.8476, 0.4238, 0.4238, -1.0]
    assert nacl_water.ions(charges) == (["Na"], ["Cl"])
    # A polyatomic anion is found from the charges
    s = Species(["Li", "B", "F", "F", "F", "F"], [[0], [1, 2, 3, 4, 5]])
    assert s.ions() == (["Li"], [])
    assert s.ions([1.0, 0.8, -0.45, -0.45, -0.45, -0.45]) == (["Li"], ["BF4"])


# ----------------------------------------------------------------------
# Radial distribution functions
# ----------------------------------------------------------------------
def test_lattice_coordination_numbers():
    """A simple cubic lattice: 6 neighbours at a, 12 at √2 a, 8 at √3 a."""
    xyz, box = cubic_lattice()
    n = len(xyz)
    rdf = RDF(np.arange(n), np.arange(n), np.arange(n), r_max=8.0, dr=0.05)
    rdf.add_frame(xyz, box=box, volume=box[0] ** 3)
    res = rdf.result(n_blocks=1)
    assert [coordination_number(res, c)[0] for c in (3.1, 4.3, 5.3)] == [
        6.0,
        18.0,
        26.0,
    ]
    assert res["g_err"] is None and res["n_blocks"] == 1


def test_ideal_gas():
    """Uncorrelated points: g(r) = 1, n(r) = ρ 4/3 π r^3, with honest errors."""
    rng = np.random.default_rng(1)
    L, n = 20.0, 1000
    box = np.array([L] * 3 + [90.0] * 3)
    rdf = RDF(np.arange(n), np.arange(n), np.arange(n), r_max=8.0, dr=0.1)
    for _ in range(20):
        rdf.add_frame(rng.uniform(0, L, (n, 3)), box=box, volume=L**3)
    res = rdf.result(n_blocks=5)
    assert res["g"][20:].mean() == pytest.approx(1.0, abs=0.01)
    cn, err = coordination_number(res, 5.0)
    expected = (n - 1) / L**3 * 4 / 3 * math.pi * 5.0**3
    assert abs(cn - expected) < 4 * err
    # The cross g(r) of two halves
    half = RDF(np.arange(500), np.arange(500, 1000), np.arange(n), r_max=8.0, dr=0.1)
    for _ in range(10):
        half.add_frame(rng.uniform(0, L, (n, 3)), box=box, volume=L**3)
    assert half.result()["g"][20:].mean() == pytest.approx(1.0, abs=0.02)


def test_pairs_within_a_molecule_are_left_out():
    """Two-atom molecules: the partner at 1 Å is not a neighbour."""
    xyz = np.array([[0, 0, 0], [1, 0, 0], [5, 0, 0], [6, 0, 0]], float)
    molecules = np.array([0, 0, 1, 1])
    rdf = RDF(np.arange(4), np.arange(4), molecules, r_max=6.5, dr=0.5)
    rdf.add_frame(xyz)
    res = rdf.result()
    assert res["g"] is None  # no cell, so only the counts
    assert coordination_number(res, 3.0)[0] == 0.0
    assert coordination_number(res, 6.5)[0] == 2.0  # 4, 5 (and 5, 6) Å away


def test_peak_and_minimum():
    r = np.arange(0.01, 10.0, 0.02)
    g = np.where(r < 2.5, 0.0, 1.0 + 2.0 * np.exp(-(((r - 3.0) / 0.3) ** 2)))
    g -= 0.6 * np.exp(-(((r - 4.0) / 0.3) ** 2)) * (r > 2.5)
    g += 0.4 * np.exp(-(((r - 5.2) / 0.4) ** 2)) * (r > 2.5)
    g -= 0.2 * np.exp(-(((r - 6.4) / 0.4) ** 2)) * (r > 2.5)
    pm = peak_and_minimum(r, g)
    assert pm["peak"] == pytest.approx(3.0, abs=0.05)
    assert pm["minimum"] == pytest.approx(4.0, abs=0.1)
    assert second_minimum(r, g, pm["minimum"]) == pytest.approx(6.4, abs=0.15)
    # A g(r) that never rises above 1 has no peak
    assert peak_and_minimum(r, np.full_like(r, 0.5))["peak"] is None


def test_minimum_of_a_noisy_sparse_pair():
    """A few ions give single-bin spikes above 1 everywhere; they must not be
    taken for the next peak (run 12197: 9 NaCl, the minimum fell at 2.77 Å)."""
    rng = np.random.default_rng(3)
    r = np.arange(0.01, 8.0, 0.02)
    g = np.where(r < 2.4, 0.0, 0.3)
    g += 1.5 * np.exp(-(((r - 2.7) / 0.12) ** 2))
    g += np.where(r > 4.0, 1.0 + 1.2 * np.exp(-(((r - 4.8) / 0.4) ** 2)), 0.0)
    noisy = g + np.where(
        (r > 2.9) & (r < 4.0), rng.choice([0.0, 1.5], r.size, p=[0.8, 0.2]), 0
    )
    pm = peak_and_minimum(r, noisy)
    assert pm["peak"] == pytest.approx(2.7, abs=0.05)
    assert 3.0 < pm["minimum"] < 4.0


# ----------------------------------------------------------------------
# Ion pairing
# ----------------------------------------------------------------------
@pytest.mark.parametrize("periodic", [True, False])
def test_ion_pair_classes(periodic):
    cations = np.array(
        [[0, 0, 0], [10, 0, 0], [20, 0, 0], [30, 0, 0], [30, 3, 0]], float
    )
    anions = np.array([[2.5, 0, 0], [40, 40, 40], [24.5, 0, 0], [31.5, 1.5, 0]], float)
    box = np.array([100.0] * 3 + [90.0] * 3) if periodic else None
    classes = classify_ion_pairs(cations, anions, 3.0, 5.0, box=box)
    assert classes.tolist() == ["CIP", "free", "SSIP", "AGG", "AGG"]


def test_two_anions_in_contact_is_an_aggregate():
    classes = classify_ion_pairs(
        np.array([[0, 0, 0.0]]), np.array([[2.5, 0, 0], [-2.5, 0, 0.0]]), 3.0, 5.0
    )
    assert classes.tolist() == ["AGG"]


def test_polyatomic_anion_counts_once():
    """Li touching two F of one BF4 is one contact: CIP, not AGG."""
    cation = np.array([[0, 0, 0.0]])
    bf4 = np.array([[3.0, 0, 0], [2.0, 0.6, 0], [2.0, -0.6, 0], [3.0, 1.2, 0]])
    classes = classify_ion_pairs(cation, bf4, 2.5, 5.0, anion_ids=[0, 0, 0, 0])
    assert classes.tolist() == ["CIP"]


def test_pairing_fractions():
    acc = PairingAccumulator()
    for _ in range(10):
        acc.add_frame(np.array(["CIP", "free", "free", "SSIP"], dtype=object))
    r = acc.result(n_blocks=5)
    assert r["fractions"] == {"free": 0.5, "SSIP": 0.25, "CIP": 0.25, "AGG": 0.0}
    assert r["errors"]["free"] == pytest.approx(0.0)


# ----------------------------------------------------------------------
# Solvation shells
# ----------------------------------------------------------------------
def test_shells_and_speciation():
    """A cation with two waters (one oxygen each) inside 3 Å and one outside,
    and an anion within 3 Å in half the frames."""
    water_o = np.array([[2.0, 0, 0], [-2.0, 0, 0], [0, 6.0, 0]])
    acc = ShellAccumulator(
        np.array([0]),
        np.array([0]),
        {
            "H2O": (np.array([1, 2, 3]), np.array([10, 11, 12]), 3.0),
            "Cl": (np.array([4]), np.array([20]), 3.5),
        },
    )
    for t in range(10):
        x = np.zeros((5, 3))
        x[1:4] = water_o
        x[4] = [0, 0, 3.0] if t % 2 == 0 else [0, 0, 8.0]
        acc.add_frame(x)
    r = acc.result(n_blocks=1, dt=0.1)
    assert r["coordination"] == {"H2O": 2.0, "Cl": 0.5}
    assert r["pairing"] == {"H2O": 1.0, "Cl": 0.5}
    shells = {tuple(sorted(c.items())): f for c, f in r["speciation"]}
    assert shells == {(("H2O", 2),): 0.5, (("Cl", 1), ("H2O", 2)): 0.5}
    # The waters never leave
    assert r["residence"]["H2O"] is None


def test_residence_time():
    # Each pair stays 10 frames, then is replaced
    members = [{(0, t // 10)} for t in range(200)]
    tau = residence_time(members, 0.1)
    assert 0.4 < tau < 0.8
    assert residence_time([set()] * 10, 0.1) is None
    assert residence_time([{(0, 1)}] * 3, 0.1) is None


# ----------------------------------------------------------------------
# Parsing the parameters
# ----------------------------------------------------------------------
def test_parse_pairs():
    assert parse_pairs("O(H2O)-O(H2O), Na-Cl , Li - select: name O") == [
        ("O(H2O)", "O(H2O)"),
        ("Na", "Cl"),
        ("Li", "select: name O"),
    ]
    assert parse_pairs("") == []
    with pytest.raises(AnalysisError, match="joined by"):
        split_pair("NaCl")


def test_parse_cutoffs():
    assert parse_cutoffs("Na-O(H2O): 3.2, Na-Cl: 3.5 5.6") == {
        ("Na", "O(H2O)"): [3.2],
        ("Na", "Cl"): [3.5, 5.6],
    }
    for bad in ("Na-Cl 3.5", "Na-Cl: x", "Na-Cl: 1 2 3", "Na-Cl: -1"):
        with pytest.raises(AnalysisError):
            parse_cutoffs(bad)
