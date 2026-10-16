# -*- coding: utf-8 -*-

"""The structural analysis of the Liquid Structure step, independent of SEAMM.

- Species and sites: each molecule is named by its formula, and a *site* is an
  element of a species ("O(H2O)"), an element in every species ("O"), or an
  MDAnalysis selection ("select: name O and resname 1").
- Radial distribution functions accumulated frame by frame, with errors from
  blocks of frames; the first peak and minimum; coordination numbers.
- The classification of cation-anion pairs as free, solvent-separated (SSIP),
  contact (CIP) or aggregated (AGG), from distances only, so it also works without
  a periodic cell.
"""

from collections import Counter
import math
import re

import numpy as np

#: The order of elements in a formula after C and H (Hill order)
_SITE = re.compile(r"^\s*([A-Z][a-z]?)\s*(?:\(\s*([^()]+?)\s*\))?\s*$")


class AnalysisError(ValueError):
    """A problem with what was asked for."""


# ----------------------------------------------------------------------
# Species and sites
# ----------------------------------------------------------------------
def hill_formula(symbols):
    """The formula of a molecule in Hill order: C, then H, then alphabetical; with
    no carbon, all alphabetical."""
    counts = Counter(symbols)
    if "C" in counts:
        order = ["C"] + (["H"] if "H" in counts else [])
        order += sorted(e for e in counts if e not in ("C", "H"))
    else:
        order = sorted(counts)
    return "".join(e + (str(counts[e]) if counts[e] > 1 else "") for e in order)


class Species:
    """The species of a system: the formula of each molecule and of each atom.

    Parameters
    ----------
    symbols : [str]
        The element of each atom.
    molecules : [[int]]
        The atoms (0-based) of each molecule.
    """

    def __init__(self, symbols, molecules):
        self.symbols = np.array(symbols)
        self.n_atoms = len(symbols)
        self.molecule_of_atom = np.full(self.n_atoms, -1, dtype=int)
        self.molecule_species = []
        for i, molecule in enumerate(molecules):
            self.molecule_of_atom[molecule] = i
            self.molecule_species.append(hill_formula(self.symbols[molecule].tolist()))
        # Atoms in no molecule are their own
        for i in np.where(self.molecule_of_atom < 0)[0]:
            self.molecule_of_atom[i] = len(self.molecule_species)
            self.molecule_species.append(str(self.symbols[i]))
        self.atom_species = np.array(
            [self.molecule_species[m] for m in self.molecule_of_atom]
        )
        self.molecules = [
            np.where(self.molecule_of_atom == m)[0]
            for m in range(len(self.molecule_species))
        ]

    @property
    def names(self):
        """The species, in order of appearance."""
        return list(dict.fromkeys(self.molecule_species))

    def count(self, species):
        return self.molecule_species.count(species)

    def site(self, text, universe=None):
        """The atoms (0-based indices) of a site.

        Parameters
        ----------
        text : str
            "O(H2O)", an element of a species; "O", the element in every species;
            "H2O", every atom of a species; or "select: <MDAnalysis selection>".
        universe : MDAnalysis.Universe
            Needed for a selection.

        Returns
        -------
        numpy.ndarray of int
        """
        text = text.strip()
        if text.startswith("select:"):
            if universe is None:
                raise AnalysisError(f"A selection needs the trajectory: {text}")
            return universe.select_atoms(text[len("select:") :].strip()).indices
        if text in self.names:
            return np.where(self.atom_species == text)[0]
        m = _SITE.match(text)
        if m is None:
            raise AnalysisError(
                f"'{text}' is not a site: give an element of a species, e.g. "
                "'O(H2O)', an element, a species, or 'select: <selection>'."
            )
        element, species = m.groups()
        mask = self.symbols == element
        if species is not None:
            if species not in self.names:
                raise AnalysisError(
                    f"There is no species {species}; the species are "
                    f"{', '.join(self.names)}."
                )
            mask &= self.atom_species == species
        atoms = np.where(mask)[0]
        if atoms.size == 0:
            raise AnalysisError(f"No atoms are '{text}'.")
        return atoms

    def heavy_atoms(self, species):
        """The atoms of a species other than hydrogens (all of them if it has only
        hydrogens): its coordinating atoms, for shells and ion pairing."""
        atoms = np.where(self.atom_species == species)[0]
        if atoms.size == 0:
            raise AnalysisError(f"There is no species {species}.")
        heavy = atoms[self.symbols[atoms] != "H"]
        return heavy if heavy.size > 0 else atoms

    def site_label(self, atoms):
        """A short label for a set of atoms: 'O(H2O)' if they are one element of
        one species, else the element or the species."""
        elements = set(self.symbols[atoms].tolist())
        species = set(self.atom_species[atoms].tolist())
        if len(elements) == 1 and len(species) == 1:
            (e,) = elements
            (s,) = species
            return e if e == s else f"{e}({s})"
        if len(elements) == 1:
            return next(iter(elements))
        if len(species) == 1:
            return next(iter(species))
        return "selection"

    def element_pairs(self):
        """Every pair of elements, each pair once."""
        elements = sorted(set(self.symbols.tolist()))
        return [(a, b) for i, a in enumerate(elements) for b in elements[i:]]

    def ions(self, charges=None):
        """The monatomic species, split into cations and anions by their charges
        (or, without charges, by element: the alkali and alkaline-earth metals are
        cations, the halogens anions); and the polyatomic ions (anions such as
        BF4) when the charges show them.

        Returns
        -------
        ([str], [str])
            The cation and anion species.
        """
        cations, anions = [], []
        for name in self.names:
            atoms = np.where(self.atom_species == name)[0]
            n = self.count(name)
            if charges is not None:
                q = float(np.sum(np.asarray(charges)[atoms])) / max(n, 1)
                if q > 0.5:
                    cations.append(name)
                elif q < -0.5:
                    anions.append(name)
            elif len(atoms) == n:  # monatomic
                if name in CATIONS:
                    cations.append(name)
                elif name in ANIONS:
                    anions.append(name)
        return cations, anions


#: Monatomic species taken as cations or anions when there are no charges
CATIONS = {"Li", "Na", "K", "Rb", "Cs", "Mg", "Ca", "Sr", "Ba", "Zn"}
ANIONS = {"F", "Cl", "Br", "I"}


# ----------------------------------------------------------------------
# Radial distribution functions
# ----------------------------------------------------------------------
class RDF:
    """A radial distribution function accumulated over frames.

    Pairs within the same molecule are left out.

    Parameters
    ----------
    a, b : numpy.ndarray of int
        The atoms of the two sites.
    molecule_of_atom : numpy.ndarray of int
        The molecule of each atom.
    r_max : float
        The range, Å.
    dr : float
        The width of the bins, Å.
    """

    def __init__(self, a, b, molecule_of_atom, r_max=12.0, dr=0.02):
        self.a = np.asarray(a)
        self.b = np.asarray(b)
        self.same = np.array_equal(np.sort(self.a), np.sort(self.b))
        self.mol_a = molecule_of_atom[self.a]
        self.mol_b = molecule_of_atom[self.b]
        self.r_max = float(r_max)
        self.dr = float(dr)
        self.n_bins = int(round(self.r_max / self.dr))
        self.edges = np.linspace(0.0, self.n_bins * self.dr, self.n_bins + 1)
        self.frames = []  # per frame: normalized g(r) and the counts per a atom

    def add_frame(self, positions, box=None, volume=None):
        """Add a frame.

        Parameters
        ----------
        positions : numpy.ndarray(n_atoms, 3)
            All the positions, Å.
        box : numpy.ndarray(6) or None
            The MDAnalysis box (a, b, c, alpha, beta, gamma), or None without a cell.
        volume : float or None
            The volume of the cell, Å^3; needed for g(r), not for the counts.
        """
        from MDAnalysis.lib.distances import capped_distance

        pa = positions[self.a]
        pb = positions[self.b]
        pairs, distances = capped_distance(
            pa, pb, self.r_max, box=box, return_distances=True
        )
        keep = self.mol_a[pairs[:, 0]] != self.mol_b[pairs[:, 1]]
        hist, _ = np.histogram(distances[keep], bins=self.edges)
        # The neighbours of each a atom in each shell. With the same atoms in both
        # sites each pair is found twice, as (i, j) and (j, i), which is right:
        # each is a neighbour of the other.
        counts = hist.astype(float) / len(self.a)
        g = None
        if volume is not None:
            n_b = len(self.b)
            shell = 4.0 / 3.0 * math.pi * (self.edges[1:] ** 3 - self.edges[:-1] ** 3)
            density = (n_b - 1 if self.same else n_b) / volume
            g = counts / (density * shell)
            self.density = density
        self.frames.append((g, counts))

    @property
    def r(self):
        """The centres of the bins, Å."""
        return 0.5 * (self.edges[1:] + self.edges[:-1])

    def result(self, n_blocks=5):
        """The g(r) and running coordination number, with errors from blocks.

        Returns
        -------
        dict
            ``r``; ``g`` and ``g_err`` (None without a cell); ``n`` and ``n_err``,
            the running coordination number (neighbours within r); ``n_frames``
            and ``n_blocks``.
        """
        if len(self.frames) == 0:
            raise AnalysisError("No frames were analysed.")
        counts = np.array([c for _, c in self.frames])
        n = np.cumsum(counts, axis=1)
        result = {"r": self.r, "r_upper": self.edges[1:], "n_frames": len(self.frames)}
        result["n"], result["n_err"], used = _mean_and_block_error(n, n_blocks)
        if self.frames[0][0] is not None:
            g = np.array([x for x, _ in self.frames])
            result["g"], result["g_err"], used = _mean_and_block_error(g, n_blocks)
        else:
            result["g"] = result["g_err"] = None
        result["n_blocks"] = used
        return result


def _mean_and_block_error(rows, n_blocks):
    """The mean over frames and its standard error from the means of blocks."""
    mean = rows.mean(axis=0)
    n_blocks = min(int(n_blocks), rows.shape[0])
    if n_blocks < 2:
        return mean, None, 1
    blocks = np.array([b.mean(axis=0) for b in np.array_split(rows, n_blocks)])
    err = blocks.std(axis=0, ddof=1) / math.sqrt(n_blocks)
    return mean, err, n_blocks


def _smooth(g, n):
    """A moving average over n bins, normalized at the ends rather than padded
    with zeros."""
    if n <= 1:
        return g
    kernel = np.ones(n)
    return np.convolve(g, kernel, mode="same") / np.convolve(
        np.ones_like(g), kernel, mode="same"
    )


def _first_run(mask, n):
    """The first index where ``mask`` is true for ``n`` bins in a row, or None."""
    if mask.size < n:
        return None
    runs = np.convolve(mask.astype(int), np.ones(n, dtype=int), mode="valid")
    hits = np.where(runs == n)[0]
    return int(hits[0]) if hits.size > 0 else None


def peak_and_minimum(r, g, width=0.2, sustain=0.2):
    """The first peak and first minimum of a g(r).

    The g(r) is smoothed over ``width`` Å, and a crossing of 1 counts only when it
    lasts ``sustain`` Å, so that the noise of a sparse pair (a few ions) is not
    taken for structure. The first peak is the highest point (of the unsmoothed
    g(r)) of the first stretch above 1; the first minimum the lowest point of the
    smoothed g(r) between it and the next stretch above 1.

    Returns
    -------
    dict
        ``peak`` and ``peak_height``, ``minimum`` and ``minimum_height``; None
        where they cannot be found.
    """
    r = np.asarray(r)
    g = np.asarray(g, dtype=float)
    result = {
        "peak": None,
        "peak_height": None,
        "minimum": None,
        "minimum_height": None,
    }
    if r.size < 3:
        return result
    dr = float(r[1] - r[0])
    smooth = max(1, int(round(width / dr)))
    n = max(1, int(round(sustain / dr)))
    gs = _smooth(g, smooth)
    start = _first_run(gs > 1.0, n)
    if start is None:
        return result
    below = _first_run(gs[start:] < 1.0, n)
    end = start + below if below is not None else len(gs)
    i_peak = start + int(np.argmax(g[start:end]))
    result["peak"] = float(r[i_peak])
    result["peak_height"] = float(g[i_peak])
    if below is None:
        return result
    rise = _first_run(gs[end:] > 1.0, n)
    stop = end + rise if rise is not None else len(gs)
    i_min = i_peak + int(np.argmin(gs[i_peak:stop]))
    if i_min >= len(gs) - max(1, smooth // 2):
        return result  # still falling at the end of the range
    result["minimum"] = float(r[i_min])
    result["minimum_height"] = float(gs[i_min])
    return result


def second_minimum(r, g, first, width=0.2, sustain=0.2):
    """The second minimum of a g(r): the lowest point between the peak after the
    first minimum and the next stretch above 1, or None."""
    r = np.asarray(r)
    if first is None:
        return None
    i0 = int(np.searchsorted(r, first))
    rest = peak_and_minimum(r[i0:], np.asarray(g)[i0:], width=width, sustain=sustain)
    return rest["minimum"]


def coordination_number(rdf_result, cutoff):
    """The coordination number within ``cutoff`` and its error, from the running
    coordination number: the neighbours in the bins wholly within the cutoff."""
    upper = rdf_result["r_upper"]
    i = int(np.searchsorted(upper, cutoff + 1.0e-9, side="right")) - 1
    if i < 0:
        return 0.0, 0.0
    i = min(i, len(upper) - 1)
    n = float(rdf_result["n"][i])
    err = None if rdf_result["n_err"] is None else float(rdf_result["n_err"][i])
    return n, err


# ----------------------------------------------------------------------
# Ion pairing
# ----------------------------------------------------------------------
PAIR_CLASSES = ("free", "SSIP", "CIP", "AGG")


def classify_ion_pairs(cations, anions, r1, r2, box=None, anion_ids=None):
    """Classify each cation of one frame.

    Parameters
    ----------
    cations : numpy.ndarray(n, 3)
        The positions of the cations.
    anions : numpy.ndarray(m, 3)
        The positions of the anions' atoms: one row per anion for monatomic anions,
        or every atom of polyatomic ones, with ``anion_ids`` giving the anion of
        each row. An anion is within a cutoff when any of its atoms is.
    r1 : float
        The contact cutoff, Å (the first minimum of their g(r)).
    r2 : float
        The solvent-separated cutoff, Å (the second minimum).
    box : numpy.ndarray(6) or None

    Returns
    -------
    numpy.ndarray of str
        For each cation: "free" (no anion within r2), "SSIP" (none within r1, one
        or more within r2), "CIP" (exactly one within r1, in contact with no other
        cation) or "AGG" (an anion within r1 that is also within r1 of another
        cation, or more than one within r1).
    """
    from MDAnalysis.lib.distances import capped_distance

    n_c = len(cations)
    classes = np.full(n_c, "free", dtype=object)
    if n_c == 0 or len(anions) == 0:
        return classes
    if anion_ids is None:
        anion_ids = np.arange(len(anions))
    anion_ids = np.asarray(anion_ids)
    n_a = int(anion_ids.max()) + 1
    pairs, d = capped_distance(cations, anions, max(r1, r2), box=box)
    # One entry per cation and anion (not per atom of the anion)
    contact = np.unique(
        np.column_stack([pairs[d < r1, 0], anion_ids[pairs[d < r1, 1]]]), axis=0
    ).reshape(-1, 2)
    near = np.unique(
        np.column_stack([pairs[d < r2, 0], anion_ids[pairs[d < r2, 1]]]), axis=0
    ).reshape(-1, 2)
    n_contact = np.bincount(contact[:, 0], minlength=n_c)
    n_near = np.bincount(near[:, 0], minlength=n_c)
    cations_per_anion = np.bincount(contact[:, 1], minlength=n_a)
    classes[n_near > 0] = "SSIP"
    shared = np.zeros(n_c, dtype=bool)
    for c, a in contact:
        if cations_per_anion[a] > 1:
            shared[c] = True
    classes[n_contact == 1] = "CIP"
    classes[(n_contact > 1) | shared] = "AGG"
    return classes


class PairingAccumulator:
    """The fractions of cations in each pairing class, over frames, with errors
    from blocks."""

    def __init__(self):
        self.frames = []

    def add_frame(self, classes):
        n = max(len(classes), 1)
        self.frames.append([np.count_nonzero(classes == c) / n for c in PAIR_CLASSES])

    def result(self, n_blocks=5):
        rows = np.array(self.frames)
        mean, err, used = _mean_and_block_error(rows, n_blocks)
        return {
            "fractions": dict(zip(PAIR_CLASSES, mean.tolist())),
            "errors": None if err is None else dict(zip(PAIR_CLASSES, err.tolist())),
            "n_frames": len(self.frames),
            "n_blocks": used,
        }


# ----------------------------------------------------------------------
# Solvation shells
# ----------------------------------------------------------------------
class ShellAccumulator:
    """The solvation shells of the molecules of one solute species.

    A solvent molecule is in a solute's shell when any of its atoms is within the
    cutoff of any of the solute's atoms.

    Parameters
    ----------
    solute_atoms : numpy.ndarray of int
        The atoms of the solute molecules.
    solute_ids : numpy.ndarray of int
        The solute molecule (0, 1, ...) of each of those atoms.
    solvents : {str: (numpy.ndarray of int, numpy.ndarray of int, float)}
        For each solvent species: its atoms, the molecule of each, and the cutoff
        (Å).
    """

    def __init__(self, solute_atoms, solute_ids, solvents):
        self.solute_atoms = np.asarray(solute_atoms)
        self.solute_ids = np.asarray(solute_ids)
        self.n_solutes = int(self.solute_ids.max()) + 1
        self.solvents = solvents
        self.names = list(solvents)
        self.counts = []  # per frame: (n_solutes, n_solvents)
        self.members = {name: [] for name in self.names}  # per frame: set of pairs

    def add_frame(self, positions, box=None):
        from MDAnalysis.lib.distances import capped_distance

        p_solute = positions[self.solute_atoms]
        counts = np.zeros((self.n_solutes, len(self.names)), dtype=int)
        for k, name in enumerate(self.names):
            atoms, molecules, cutoff = self.solvents[name]
            pairs = capped_distance(
                p_solute, positions[atoms], cutoff, box=box, return_distances=False
            )
            pairs = np.unique(
                np.column_stack([self.solute_ids[pairs[:, 0]], molecules[pairs[:, 1]]]),
                axis=0,
            ).reshape(-1, 2)
            counts[:, k] = np.bincount(pairs[:, 0], minlength=self.n_solutes)
            self.members[name].append({(int(a), int(b)) for a, b in pairs})
        self.counts.append(counts)

    def result(self, n_blocks=5, dt=None, top=10):
        """The coordination numbers, shell compositions, pairing and (with ``dt``,
        the interval between frames in ps) residence times.

        Returns
        -------
        dict
            ``coordination`` and ``coordination_err``: {solvent: number};
            ``pairing``: {solvent: fraction of solutes with at least one};
            ``speciation``: [(composition, fraction)], most common first, the
            composition as {solvent: count}; ``residence``: {solvent: ps} or
            None.
        """
        counts = np.array(self.counts)  # frames x solutes x solvents
        cn = counts.mean(axis=1)  # frames x solvents
        mean, err, used = _mean_and_block_error(cn, n_blocks)
        paired = (counts > 0).mean(axis=1)
        compositions = Counter()
        for frame in counts:
            for row in frame:
                compositions[tuple(int(x) for x in row)] += 1
        total = sum(compositions.values())
        speciation = [
            (
                {n: c for n, c in zip(self.names, comp) if c > 0},
                count / total,
            )
            for comp, count in compositions.most_common(top)
        ]
        result = {
            "coordination": dict(zip(self.names, mean.tolist())),
            "coordination_err": (
                None if err is None else dict(zip(self.names, err.tolist()))
            ),
            "pairing": dict(zip(self.names, paired.mean(axis=0).tolist())),
            "speciation": speciation,
            "n_frames": counts.shape[0],
            "n_blocks": used,
            "residence": None,
        }
        if dt is not None:
            result["residence"] = {
                name: residence_time(self.members[name], dt) for name in self.names
            }
        return result


def residence_time(members, dt):
    """The residence time of solvent molecules in shells, ps: the time at which the
    normalized autocorrelation of the occupancy, <h(0) h(t)> / <h>, first falls
    below 1/e (interpolated), or None if it does not within half the run.

    Parameters
    ----------
    members : [set of (solute, solvent molecule)]
        The pairs in contact in each frame.
    dt : float
        The interval between frames, ps.
    """
    n = len(members)
    if n < 4:
        return None
    pairs = sorted(set().union(*members))
    if len(pairs) == 0:
        return None
    index = {p: i for i, p in enumerate(pairs)}
    h = np.zeros((len(pairs), n), dtype=float)
    for t, frame in enumerate(members):
        for p in frame:
            h[index[p], t] = 1.0
    max_lag = n // 2
    occupancy = h.sum()
    if occupancy == 0:
        return None
    c = np.array(
        [np.sum(h[:, : n - lag] * h[:, lag:]) / (n - lag) for lag in range(max_lag)]
    )
    c /= c[0]
    target = 1.0 / math.e
    below = np.where(c < target)[0]
    if below.size == 0:
        return None
    i = below[0]
    if i == 0:
        return 0.0
    # Linear interpolation between lags i-1 and i
    frac = (c[i - 1] - target) / (c[i - 1] - c[i])
    return float((i - 1 + frac) * dt)


# ----------------------------------------------------------------------
# Parsing the parameters
# ----------------------------------------------------------------------
def split_pair(text):
    """Two sites joined by ' - ' or, when neither contains a hyphen, '-'."""
    text = text.strip()
    if " - " in text:
        a, b = text.split(" - ", 1)
    else:
        depth = 0
        cut = None
        for i, c in enumerate(text):
            if c == "(":
                depth += 1
            elif c == ")":
                depth -= 1
            elif c == "-" and depth == 0:
                cut = i
                break
        if cut is None:
            raise AnalysisError(f"'{text}' is not a pair of sites joined by '-'.")
        a, b = text[:cut], text[cut + 1 :]
    a, b = a.strip(), b.strip()
    if a == "" or b == "":
        raise AnalysisError(f"'{text}' is not a pair of sites joined by '-'.")
    return a, b


def parse_pairs(text):
    """The pairs in a comma-separated list: [(site, site)]."""
    return [split_pair(item) for item in (text or "").split(",") if item.strip()]


def parse_cutoffs(text):
    """Fixed cutoffs: 'Na-O(H2O): 3.2, Na-Cl: 3.5 5.6' -> {(a, b): [3.2]...}.

    The key is the pair of sites in the order given; one value is a cutoff, two a
    contact and solvent-separated cutoff for ion pairing.
    """
    result = {}
    for item in (text or "").split(","):
        if item.strip() == "":
            continue
        if ":" not in item:
            raise AnalysisError(f"'{item.strip()}' should be 'pair: cutoff'.")
        pair, values = item.rsplit(":", 1)
        try:
            numbers = [float(x) for x in values.split()]
        except ValueError:
            raise AnalysisError(f"'{values.strip()}' are not cutoffs in Å.")
        if len(numbers) not in (1, 2) or any(x <= 0 for x in numbers):
            raise AnalysisError(f"'{item.strip()}' needs one or two positive cutoffs.")
        result[split_pair(pair)] = numbers
    return result


def find_cutoff(cutoffs, a, b):
    """The fixed cutoffs for a pair, in either order, or None."""
    return cutoffs.get((a, b), cutoffs.get((b, a)))
