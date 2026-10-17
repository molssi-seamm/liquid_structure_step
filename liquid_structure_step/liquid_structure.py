# -*- coding: utf-8 -*-

"""Non-graphical part of the Liquid Structure step in a SEAMM flowchart

The structure of liquids and solutions from the positions sampled in MD: partial
radial distribution functions and coordination numbers, the classification of ion
pairs (free, solvent-separated, contact, aggregated), and the solvation shells of
solutes (coordination numbers per solvent, the shells' composition, pairing and
residence times). The analysis itself is in :mod:`liquid_structure_step.analysis`.

The MD is run before this step -- in a Loop for independent runs -- or read from
files; the step analyzes it. The errors are from the independent runs, or from
blocks of a single run.
"""

import importlib.resources
import logging
from pathlib import Path
import warnings

import numpy as np
from seamm_md_analysis import Estimate, open_md
from seamm_md_analysis.requirements import Need, Requirements
from seamm_md_analysis.stats import t_factor
from seamm_md_analysis.step_support import (
    check_md_flowchart,
    md_inputs,
    print_details,
    read_files,
    report_state,
    size_warnings,
    state_data,
)
from tabulate import tabulate

import liquid_structure_step
from liquid_structure_step.analysis import (
    AnalysisError,
    PAIR_CLASSES,
    PairingAccumulator,
    RDF,
    ShellAccumulator,
    Species,
    classify_ion_pairs,
    coordination_by_block,
    coordination_number,
    find_cutoff,
    parse_cutoffs,
    parse_pairs,
    peak_and_minimum,
    second_minimum,
    split_blocks,
)
import molsystem
import seamm
import seamm_util
from seamm_util import Q_
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

logger = logging.getLogger(__name__)
# MDAnalysis warns about optional packages (netCDF4) in the step's output
logging.getLogger("MDAnalysis").setLevel(logging.ERROR)
printer = printing.getPrinter("Liquid Structure")

# Add this module's properties to the standard properties
path = importlib.resources.files("liquid_structure_step") / "data"
csv_file = path / "properties.csv"
if path.exists():
    molsystem.add_properties_from_file(csv_file)

#: The interval between frames that the residence times need
RESIDENCE_INTERVAL = Q_(0.1, "ps")


#: The heading of the columns of errors in the tables
PM = "± (95%)"


def cells(estimate):
    """A value and its 95% interval as text for a table; empty without a value.

    A value without an error estimate is given to four significant figures; one
    whose error is negligible (the same in every sample) as Estimate gives it.
    """
    if estimate is None or estimate.value is None or not np.isfinite(estimate.value):
        return "", ""
    if estimate.ci95 is None:
        return f"{estimate.value:.4g}", ""
    return estimate.formatted()


def yes(value):
    """Whether a boolean parameter is set, whether given as a bool or as text: the
    values for the description are 'yes' or 'no'."""
    if isinstance(value, str):
        return value.strip().lower() in ("yes", "true", "1")
    return bool(value)


def stderr_of(samples):
    """The standard error of the finite values among samples, and how many there
    are; the error is None with fewer than two."""
    values = np.array([np.nan if x is None else x for x in samples], dtype=float)
    e = Estimate.from_samples(values)
    return e.stderr, e.n


class LiquidStructure(seamm.Node):
    """
    The non-graphical part of a Liquid Structure step in a flowchart.

    Attributes
    ----------
    parameters : LiquidStructureParameters
        The control parameters for Liquid Structure.

    See Also
    --------
    TkLiquidStructure,
    LiquidStructure, LiquidStructureParameters
    """

    def __init__(
        self,
        flowchart=None,
        title="Liquid Structure",
        extension=None,
        logger=logger,
    ):
        """A step for Liquid Structure in a SEAMM flowchart.

        Parameters
        ----------
        flowchart: seamm.Flowchart
            The non-graphical flowchart that contains this step.
        title: str
            The name displayed in the flowchart.
        extension: None
            Not yet implemented
        logger : Logger = logger
            The logger to use and pass to parent classes
        """
        logger.debug(f"Creating Liquid Structure {self}")
        super().__init__(
            flowchart=flowchart,
            title="Liquid Structure",
            extension=extension,
            module=__name__,
            logger=logger,
        )

        self._metadata = liquid_structure_step.metadata
        self.parameters = liquid_structure_step.LiquidStructureParameters()
        self._reset()

    @property
    def version(self):
        """The semantic version of this module."""
        return liquid_structure_step.__version__

    @property
    def git_revision(self):
        """The git version of this module."""
        return liquid_structure_step.__git_revision__

    def _reset(self):
        """Clear the data of earlier runs, e.g. from an earlier loop iteration."""
        self.paths = []  # The md.h5md file of each run
        self._source = "blocks"  # What the errors are from: "runs" or "blocks"
        self._n_samples = 1  # The number of runs or blocks

    def create_parser(self):
        """Setup the command-line / config file parser"""
        parser_name = "liquid-structure-step"
        parser = seamm_util.getParser()

        parser_exists = parser.exists(parser_name)
        super().create_parser(name=parser_name)
        if not parser_exists:
            parser.add_argument(
                parser_name,
                "--html",
                action="store_true",
                help="whether to write out html files for graphs, etc.",
            )
        return self.next()

    def description_text(self, P=None, short=False):
        """Create the text description of what this step will do.

        Parameters
        ----------
        P: dict
            An optional dictionary of the current values of the control
            parameters.
        short : bool = False
            Not used: the description is always short.

        Returns
        -------
        str
            A description of the current step.
        """
        if P is None:
            P = self.parameters.values_to_dict()

        what = "radial distribution functions and coordination numbers"
        if P["pairs"] == "the pairs given":
            what += f" of {P['pair list']}"
        else:
            what += " of every pair of elements"
        extras = []
        if yes(P["ion pairing"]):
            extras.append("the pairing of the ions")
        if yes(P["solvation shells"]):
            extras.append(
                "the solvation shells of "
                + (
                    P["solute list"]
                    if P["solutes"] == "the solutes given"
                    else "the monatomic ions"
                )
                + (", with residence times" if yes(P["residence times"]) else "")
            )
        if extras:
            what += ", and " + " and ".join(extras)

        source = P["md source"]
        text = f"Analyze the structure: {what}, from "
        if source == read_files:
            text += f"the MD in {P['md files']}."
        else:
            text += f"{source} earlier in the flowchart."
        return self.header + "\n" + __(text, indent=4 * " ").__str__()

    def run(self):
        """Run a Liquid Structure step.

        Returns
        -------
        seamm.Node
            The next node object in the flowchart.
        """
        next_node = super().run(printer)
        self._reset()

        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )

        printer.important(__(self.description_text(P), indent=self.indent))

        # Check what is asked for before reading the MD
        try:
            parse_cutoffs(P["cutoffs"])
            if P["pairs"] == "the pairs given":
                if len(parse_pairs(P["pair list"])) == 0:
                    raise AnalysisError("No pairs were given.")
        except AnalysisError as e:
            raise RuntimeError(f"'{self.title}': {e}") from None

        self.references.cite(
            raw=self._bibliography["MDAnalysis2011"],
            alias="MDAnalysis2011",
            module="liquid_structure_step",
            level=1,
            note="MDAnalysis, which reads the trajectory and finds the distances.",
        )
        self.references.cite(
            raw=self._bibliography["MDAnalysis2016"],
            alias="MDAnalysis2016",
            module="liquid_structure_step",
            level=1,
            note="MDAnalysis, which reads the trajectory and finds the distances.",
        )

        self.paths = md_inputs(self, P, name="all/position")
        printer.job(
            __(
                f"Analyzing {len(self.paths)} MD run"
                f"{'s' if len(self.paths) > 1 else ''}: "
                + ", ".join(str(p) for p in self.paths),
                indent=self.indent + 4 * " ",
            )
        )
        if len(self.paths) == 0:
            raise RuntimeError("No MD run finished, so there is nothing to analyze.")

        _, configuration = self.get_system_configuration()
        self.analyze(P=P, configuration=configuration)

        return next_node

    # ------------------------------------------------------------------
    # The analysis
    # ------------------------------------------------------------------
    def analyze(self, indent="", P=None, configuration=None, **kwargs):
        """Analyze, print, plot and store the structure.

        Parameters
        ----------
        P : dict
            The values of the parameters.
        configuration : molsystem._Configuration
            The configuration that was simulated: its atoms are in the order of the
            atoms in the MD files.
        """
        if configuration is None:
            raise RuntimeError("There is no structure to analyze.")
        species = Species(
            configuration.atoms.symbols, configuration.find_molecules(as_indices=True)
        )
        cutoffs = parse_cutoffs(P["cutoffs"])
        nruns = len(self.paths)

        # Look at the first file for the cell, the charges and the frames
        with open_md(self.paths[0]) as md:
            periodic = md.periodic
            charges = md.charges("all")
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                universe = md.universe(configuration)
            n_file = universe.atoms.n_atoms
        if n_file != species.n_atoms:
            raise RuntimeError(
                f"The MD has {n_file} atoms but the structure {species.n_atoms}."
            )

        # The pairs of sites
        try:
            if P["pairs"] == "the pairs given":
                site_pairs = parse_pairs(P["pair list"])
            else:
                site_pairs = species.element_pairs()
            pairs = [
                (a, b, species.site(a, universe), species.site(b, universe))
                for a, b in site_pairs
            ]
        except AnalysisError as e:
            raise RuntimeError(f"'{self.title}': {e}") from None

        # The ions, for the pairing and as default solutes
        cations, anions = species.ions(charges)
        solutes = []
        if P["solvation shells"]:
            if P["solutes"] == "the solutes given":
                solutes = [x.strip() for x in P["solute list"].split(",") if x.strip()]
                for name in solutes:
                    if name not in species.names:
                        raise RuntimeError(
                            f"'{self.title}': there is no solute {name}; the species "
                            f"are {', '.join(species.names)}."
                        )
            else:
                solutes = [
                    name
                    for name in cations + anions
                    if np.count_nonzero(species.atom_species == name)
                    == species.count(name)
                ]

        # The radial distribution functions needed: those asked for, and for the
        # cutoffs of the ion pairing and shells, cation-anion and solute-solvent
        r_max = P["r max"].m_as("Å")
        dr = P["bin width"].m_as("Å")
        rdfs = {}
        for a, b, atoms_a, atoms_b in pairs:
            rdfs[(a, b)] = RDF(atoms_a, atoms_b, species.molecule_of_atom, r_max, dr)
        ion_pairs = []
        if P["ion pairing"]:
            ion_pairs = [(c, a) for c in cations for a in anions]
        shell_pairs = [(s, x) for s in solutes for x in species.names if x != s]
        # These use the heavy atoms of each species, so that a molecule is in a
        # shell when its coordinating atoms are, not merely its hydrogens
        for a, b in ion_pairs + shell_pairs:
            if (a, b) not in rdfs:
                rdfs[(a, b)] = RDF(
                    species.heavy_atoms(a),
                    species.heavy_atoms(b),
                    species.molecule_of_atom,
                    r_max,
                    dr,
                )

        # The range is at most half the shortest width of the cell
        if periodic:
            for x, box, volume in self._frames(configuration, periodic, first=True):
                r_max = min(r_max, 0.5 * _shortest_width(box))
            rdfs = {
                key: RDF(rdf.a, rdf.b, species.molecule_of_atom, r_max, dr)
                for key, rdf in rdfs.items()
            }

        # The first pass over the frames: the radial distribution functions
        n_frames = 0
        for x, box, volume in self._frames(configuration, periodic):
            for rdf in rdfs.values():
                rdf.add_frame(x, box=box, volume=volume)
            n_frames += 1
        dt = self._dt
        run_frames = self._run_frames

        # The independent samples for the errors: the runs, or blocks of one run
        if nruns > 1:
            self._source = "runs"
            blocks = run_frames
        else:
            self._source = "blocks"
            blocks = max(1, int(P["number of blocks"]))
        self._n_samples = len(split_blocks(range(n_frames), blocks))

        text = f"Analyzed {n_frames} frames from {nruns} run{'s' if nruns > 1 else ''}"
        if dt is not None:
            text += f", every {dt:.3g} ps"
        text += ". " + self._error_note()
        printer.normal(__(text, indent=self.indent + 4 * " "))
        if not periodic:
            printer.normal(
                __(
                    "The system is not periodic, so there are no radial "
                    "distribution functions: the coordination numbers and ion "
                    "pairing use the cutoffs given.",
                    indent=self.indent + 4 * " ",
                )
            )
        printer.normal("")

        results = {}
        rdf_results = {key: rdf.result(blocks) for key, rdf in rdfs.items()}
        # The ion pairing and shells read the frames again, with the cutoffs
        frames = lambda: self._frames(configuration, periodic)  # noqa: E731
        results["coordination numbers"] = self._coordination(
            pairs, rdf_results, cutoffs, periodic
        )
        if ion_pairs:
            results["ion pairing"] = self._ion_pairing(
                ion_pairs, rdf_results, cutoffs, species, frames, blocks
            )
        if solutes:
            results["solvation shells"] = self._shells(
                P,
                solutes,
                rdf_results,
                cutoffs,
                species,
                frames,
                blocks,
                dt,
                run_frames,
            )

        state = report_state(self, printer, self.paths)

        self._plot(pairs, rdf_results)

        # The results, for variables, tables and properties
        self.model = self._find_model()
        data = {
            "n_frames": n_frames,
            "nruns": nruns,
            "error source": self._source,
            "n_samples": self._n_samples,
        }
        data.update(results)
        if state is not None:
            data.update(state_data(state))
        if P["store"] and configuration is not None:
            self._store(configuration, results)
        self.store_results(configuration=configuration, data=data)

    def _estimate(self, value, stderr=None):
        """An Estimate from a value and its standard error over the samples, or
        None without a value."""
        if value is None:
            return None
        n = self._n_samples if stderr is not None else 1
        return Estimate(float(value), stderr, n, source=self._source)

    def _error_note(self):
        """What the ± in the tables is."""
        return Estimate(0.0, 1.0, self._n_samples, source=self._source).describe_error()

    def _coordination(self, pairs, rdf_results, cutoffs, periodic):
        """The peaks, minima and coordination numbers of the pairs asked for."""
        peaks = []
        coordination = []
        by_sample = []
        results = {}
        notes = []
        for a, b, _, _ in pairs:
            label = f"{a}-{b}"
            res = rdf_results[(a, b)]
            fixed = find_cutoff(cutoffs, a, b)
            entry = {}
            pm = {"peak": None, "peak_height": None, "minimum": None}
            if res["g"] is not None:
                pm = peak_and_minimum(res["r"], res["g"])
                # The errors from the peak and minimum of each run or block
                by_block = [peak_and_minimum(res["r"], g) for g in res["g_blocks"]]
                row = [label]
                for key, name in (
                    ("peak", "peak"),
                    ("peak_height", "peak height"),
                    ("minimum", "minimum"),
                ):
                    err, _ = stderr_of([x[key] for x in by_block])
                    if pm[key] is None:
                        err = None
                    entry[name] = pm[key]
                    if err is not None:
                        entry[f"{name},stderr"] = err
                    row.extend(cells(self._estimate(pm[key], err)))
                peaks.append(row)
            cutoff = fixed[0] if fixed else pm["minimum"]
            if cutoff is not None:
                n, err = coordination_number(res, cutoff)
                by_sample.append([label, *coordination_by_block(res, cutoff)])
                entry.update({"cutoff": cutoff, "cutoff given": bool(fixed), "CN": n})
                if err is not None:
                    entry["CN,stderr"] = err
                if (
                    not fixed
                    and pm["peak"] is not None
                    and cutoff > 1.6 * pm["peak"] + 0.5
                ):
                    notes.append(
                        f"The first minimum of {label}, {cutoff:.2f} Å, is far beyond "
                        f"its peak at {pm['peak']:.2f} Å, so it may be a later "
                        "minimum (a solvent-shared peak, for example). Consider a "
                        "fixed cutoff."
                    )
                coordination.append(
                    [
                        label,
                        f"{cutoff:.2f}" + (" (given)" if fixed else ""),
                        *cells(self._estimate(n, err)),
                    ]
                )
            elif not periodic:
                notes.append(
                    f"{label}: no cutoff was given, so no coordination number."
                )
            else:
                notes.append(
                    f"{label}: the g(r) has no clear first minimum, so there is no "
                    "coordination number. Give a fixed cutoff for one."
                )
            results[label] = entry

        if peaks:
            self._print_table(
                ["Pair", "Peak (Å)", PM, "g(peak)", PM, "Minimum (Å)", PM],
                peaks,
                "The first peak and minimum of g(r)",
            )
        if coordination:
            self._print_table(
                ["Pair", "Cutoff (Å)", "CN", PM],
                coordination,
                "Coordination numbers",
                align=["left", "left", "decimal", "decimal"],
                note=(
                    "CN is the number of the second site within the cutoff of the "
                    "first, the first minimum of g(r) unless given."
                ),
            )
        for note in notes:
            printer.normal(__(note, indent=8 * " "))
        if notes:
            printer.normal("")
        self._print_samples(by_sample)
        return results

    def _print_samples(self, by_sample):
        """The coordination numbers in each run or block, in the step's own output
        (step.out)."""
        if len(by_sample) == 0 or self._n_samples < 2:
            return
        what = "Run" if self._source == "runs" else "Block"
        headers = ["Pair"] + [f"{what} {i}" for i in range(1, self._n_samples + 1)]
        rows = [[row[0], *(f"{x:.3f}" for x in row[1:])] for row in by_sample]
        text = tabulate(
            rows,
            headers=headers,
            tablefmt="rounded_outline",
            disable_numparse=True,
            colalign=["left"] + ["decimal"] * self._n_samples,
        )
        length = len(text.splitlines()[0])
        title = f"Coordination numbers by {what.lower()}".center(length)
        print_details(
            printer,
            str(__("\n" + title + "\n" + text + "\n", indent=8 * " ", wrap=False)),
        )

    def _frames(self, configuration, periodic, first=False):
        """The frames of every run: positions (Å), box and volume (Å^3), the box
        and volume None without a cell. Sets ``self._dt``, the interval between
        frames in ps, and ``self._run_frames``, the number of frames of each run."""
        dts = []
        run_frames = []
        for path in self.paths:
            with open_md(path) as md:
                # The interval between frames from the file's own time axis:
                # MDAnalysis's H5MD reader does not find it
                times = md.first["particles/all/position/time"]
                if times.shape[0] > 1:
                    unit = times.attrs.get("unit", "ps")
                    if isinstance(unit, bytes):
                        unit = unit.decode()
                    dts.append(Q_(float(times[1] - times[0]), unit).m_as("ps"))
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    u = md.universe(configuration)
                run_frames.append(0)
                for ts in u.trajectory:
                    box = ts.dimensions if periodic else None
                    volume = ts.volume if periodic else None
                    run_frames[-1] += 1
                    yield ts.positions, box, volume
                    if first:
                        break
            if first:
                break
        self._dt = float(np.median(dts)) if dts else None
        self._run_frames = run_frames

    def _ion_pairing(self, ion_pairs, rdf_results, cutoffs, species, frames, blocks):
        """Classify the cations of each cation-anion pair in every frame."""
        rows = []
        results = {}
        notes = []
        for cation, anion in ion_pairs:
            label = f"{cation}-{anion}"
            res = rdf_results[(cation, anion)]
            fixed = find_cutoff(cutoffs, cation, anion)
            r1 = r2 = None
            if fixed:
                r1 = fixed[0]
                r2 = fixed[1] if len(fixed) > 1 else None
            if res["g"] is not None:
                pm = peak_and_minimum(res["r"], res["g"])
                if r1 is None:
                    r1 = pm["minimum"]
                if r2 is None:
                    r2 = second_minimum(res["r"], res["g"], r1)
            if r1 is None:
                results[label] = {"error": "no contact cutoff"}
                notes.append(
                    f"{label}: the g(r) has no clear first minimum, so the pairs "
                    f"cannot be classified. Give the cutoffs, e.g. '{label}: 3.5 5.6'."
                )
                continue
            if r2 is None:
                # No solvent-separated class without a second cutoff
                r2 = r1
                notes.append(
                    f"{label}: the g(r) has no clear second minimum, so solvent-"
                    "separated pairs (SSIP) are not told apart from free cations. "
                    f"Give two cutoffs to separate them, e.g. '{label}: {r1:.1f} 5.6'."
                )
            cat_atoms = species.heavy_atoms(cation)
            an_atoms = species.heavy_atoms(anion)
            an_ids = np.unique(species.molecule_of_atom[an_atoms], return_inverse=True)[
                1
            ]
            acc = PairingAccumulator()
            for x, box, _ in frames():
                acc.add_frame(
                    classify_ion_pairs(
                        x[cat_atoms], x[an_atoms], r1, r2, box=box, anion_ids=an_ids
                    )
                )
            r = acc.result(blocks)
            entry = {"r1": r1, "r2": r2, "fractions": r["fractions"]}
            if r["errors"] is not None:
                entry["fractions,stderr"] = r["errors"]
            results[label] = entry
            for i, c in enumerate(PAIR_CLASSES):
                err = None if r["errors"] is None else r["errors"][c]
                rows.append(
                    [
                        label if i == 0 else "",
                        f"{r1:.2f}" if i == 0 else "",
                        f"{r2:.2f}" if i == 0 and r2 != r1 else "",
                        c,
                        *cells(self._estimate(r["fractions"][c], err)),
                    ]
                )
        if rows:
            self._print_table(
                ["Pair", "r1 (Å)", "r2 (Å)", "Class", "Fraction", PM],
                rows,
                "Ion pairing (fraction of the cations)",
                align=["left", "decimal", "decimal", "left", "decimal", "decimal"],
                note=(
                    "free: no anion within r2; SSIP: none within r1, one or more "
                    "within r2; CIP: one anion within r1, touching no other cation; "
                    "AGG: an anion within r1 shared with another cation, or more "
                    "than one."
                ),
            )
        for note in notes:
            printer.normal(__(note, indent=8 * " "))
        if notes:
            printer.normal("")
        return results

    def _shells(
        self,
        P,
        solutes,
        rdf_results,
        cutoffs,
        species,
        frames,
        blocks,
        dt,
        run_frames,
    ):
        """The solvation shells of each solute."""
        results = {}
        for solute in solutes:
            solute_atoms = species.heavy_atoms(solute)
            solute_ids = np.unique(
                species.molecule_of_atom[solute_atoms], return_inverse=True
            )[1]
            solvents = {}
            skipped = []
            for name in species.names:
                if name == solute:
                    continue
                fixed = find_cutoff(cutoffs, solute, name)
                cutoff = fixed[0] if fixed else None
                res = rdf_results.get((solute, name))
                if cutoff is None and res is not None and res["g"] is not None:
                    cutoff = peak_and_minimum(res["r"], res["g"])["minimum"]
                if cutoff is None:
                    skipped.append(name)
                    continue
                atoms = species.heavy_atoms(name)
                solvents[name] = (atoms, species.molecule_of_atom[atoms], cutoff)
            if len(solvents) == 0:
                results[solute] = {"error": "no cutoffs for any solvent"}
                printer.normal(
                    __(
                        f"There is no cutoff for any species around {solute}, so it "
                        "has no solvation shell. Give fixed cutoffs for it.\n",
                        indent=8 * " ",
                    )
                )
                continue
            acc = ShellAccumulator(solute_atoms, solute_ids, solvents)
            for x, box, _ in frames():
                acc.add_frame(x, box=box)
            residence_dt = dt if P["residence times"] else None
            r = acc.result(blocks, dt=residence_dt, runs=run_frames)

            entry = {
                "cutoffs": {name: v[2] for name, v in solvents.items()},
                "coordination": r["coordination"],
                "pairing": r["pairing"],
            }
            if r["coordination_err"] is not None:
                entry["coordination,stderr"] = r["coordination_err"]
                entry["pairing,stderr"] = r["pairing_err"]

            headers = ["Solvent", "Cutoff (Å)", "CN", PM, "Paired", PM]
            residence = r["residence"]
            if residence is not None:
                headers.extend(["Residence (ps)", PM])
                entry["residence"] = residence
                entry["residence,stderr"] = {}
            rows = []
            for name in solvents:
                err = r["coordination_err"]
                p_err = r["pairing_err"]
                row = [
                    name,
                    f"{entry['cutoffs'][name]:.2f}",
                    *cells(
                        self._estimate(
                            r["coordination"][name], None if err is None else err[name]
                        )
                    ),
                    *cells(
                        self._estimate(
                            r["pairing"][name], None if p_err is None else p_err[name]
                        )
                    ),
                ]
                if residence is not None:
                    t_err, n = stderr_of(r["residence_blocks"][name])
                    entry["residence,stderr"][name] = t_err
                    t = residence[name]
                    estimate = None
                    if t is not None:
                        estimate = Estimate(
                            t, t_err, n if t_err is not None else 1, source=self._source
                        )
                    row.extend(cells(estimate))
                rows.append(row)
            note = (
                "CN is the number of molecules of each species in the shell; Paired, "
                "the fraction of the solutes with at least one."
            )
            if residence is not None:
                note += (
                    " The residence time is when the autocorrelation of the "
                    "occupancy of the shell falls to 1/e; blank if it does not "
                    "within half a run or block."
                )
            self._print_table(
                headers,
                rows,
                f"The solvation shell of {solute}",
                note=note,
            )

            # The most common shells
            rows = []
            entry["speciation"] = []
            for composition, fraction, err in r["speciation"][:5]:
                shell = ", ".join(
                    f"{count} {name}" if count > 1 else name
                    for name, count in composition.items()
                )
                rows.append(
                    [
                        f"{solute}({shell or 'none'})",
                        *cells(self._estimate(fraction, err)),
                    ]
                )
            for composition, fraction, err in r["speciation"]:
                item = {"composition": composition, "fraction": fraction}
                if err is not None:
                    item["fraction,stderr"] = err
                entry["speciation"].append(item)
            self._print_table(
                ["Shell", "Fraction", PM],
                rows,
                f"The most common shells of {solute}",
            )
            results[solute] = entry

            notes = []
            if skipped:
                notes.append(
                    f"No cutoff for {', '.join(skipped)} around {solute}, so they are "
                    "left out."
                )
            if residence is not None and any(
                residence[name] is not None
                and entry["residence,stderr"][name] is None
                and self._n_samples > 1
                for name in solvents
            ):
                notes.append(
                    f"Some residence times have no error: they are too long to be "
                    f"found in each of the {self._source}. Use more or longer runs."
                )
            if P["residence times"] and dt is not None and dt > 0.15:
                notes.append(
                    f"The frames are {dt:.2g} ps apart; the residence times need "
                    f"frames every {RESIDENCE_INTERVAL:~P} or so."
                )
            for note in notes:
                printer.normal(__(note, indent=8 * " "))
            if notes:
                printer.normal("")
        return results

    def _print_table(self, headers, rows, title, note=None, align=None):
        """Print a table, its title above and an optional note below.

        Parameters
        ----------
        headers : [str]
            The headings of the columns.
        rows : [[str]]
            The rows of the table, as text.
        title : str
            The title, centred above the table.
        note : str
            A note for under the table.
        align : [str]
            The alignment of each column; by default the first left and the rest
            on the decimal point.
        """
        if align is None:
            align = ["left"] + ["decimal"] * (len(headers) - 1)
        text = tabulate(
            rows,
            headers=headers,
            tablefmt="rounded_outline",
            disable_numparse=True,
            colalign=align,
        )
        length = len(text.splitlines()[0])
        printer.normal(
            __(
                "\n" + title.center(length) + "\n" + text,
                indent=8 * " ",
                wrap=False,
                dedent=False,
            )
        )
        if note is not None:
            printer.normal(__(note, indent=8 * " "))
        printer.normal("")

    def _plot(self, pairs, rdf_results):
        """A graph of g(r) and the running coordination number for each pair."""
        for a, b, _, _ in pairs:
            res = rdf_results[(a, b)]
            r = res["r"].tolist()
            figure = self.create_figure(
                module_path=("seamm",),
                template="line.graph_template",
                title=f"{a}-{b}",
            )
            plot = figure.add_plot("rdf")
            x_axis = plot.add_axis("x", label="r (Å)")
            if res["g"] is not None:
                y_axis = plot.add_axis("y", label="g(r)", anchor=x_axis)
                x_axis.anchor = y_axis
                if res["g_err"] is not None:
                    # The 95% confidence interval, Student's t times the stderr
                    half = t_factor(res["n_blocks"]) * res["g_err"]
                    upper = (res["g"] + half).tolist()
                    lower = (res["g"] - half).tolist()
                    plot.add_trace(
                        x_axis=x_axis,
                        y_axis=y_axis,
                        name="95% interval",
                        x=r + r[::-1],
                        xlabel="r",
                        xunits="Å",
                        y=upper + lower[::-1],
                        ylabel="95% interval",
                        yunits="",
                        color="rgba(211,211,211,0.5)",
                        fill="toself",
                    )
                plot.add_trace(
                    x_axis=x_axis,
                    y_axis=y_axis,
                    name="g(r)",
                    x=r,
                    xlabel="r",
                    xunits="Å",
                    y=res["g"].tolist(),
                    ylabel="g(r)",
                    yunits="",
                    color="#4dbd74",
                )
                y2_axis = plot.add_axis(
                    "y", label="n(r)", anchor=x_axis, overlaying="y", side="right"
                )
            else:
                y2_axis = plot.add_axis("y", label="n(r)", anchor=x_axis)
                x_axis.anchor = y2_axis
            plot.add_trace(
                x_axis=x_axis,
                y_axis=y2_axis,
                name="n(r)",
                x=r,
                xlabel="r",
                xunits="Å",
                y=res["n"].tolist(),
                ylabel="n(r)",
                yunits="",
                color="black",
            )
            figure.grid_plots("rdf")
            name = f"{a}-{b}".replace("(", "_").replace(")", "").replace(" ", "_")
            name = "".join(c for c in name if c.isalnum() or c in "_-")
            path = Path(self.directory) / f"rdf_{name}.graph"
            figure.dump(path)
            if "html" in self.options and self.options["html"]:
                figure.template = "line.html_template"
                figure.dump(path.with_suffix(".html"))

    def _store(self, configuration, results):
        """Store the coordination numbers and pairing as properties (JSON)."""
        properties = configuration.properties
        values = {
            "coordination numbers": {
                k: v.get("CN") for k, v in results["coordination numbers"].items()
            },
        }
        if "ion pairing" in results:
            values["ion pairing"] = {
                k: v.get("fractions") for k, v in results["ion pairing"].items()
            }
        for key, value in values.items():
            name = f"{key}#Liquid Structure#{self.model}"
            if not properties.exists(name):
                properties.add(
                    name,
                    "json",
                    description=f"The {key} of the liquid ({self.model})",
                    noerror=True,
                )
            properties.put(name, value)

    def _find_model(self):
        """The model (forcefield or model chemistry) of the MD, for properties."""
        try:
            if self.variable_exists("_model_chemistry"):
                return self.get_variable("_model_chemistry")["level"]
            ff = self.get_variable("_forcefield")
            if ff == "OpenKIM":
                return "OpenKIM/" + self.get_variable("_OpenKIM_Potential")
            if ff == "PyTorch":
                return "PyTorch/" + self.get_variable("_pytorch_model")
            return ff.current_forcefield
        except Exception:
            return "unknown"

    # ------------------------------------------------------------------
    # What the MD must provide
    # ------------------------------------------------------------------
    def md_requirements(self, P=None):
        """What the analysis needs from the MD: positions every ps or so (every
        0.1 ps for the residence times). The structure is a static property, so
        the dynamics need not be undisturbed."""
        if P is None:
            interval = self._static("frame interval")
            residence = self._static("residence times")
            shells = self._static("solvation shells")
        else:
            interval = P["frame interval"]
            residence = P["residence times"]
            shells = P["solvation shells"]
        if interval is None or isinstance(interval, str):
            interval = Q_(1.0, "ps")
        elif not hasattr(interval, "m_as"):
            interval = Q_(float(interval), "ps")
        why = "the structure"
        if yes(residence) and yes(shells):
            if interval > RESIDENCE_INTERVAL:
                interval = RESIDENCE_INTERVAL
                why = "the residence times in the solvation shells"
        return Requirements(
            needs=[Need("all/position", max_time=interval, why=why)],
            periodic=False,
            dynamics=False,
        )

    def _static(self, key):
        """A parameter's value without evaluating variables, or None."""
        parameter = self.parameters[key]
        if parameter.is_expr:
            return None
        return parameter.get()

    def check_flowchart(self, preceding):
        """Check, before the job runs, that the MD will sample what is needed and
        that the pairs and cutoffs make sense."""
        problems = []
        for key, parse in (("cutoffs", parse_cutoffs), ("pair list", parse_pairs)):
            value = self._static(key)
            if isinstance(value, str):
                try:
                    parse(value)
                except AnalysisError as e:
                    problems.append(f"'{self.title}': {e}")
        source = self._static("md source")
        if source is not None:
            problems.extend(
                check_md_flowchart(self, preceding, self.md_requirements(), source)
            )
        return problems

    def check_configuration(self, configuration):
        """Called by an MD step before it runs MD that this step will analyze.

        Any system can be analyzed; without a cell, only with cutoffs given. A
        liquid in a periodic cell should be large enough.

        Returns
        -------
        [str]
            Warnings about the structure, e.g. that it is small.
        """
        if configuration.periodicity != 3:
            return []
        return [f"'{self.title}': {w}" for w in size_warnings(configuration)]


def _shortest_width(dimensions):
    """The shortest perpendicular width of a cell given as (a, b, c, alpha, beta,
    gamma)."""
    from MDAnalysis.lib.mdamath import triclinic_vectors

    a, b, c = triclinic_vectors(dimensions)
    volume = abs(np.dot(a, np.cross(b, c)))
    return min(
        volume / np.linalg.norm(np.cross(b, c)),
        volume / np.linalg.norm(np.cross(c, a)),
        volume / np.linalg.norm(np.cross(a, b)),
    )
