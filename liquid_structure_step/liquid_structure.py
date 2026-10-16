# -*- coding: utf-8 -*-

"""Non-graphical part of the Liquid Structure step in a SEAMM flowchart

The structure of liquids and solutions from the positions sampled in MD: partial
radial distribution functions and coordination numbers, the classification of ion
pairs (free, solvent-separated, contact, aggregated), and the solvation shells of
solutes (coordination numbers per solvent, the shells' composition, pairing and
residence times). The analysis itself is in :mod:`liquid_structure_step.analysis`.
"""

import importlib.resources
import logging
import math
from pathlib import Path
import sys
import traceback
import warnings

import numpy as np
from seamm_md_analysis import open_md
from seamm_md_analysis.requirements import Need, Requirements
from seamm_md_analysis.step_support import (
    check_md_flowchart,
    md_inputs,
    read_files,
    request_sampling,
    run_own_md,
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
    coordination_number,
    find_cutoff,
    parse_cutoffs,
    parse_pairs,
    peak_and_minimum,
    second_minimum,
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
job = printing.getPrinter()
printer = printing.getPrinter("Liquid Structure")

# Add this module's properties to the standard properties
path = importlib.resources.files("liquid_structure_step") / "data"
csv_file = path / "properties.csv"
if path.exists():
    molsystem.add_properties_from_file(csv_file)

#: The interval between frames that the residence times need
RESIDENCE_INTERVAL = Q_(0.1, "ps")


def fmt(value, err=None, digits=3):
    """A value and its error as text."""
    if value is None:
        return "", ""
    if err is None or not math.isfinite(err):
        return f"{value:.{digits}f}", ""
    return f"{value:.{digits}f}", f"{err:.{digits}f}"


class LiquidStructure(seamm.Node):
    """
    The non-graphical part of a Liquid Structure step in a flowchart.

    Attributes
    ----------
    subflowchart : seamm.Flowchart
        The MD to run, when this step runs its own MD.

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
        namespace="org.molssi.seamm",
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
        namespace : str
            The namespace for the plug-ins of the subflowchart
        extension: None
            Not yet implemented
        logger : Logger = logger
            The logger to use and pass to parent classes
        """
        logger.debug(f"Creating Liquid Structure {self}")
        self.subflowchart = seamm.Flowchart(
            parent=self, name="Liquid Structure", namespace=namespace
        )

        super().__init__(
            flowchart=flowchart,
            title="Liquid Structure",
            extension=extension,
            module=__name__,
            logger=logger,
        )

        self._metadata = liquid_structure_step.metadata
        self.parameters = liquid_structure_step.LiquidStructureParameters()
        self._file_handler = None
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

    def create_parser(self):
        """Setup the command-line / config file parser"""
        parser_name = "liquid-structure-step"
        parser = seamm_util.getParser()

        self.subflowchart._parser = self.flowchart._parser

        parser_exists = parser.exists(parser_name)
        super().create_parser(name=parser_name)
        if not parser_exists:
            parser.add_argument(
                parser_name,
                "--html",
                action="store_true",
                help="whether to write out html files for graphs, etc.",
            )

        self.subflowchart.reset_visited()
        node = self.subflowchart.get_node("1").next()
        while node is not None:
            node = node.create_parser()

        return self.next()

    def description_text(self, P=None, short=False):
        """Create the text description of what this step will do.

        Parameters
        ----------
        P: dict
            An optional dictionary of the current values of the control
            parameters.
        short : bool = False
            Leave out the description of the sub-flowchart.

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
        if P["ion pairing"]:
            extras.append("the pairing of the ions")
        if P["solvation shells"]:
            extras.append(
                "the solvation shells of "
                + (
                    P["solute list"]
                    if P["solutes"] == "the solutes given"
                    else "the monatomic ions"
                )
                + (", with residence times" if P["residence times"] else "")
            )
        if extras:
            what += ", and " + " and ".join(extras)

        source = P.get("md source", run_own_md)
        if source == run_own_md:
            text = f"Analyze the structure: {what}, from {P['nruns']} runs of the MD "
            text += "below.\n\n"
        else:
            text = f"Analyze the structure: {what}, from "
            if source == read_files:
                text += f"the MD in {P['md files']}.\n\n"
            else:
                text += f"{source} earlier in the flowchart.\n\n"
            short = True

        if not short:
            self.subflowchart.root_directory = self.flowchart.root_directory
            node = self.subflowchart.get_node("1").next()
            while node is not None:
                try:
                    text += __(node.description_text()).__str__()
                except Exception as e:
                    self.logger.critical(
                        f"Error describing the liquid structure flowchart: {e} in "
                        f"{node}"
                    )
                    raise
                text += "\n"
                node = node.next()

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

        printer.important(__(self.description_text(P, short=True), indent=self.indent))

        # Check what is asked for before any MD runs
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

        job_handler = None
        out_handler = None
        job_level = None
        out_level = None
        for handler in job.handlers:
            if (
                isinstance(handler, logging.FileHandler)
                and "job.out" in handler.baseFilename
            ):
                job_handler = handler
                job_level = job_handler.level
                job_handler.setLevel(printing.JOB)
            elif isinstance(handler, logging.StreamHandler):
                out_handler = handler
                out_level = out_handler.level
                out_handler.setLevel(printing.JOB)

        if P["md source"] == run_own_md:
            self._run_own_md(P)
        else:
            self.paths = md_inputs(self, P, name="all/position")
            printer.job(
                __(
                    f"Analyzing {len(self.paths)} MD run"
                    f"{'s' if len(self.paths) > 1 else ''}: "
                    + ", ".join(str(p) for p in self.paths),
                    indent=self.indent + 4 * " ",
                )
            )

        if job_handler is not None:
            job_handler.setLevel(job_level)
        if out_handler is not None:
            out_handler.setLevel(out_level)

        if len(self.paths) == 0:
            raise RuntimeError("No MD run finished, so there is nothing to analyze.")

        _, configuration = self.get_system_configuration()
        self.analyze(P=P, configuration=configuration)

        return next_node

    def _run_own_md(self, P):
        """Run the MD in the sub-flowchart several times, keeping each run's file."""
        changes, missing = request_sampling(self.subflowchart, self.md_requirements(P))
        for change in changes:
            printer.job(__(f"Note: {change}.", indent=self.indent + 4 * " "))
        if len(missing) > 0:
            raise RuntimeError(
                "The MD in the sub-flowchart cannot provide what is needed: "
                + " ".join(missing)
            )

        first_node = self.subflowchart.get_node("1").next()
        if first_node is None:
            raise RuntimeError(
                "The sub-flowchart of the Liquid Structure step is empty."
            )
        node = first_node
        while node is not None:
            node.all_options = self.all_options
            node = node.next()
        self.subflowchart.executor = self.flowchart.executor

        nruns = P["nruns"]
        fmt_run = f"0{len(str(nruns))}"
        for run in range(1, nruns + 1):
            run_dir = Path(self.directory) / f"run_{run:{fmt_run}}"
            run_dir.mkdir(parents=True, exist_ok=True)

            if self._file_handler is not None:
                self._file_handler.close()
                job.removeHandler(self._file_handler)
            path = run_dir / "Run.out"
            path.unlink(missing_ok=True)
            self._file_handler = logging.FileHandler(path)
            self._file_handler.setLevel(printing.NORMAL)
            formatter = logging.Formatter(fmt="{message:s}", style="{")
            self._file_handler.setFormatter(formatter)
            job.addHandler(self._file_handler)

            self.flowchart.reset_visited()
            self.set_subids((*self._id, f"run_{run:{fmt_run}}"))

            node = first_node
            try:
                while node is not None:
                    node = node.run()
            except DeprecationWarning as e:
                printer.normal("\nDeprecation warning: " + str(e))
                traceback.print_exc(file=sys.stderr)
                traceback.print_exc(file=sys.stdout)
            except Exception as e:
                printer.job(f"Caught exception in run {run}: {str(e)}")
                with open(run_dir / "stderr.out", "a") as fd:
                    traceback.print_exc(file=fd)
                if "continue" in P["errors"]:
                    continue
                elif "exit" in P["errors"]:
                    break
                else:
                    raise
            else:
                self.paths.append(self._read_run(run, run_dir))

        if self._file_handler is not None:
            self._file_handler.close()
            job.removeHandler(self._file_handler)
            self._file_handler = None

    def _read_run(self, run, run_dir):
        """The md.h5md of one run of the sub-flowchart with the positions."""
        paths = []
        for path in sorted(run_dir.glob("**/md.h5md")):
            with open_md(path) as md:
                if md.has_particles("position"):
                    paths.append(path)
        if len(paths) == 0:
            raise RuntimeError(f"There are no positions for run {run} in {run_dir}.")
        return paths[-1]

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
        n_blocks = nruns if nruns > 1 else P["number of blocks"]
        if nruns == 1 and P["md source"] == run_own_md:
            n_blocks = 1

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

        text = f"Analyzed {n_frames} frames from {nruns} run{'s' if nruns > 1 else ''}"
        if dt is not None:
            text += f", every {dt:.3g} ps"
        if n_blocks > 1:
            text += (
                f"; the errors are the standard errors over the {nruns} runs"
                if nruns > 1
                else f"; the errors are from {n_blocks} blocks of frames"
            )
        printer.normal(__(text + ".", indent=self.indent + 4 * " "))
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
        rdf_results = {key: rdf.result(n_blocks) for key, rdf in rdfs.items()}
        # The ion pairing and shells read the frames again, with the cutoffs
        frames = lambda: self._frames(configuration, periodic)  # noqa: E731
        results["coordination numbers"] = self._coordination(
            pairs, rdf_results, cutoffs, periodic
        )
        if ion_pairs:
            results["ion pairing"] = self._ion_pairing(
                ion_pairs, rdf_results, cutoffs, species, frames, n_blocks
            )
        if solutes:
            results["solvation shells"] = self._shells(
                P,
                solutes,
                rdf_results,
                cutoffs,
                species,
                frames,
                n_blocks,
                dt,
            )

        self._plot(pairs, rdf_results)

        # The results, for variables, tables and properties
        self.model = self._find_model()
        data = {"n_frames": n_frames, "nruns": nruns}
        data.update(results)
        if P["store"] and configuration is not None:
            self._store(configuration, results)
        self.store_results(configuration=configuration, data=data)

    def _coordination(self, pairs, rdf_results, cutoffs, periodic):
        """The peaks, minima and coordination numbers of the pairs asked for."""
        table = {
            "Pair": [],
            "Peak (Å)": [],
            "g(peak)": [],
            "Minimum (Å)": [],
            "Cutoff (Å)": [],
            "CN": [],
            "±": [],
        }
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
                entry.update(
                    {
                        "peak": pm["peak"],
                        "peak height": pm["peak_height"],
                        "minimum": pm["minimum"],
                    }
                )
            cutoff = fixed[0] if fixed else pm["minimum"]
            if cutoff is not None:
                n, err = coordination_number(res, cutoff)
                entry.update({"cutoff": cutoff, "CN": n, "CN,stderr": err})
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
            elif not periodic:
                notes.append(
                    f"{label}: no cutoff was given, so no coordination number."
                )
            results[label] = entry
            table["Pair"].append(label)
            table["Peak (Å)"].append(fmt(pm["peak"], digits=2)[0])
            table["g(peak)"].append(fmt(pm["peak_height"], digits=2)[0])
            table["Minimum (Å)"].append(fmt(pm["minimum"], digits=2)[0])
            table["Cutoff (Å)"].append(
                fmt(cutoff, digits=2)[0] + (" (given)" if fixed else "")
            )
            v, e = fmt(entry.get("CN"), entry.get("CN,stderr"), digits=2)
            table["CN"].append(v)
            table["±"].append(e)

        self._print_table(table, "Radial distribution functions")
        for note in notes:
            printer.normal(__(note, indent=8 * " "))
        if notes:
            printer.normal("")
        return results

    def _frames(self, configuration, periodic, first=False):
        """The frames of every run: positions (Å), box and volume (Å^3), the box
        and volume None without a cell. Sets ``self._dt``, the interval between
        frames in ps."""
        dts = []
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
                for ts in u.trajectory:
                    box = ts.dimensions if periodic else None
                    volume = ts.volume if periodic else None
                    yield ts.positions, box, volume
                    if first:
                        break
            if first:
                break
        self._dt = float(np.median(dts)) if dts else None

    def _ion_pairing(self, ion_pairs, rdf_results, cutoffs, species, frames, n_blocks):
        """Classify the cations of each cation-anion pair in every frame."""
        table = {"Pair": [], "r1 (Å)": [], "r2 (Å)": []}
        for c in PAIR_CLASSES:
            table[c] = []
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
            r = acc.result(n_blocks)
            results[label] = {"r1": r1, "r2": r2, **r}
            table["Pair"].append(label)
            table["r1 (Å)"].append(f"{r1:.2f}")
            table["r2 (Å)"].append(f"{r2:.2f}" if r2 != r1 else "")
            for c in PAIR_CLASSES:
                v, e = fmt(
                    r["fractions"][c],
                    None if r["errors"] is None else r["errors"][c],
                )
                table[c].append(v + (f" ± {e}" if e else ""))
        if table["Pair"]:
            self._print_table(table, "Ion pairing (fraction of cations)")
            printer.normal(
                __(
                    "free: no anion within r2; SSIP: none within r1, one or more "
                    "within r2; CIP: one anion within r1, touching no other cation; "
                    "AGG: an anion within r1 shared with another cation, or more "
                    "than one.",
                    indent=8 * " ",
                )
            )
            for note in notes:
                printer.normal(__(note, indent=8 * " "))
            printer.normal("")
        return results

    def _shells(self, P, solutes, rdf_results, cutoffs, species, frames, n_blocks, dt):
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
                continue
            acc = ShellAccumulator(solute_atoms, solute_ids, solvents)
            for x, box, _ in frames():
                acc.add_frame(x, box=box)
            residence_dt = dt if P["residence times"] else None
            r = acc.result(n_blocks, dt=residence_dt)
            r["cutoffs"] = {name: v[2] for name, v in solvents.items()}
            results[solute] = r

            table = {"Solvent": [], "Cutoff (Å)": [], "CN": [], "±": [], "Paired": []}
            if r["residence"] is not None:
                table["Residence (ps)"] = []
            for name in solvents:
                table["Solvent"].append(name)
                table["Cutoff (Å)"].append(f"{r['cutoffs'][name]:.2f}")
                v, e = fmt(
                    r["coordination"][name],
                    (
                        None
                        if r["coordination_err"] is None
                        else r["coordination_err"][name]
                    ),
                    digits=2,
                )
                table["CN"].append(v)
                table["±"].append(e)
                table["Paired"].append(f"{r['pairing'][name]:.3f}")
                if r["residence"] is not None:
                    t = r["residence"][name]
                    table["Residence (ps)"].append("" if t is None else f"{t:.2f}")
            self._print_table(table, f"The solvation shell of {solute}")
            lines = []
            for composition, fraction in r["speciation"][:5]:
                shell = ", ".join(
                    f"{count} {name}" if count > 1 else name
                    for name, count in composition.items()
                )
                lines.append(f"{fraction:6.3f}  {solute}({shell or 'none'})")
            printer.normal(
                __(
                    "The most common shells:\n" + "\n".join(lines),
                    indent=8 * " ",
                    wrap=False,
                    dedent=False,
                )
            )
            if skipped:
                printer.normal(
                    __(
                        f"No cutoff for {', '.join(skipped)} around {solute}, so "
                        "they are left out.",
                        indent=8 * " ",
                    )
                )
            if P["residence times"] and dt is not None and dt > 0.15:
                printer.normal(
                    __(
                        f"The frames are {dt:.2g} ps apart; the residence times "
                        f"need frames every {RESIDENCE_INTERVAL:~P} or so.",
                        indent=8 * " ",
                    )
                )
            printer.normal("")
        return results

    def _print_table(self, table, title):
        text = tabulate(
            table, headers="keys", tablefmt="simple_outline", disable_numparse=True
        )
        length = len(text.splitlines()[0])
        printer.normal(
            __(
                "\n" + title.center(length) + "\n" + text + "\n",
                indent=8 * " ",
                wrap=False,
                dedent=False,
            )
        )

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
                    upper = (res["g"] + 2 * res["g_err"]).tolist()
                    lower = (res["g"] - 2 * res["g_err"]).tolist()
                    plot.add_trace(
                        x_axis=x_axis,
                        y_axis=y_axis,
                        name="95%",
                        x=r + r[::-1],
                        xlabel="r",
                        xunits="Å",
                        y=upper + lower[::-1],
                        ylabel="95%",
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
        0.1 ps for the residence times)."""
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
        if residence in (True, "yes") and shells in (True, "yes"):
            if interval > RESIDENCE_INTERVAL:
                interval = RESIDENCE_INTERVAL
                why = "the residence times in the solvation shells"
        return Requirements(
            needs=[Need("all/position", max_time=interval, why=why)],
            periodic=False,
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
        Any system can be analyzed; without a cell, only with cutoffs given."""
        return

    def set_id(self, node_id=()):
        """Sequentially number the subnodes"""
        if self.visited:
            return None
        else:
            self.visited = True
            self._id = node_id
            self.set_subids(self._id)
            return self.next()

    def set_subids(self, node_id=()):
        """Set the ids of the nodes in the subflowchart"""
        node = self.subflowchart.get_node("1").next()
        n = 1
        while node is not None:
            node = node.set_id((*node_id, str(n)))
            n += 1


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
