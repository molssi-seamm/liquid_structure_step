# -*- coding: utf-8 -*-
"""
Control parameters for the Liquid Structure step in a SEAMM flowchart
"""

import logging
import seamm
from seamm_md_analysis.step_support import source_parameters

logger = logging.getLogger(__name__)


class LiquidStructureParameters(seamm.Parameters):
    """
    The control parameters for Liquid Structure.

    The MD to analyze is run earlier in the flowchart -- in a Loop for independent
    runs -- or in other jobs (see seamm_md_analysis.step_support.source_parameters).
    """

    # Dropped with the own-MD route (phase 8): ignored in older flowcharts
    retired = ("nruns", "errors")

    parameters = {
        **source_parameters("the structure"),
        "frame interval": {
            "default": 1.0,
            "kind": "float",
            "default_units": "ps",
            "enumeration": tuple(),
            "format_string": ".2f",
            "description": "Longest interval between frames:",
            "help_text": (
                "The positions must be sampled at least this often. The structure "
                "changes slowly, so about a picosecond is enough; the residence "
                "times need frames ten times as often."
            ),
        },
        "pairs": {
            "default": "every pair of elements",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": ("every pair of elements", "the pairs given"),
            "format_string": "s",
            "description": "Radial distribution functions for:",
            "help_text": (
                "The pairs of sites for the radial distribution functions and "
                "coordination numbers."
            ),
        },
        "pair list": {
            "default": "",
            "kind": "string",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "s",
            "description": "Pairs:",
            "help_text": (
                "The pairs of sites, separated by commas, each two sites joined by "
                "'-' or ' - ': e.g. 'O(H2O)-O(H2O), Na-O(H2O), Na-Cl'. A site is an "
                "element of a species, 'O(H2O)'; an element, 'Na'; a species, "
                "'BF4'; or an MDAnalysis selection, 'select: name O'."
            ),
            "applies_when": {"pairs": "the pairs given"},
        },
        "r max": {
            "default": 12.0,
            "kind": "float",
            "default_units": "Å",
            "enumeration": tuple(),
            "format_string": ".1f",
            "description": "Range:",
            "help_text": (
                "The range of the radial distribution functions, at most half the "
                "shortest width of the cell."
            ),
        },
        "bin width": {
            "default": 0.02,
            "kind": "float",
            "default_units": "Å",
            "enumeration": tuple(),
            "format_string": ".3f",
            "description": "Bin width:",
            "help_text": "The width of the bins of the radial distribution functions.",
        },
        "cutoffs": {
            "default": "",
            "kind": "string",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "s",
            "description": "Fixed cutoffs:",
            "help_text": (
                "Cutoffs for coordination numbers and ion pairing instead of the "
                "minima of the radial distribution functions, e.g. "
                "'Na-O(H2O): 3.2, Na-Cl: 3.5' (Å). Two values for a cation-anion "
                "pair give the contact and solvent-separated cutoffs: "
                "'Na-Cl: 3.5 5.6'."
            ),
        },
        "ion pairing": {
            "default": "yes",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "s",
            "description": "Ion pairing:",
            "help_text": (
                "Classify each cation as free, solvent-separated (SSIP), in contact "
                "with one anion (CIP), or in an aggregate (AGG)."
            ),
        },
        "solvation shells": {
            "default": "yes",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "s",
            "description": "Solvation shells:",
            "help_text": (
                "Analyze the solvation shells of the solutes: the coordination "
                "number of each species, the shells' composition, and pairing. "
                "Without a periodic cell, only for the species given fixed cutoffs."
            ),
        },
        "solutes": {
            "default": "the monatomic ions",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": ("the monatomic ions", "the solutes given"),
            "format_string": "s",
            "description": "Solutes:",
            "help_text": "The species whose solvation shells are analyzed.",
            "applies_when": {"solvation shells": "yes"},
        },
        "solute list": {
            "default": "",
            "kind": "string",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "s",
            "description": "Solutes:",
            "help_text": "The solute species, separated by commas, e.g. 'Li, Na'.",
            "applies_when": {"solvation shells": "yes", "solutes": "the solutes given"},
        },
        "residence times": {
            "default": "no",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "s",
            "description": "Residence times:",
            "help_text": (
                "Also the residence times of the solvents in the solvation shells. "
                "They need frames about every 0.1 ps."
            ),
            "applies_when": {"solvation shells": "yes"},
        },
        "store": {
            "default": "no",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "s",
            "description": "Store as properties:",
            "help_text": (
                "Store the coordination numbers and ion-pairing fractions as "
                "properties of the configuration."
            ),
        },
    }

    def __init__(self, defaults={}, data=None):
        """
        Initialize the parameters, by default with the parameters defined above

        Parameters
        ----------
        defaults: dict
            A dictionary of parameters to initialize. The parameters
            above are used first and any given will override/add to them.
        data: dict
            A dictionary of keys and a subdictionary with value and units
            for updating the current, default values.

        Returns
        -------
        None
        """

        logger.debug("LiquidStructureParameters.__init__")

        super().__init__(
            defaults={**LiquidStructureParameters.parameters, **defaults}, data=data
        )
