# -*- coding: utf-8 -*-

"""
liquid_structure_step
A SEAMM plug-in for Liquid Structure
"""

# Bring up the classes so that they appear to be directly in
# the liquid_structure_step package.

from .liquid_structure import LiquidStructure  # noqa: F401, E501
from .liquid_structure_parameters import LiquidStructureParameters  # noqa: F401, E501
from .liquid_structure_step import LiquidStructureStep  # noqa: F401, E501
from .tk_liquid_structure import TkLiquidStructure  # noqa: F401, E501

from .metadata import metadata  # noqa: F401

# Handle versioneer
from ._version import get_versions

__author__ = "Paul Saxe"
__email__ = "psaxe@molssi.org"
versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
