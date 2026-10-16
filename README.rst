==============================
SEAMM Liquid Structure Plug-in
==============================

.. image:: https://img.shields.io/github/issues-pr-raw/molssi-seamm/liquid_structure_step
   :target: https://github.com/molssi-seamm/liquid_structure_step/pulls
   :alt: GitHub pull requests

.. image:: https://github.com/molssi-seamm/liquid_structure_step/workflows/CI/badge.svg
   :target: https://github.com/molssi-seamm/liquid_structure_step/actions
   :alt: Build Status

.. image:: https://codecov.io/gh/molssi-seamm/liquid_structure_step/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/molssi-seamm/liquid_structure_step
   :alt: Code Coverage

.. image:: https://github.com/molssi-seamm/liquid_structure_step/workflows/CodeQL/badge.svg
   :target: https://github.com/molssi-seamm/liquid_structure_step/security/code-scanning
   :alt: Code Quality

.. image:: https://github.com/molssi-seamm/liquid_structure_step/workflows/Release/badge.svg
   :target: https://molssi-seamm.github.io/liquid_structure_step/index.html
   :alt: Documentation Status

.. image:: https://img.shields.io/pypi/v/liquid_structure_step.svg
   :target: https://pypi.python.org/pypi/liquid_structure_step
   :alt: PyPi VERSION

A SEAMM plug-in for the structure of liquids and solutions from molecular dynamics:
radial distribution functions, coordination numbers, ion pairing and solvation shells.

* Free software: BSD-3-Clause
* Documentation: https://molssi-seamm.github.io/liquid_structure_step/index.html
* Code: https://github.com/molssi-seamm/liquid_structure_step

Features
--------

* Partial radial distribution functions with errors, their first peaks and minima, and
  coordination numbers, with fixed cutoffs where the minima mislead.
* The pairing of ions: free, solvent-separated (SSIP), contact (CIP) and aggregated
  (AGG) pairs, also for clusters.
* Solvation shells: coordination numbers per solvent, the most common shells, pairing,
  and residence times.
* Runs its own MD or analyzes MD earlier in the flowchart, from SEAMM's H5MD files.

Acknowledgements
----------------

This package was created with the `molssi-seamm/cookiecutter-seamm-plugin`_ tool, which
is based on the excellent Cookiecutter_.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`molssi-seamm/cookiecutter-seamm-plugin`: https://github.com/molssi-seamm/cookiecutter-seamm-plugin

Developed by the Molecular Sciences Software Institute (MolSSI_),
which receives funding from the `National Science Foundation`_ under
award CHE-2136142.

.. _MolSSI: https://molssi.org
.. _`National Science Foundation`: https://www.nsf.gov
