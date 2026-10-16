.. _user-guide:

**********
User Guide
**********

The Liquid Structure step analyzes the structure of a liquid or solution from the
positions of the atoms in molecular dynamics:

* **radial distribution functions** g(r) of pairs of sites, with their first peak and
  minimum, and **coordination numbers**;
* the **pairing of ions**: each cation is free, in a solvent-separated pair (SSIP), in a
  contact pair (CIP) or in an aggregate (AGG);
* the **solvation shells** of solutes: the number of each solvent in the first shell,
  the most common shells (e.g. Na(5 H2O)), the fraction of solutes paired with each
  solvent, and the residence times of the solvents in the shells.

Every value has an error, from independent runs or from blocks of a single run, and
the g(r) are graphed with their running coordination numbers.

Where the MD comes from
=======================

Like the Viscosity and Diffusivity steps, the step either runs MD itself or analyzes
MD that ran earlier:

* **run the MD in the sub-flowchart**: the step runs the MD in its own flowchart a
  number of times, as independent runs, and averages them;
* **the last MD run** earlier in the flowchart, or **the last MD run of each loop
  iteration**, as independent runs;
* **files**: the ``md.h5md`` files of earlier runs, including other jobs'
  (``job://<n>/<path>``).

The MD must sample the atoms' positions (in LAMMPS, "Sample the atomic positions" on
the trajectory tab), at least every picosecond by default ("Longest interval between
frames"); the residence times need them about every 0.1 ps. When the step runs its own
MD it turns the sampling on; for MD earlier in the flowchart, checking the flowchart
reports what is missing before the job runs.

Species and sites
=================

The molecules are named by their formulas, in Hill order: ``H2O``, ``C3H4O3`` (EC),
``Li``, ``BF4``. A *site* is

* an element of a species, ``O(H2O)``;
* an element in every species, ``O``;
* every atom of a species, ``BF4``;
* or an MDAnalysis selection, ``select: name O and resid 1``.

Pairs of sites are joined by ``-`` (or ``' - '`` when a site contains a hyphen) and
separated by commas: ``O(H2O)-O(H2O), Na-O(H2O), Na-Cl``. By default the step
analyzes every pair of elements.

Radial distribution functions and coordination numbers
======================================================

For each pair, g(r) is accumulated over the frames up to "Range" (at most half the
shortest width of the cell) in bins of "Bin width", leaving out pairs of atoms in the
same molecule. The table gives the first peak and its height, the first minimum, and
the coordination number: the number of B sites within the cutoff of an A site,
averaged over the A sites and frames.

The cutoff is the first minimum unless a **fixed cutoff** is given:
``Na-O(H2O): 3.2, Na-Cl: 3.5``. When the first minimum lies far beyond the peak -- as
when a solvent-shared ion pair adds a peak at about 5 Å -- the step suggests a fixed
cutoff.

Ion pairing
===========

The cations and anions are found from the charges (or, without charges, the monatomic
alkali and alkaline-earth ions are cations and the halides anions). For each
cation-anion pair, each cation in each frame is

* **free:** no anion within r2;
* **SSIP:** no anion within r1, at least one within r2;
* **CIP:** exactly one anion within r1, in contact with no other cation;
* **AGG:** more than one anion within r1, or one also in contact with another cation.

An anion is within a cutoff when any of its heavy atoms is. r1 and r2 are the first
and second minima of the g(r) of the cation and the anion's heavy atoms, or are given
as two fixed cutoffs: ``Na-Cl: 3.5 5.6``. Only distances are used, so the ion
pairing works for clusters too, with the cutoffs given.

Solvation shells
================

The solutes are the monatomic ions by default, or the species given. A solvent
molecule is in a solute's shell when any of its heavy atoms is within the cutoff, the
first minimum of the g(r) of the solute and the solvent's heavy atoms (or a fixed
cutoff, ``Li-C3H4O3: 2.8``). The step reports, for each solvent, the coordination
number and the fraction of solutes with at least one in their shell, and lists the
most common shells with their fractions.

The **residence time** of a solvent is the time at which the autocorrelation of the
occupancy of the shells falls to 1/e. It needs frames about every 0.1 ps; the step
says so when they are further apart.

Systems without a periodic cell
===============================

A cluster has no density, so there are no g(r): the step reports coordination numbers
and ion pairing only for the pairs given fixed cutoffs.

Results
=======

The results -- the frames and runs analyzed, the coordination numbers, the ion
pairing and the solvation shells, each a table of values with their errors -- can be
put in variables and tables, and with "Store as properties" the coordination numbers
and ion pairing are stored with the configuration. The graphs are
``rdf_<pair>.graph`` in the step's directory.

The step uses `MDAnalysis <https://www.mdanalysis.org>`_ to read the trajectory and
find the distances.

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
