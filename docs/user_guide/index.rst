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

The step analyzes MD that ran before it in the flowchart, or the ``md.h5md`` files of
earlier runs ("MD to analyze"):

* **the last MD run of each loop iteration** (the default): for independent runs, put
  the MD in a Loop -- in parallel, with each run starting from the same structure
  ("Discard the structure" in the MD sub-steps) and new velocities -- and the step,
  after the loop, analyzes the run of every iteration together. Without a loop it is
  the one MD run before the step;
* **the last MD run** earlier in the flowchart;
* **files**: the ``md.h5md`` files of earlier runs, including other jobs'
  (``job://<n>/<path>``), each an independent run.

A typical flowchart:

.. code-block:: text

    Model Chemistry (FF@...), Packmol, minimization
    NPT equilibration; Set Cell to its mean volume; NVT equilibration
    Loop over the replicas (parallel)
        LAMMPS: velocities, NVT decorrelation, NVT production
    Liquid Structure

The MD must sample the atoms' positions (in LAMMPS, "Sample the atomic positions" on
the trajectory tab), at least every picosecond by default ("Longest interval between
frames"); the residence times need them about every 0.1 ps. Checking the flowchart
before it runs reports what is missing, with a fix that turns the sampling on, and
also warns when the analysis is inside the loop holding its MD, when the loop changes
the temperature or pressure (the runs are then not replicas of one state), or when the
thermostat does not give canonical fluctuations. The structure is a static property,
so any canonical thermostat will do. Before the MD runs, the step warns if a periodic
system is small for the analysis (fewer than about 200 molecules).

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

Errors
======

The errors are from the independent runs, or, for a single run, from blocks of it
("Blocks for errors from one run", 5 by default), each analyzed on its own. Every
"± (95%)" column in the tables is the 95% confidence interval: Student's t with n − 1
degrees of freedom (2.78 for 5 runs or blocks) times the standard error, and a line
under each table says so. The peak positions and heights and the minima have errors
from the g(r) of each run or block; the coordination numbers, pairing fractions,
shells and residence times likewise. A value with no error (one run in one block) is
shown without ±. The band around g(r) in the graphs is the same 95% interval. The
coordination numbers of each run or block are in the step's own output, ``step.out``.

The state
=========

After the structure the step reports the state the runs sampled, in the table shared
by the analysis steps: temperature, pressure, volume and density with their errors.
For runs at constant volume the pressure is as measured, and the equilibrium density
at the target pressure is estimated when an NPT run before them gives the bulk
modulus. The structure belongs to the state that was run.

Results
=======

The results -- the frames and runs analyzed, the coordination numbers, the ion
pairing and the solvation shells, and the state of the runs -- can be put in variables
and tables. Each value is stored with its standard error (``CN,stderr``,
``fractions,stderr``, ``T,stderr``, ...), not the 95% interval, together with what the
errors are from ("error source": runs or blocks) and how many there are
("n_samples"). With "Store as properties" the coordination numbers and ion pairing
are stored with the configuration. The graphs are ``rdf_<pair>.graph`` in the step's
directory.

The step uses `MDAnalysis <https://www.mdanalysis.org>`_ to read the trajectory and
find the distances.

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
