=======
History
=======
2026.10.17 -- Analyze MD run before the step; labelled errors for every table
    * The step analyzes MD run before it in the flowchart -- the runs of a Loop around
      the MD, as independent replicas ("the last MD run of each loop iteration", the
      default) -- or md.h5md files; it no longer runs its own MD in a sub-flowchart.
      Older flowcharts that stored its settings for that can still be read; the check
      before the job runs says what replaces them.
    * The check before the job runs also warns when the analysis is inside the loop
      holding its MD, when the loop changes the temperature or pressure, and when the
      thermostat does not suit the analysis, with fixes; the step warns before the MD
      when the system is small.
    * Errors are the 95% confidence interval, Student's t for the number of runs or
      blocks, and the output says so; stored errors are standard errors.
    * The state of the runs -- temperature, pressure, volume and density with their
      errors, and the equilibrium density where an NPT run before them gives the bulk
      modulus -- is reported and stored. The values of each run or block are in
      step.out.
    * The peaks and minima of g(r), coordination numbers, ion pairing, solvation shells
      and residence times have 95% intervals from the runs (or blocks); residence times
      are computed per run.
    * Bugfix: the description claimed residence times, ion pairing or shells even when
      they were turned off.


2026.10.16 -- First release: the structure of liquids and solutions
    * Radial distribution functions of pairs of sites, with errors, their first peaks
      and minima, and coordination numbers, with fixed cutoffs where the minima
      mislead.
    * The pairing of ions: each cation free, in a solvent-separated (SSIP) or contact
      (CIP) pair, or in an aggregate (AGG).
    * The solvation shells of solutes: the coordination number of each solvent, the
      most common shells, pairing, and residence times.
    * Runs its own MD or analyzes the MD earlier in the flowchart, and checks before
      the job runs that the MD samples the positions often enough.
