=======
History
=======
2026.10.20.1 -- Bugfix: an empty structure is refused at once
    * Without a structure (no atoms in the current configuration), the step now stops
      at once with a clear message, before reading the MD files.

2026.10.20 -- About ten times faster for large liquids
    * Finding the pairs of atoms within the range of the radial distribution
      functions, and for the ion pairing and solvation shells, now uses a periodic
      k-d tree for rectangular cells, about ten times faster than before. For a
      9,436-atom electrolyte with every pair of elements to 12 Å, a frame took
      about 5 s, so 20,000 frames took more than a day; now about 2 hours.
      Triclinic cells are handled as before. Distances are now computed in double
      precision, so a pair exactly at the edge of a bin can rarely fall in the
      next one.

2026.10.19.1 -- Scaled ion charges
    * The cations and anions of MD with scaled charges (LAMMPS's 'Scale the charges of
      ions by') are found from their formal charges, for any factor.

2026.10.19 -- Bugfix: long runs were very slow and used too much memory
    * Every frame's g(r) was kept, about 10 kB per pair of sites per frame, so long runs
      ran out of memory: four 3 ns runs of 1 m NaCl(aq) with positions every 0.1 ps
      (120,000 frames) took 10.6 hours, mostly swapping. Only the sums over each run or
      block are kept now, with the same results. The same analysis now takes 17
      minutes in 1 GB (an hour with the radial distribution functions from every 0.1
      ps frame, as before, which reproduces the earlier results exactly).
    * The residence times took a time growing as the cube of the length of the runs;
      they are now found with FFTs, with the same results, in seconds.
    * The radial distribution functions use frames about every 'Longest interval
      between frames' even when residence times are computed; only the ion pairing and
      solvation shells use the frames every 0.1 ps. The ion pairing and all the shells
      are found in one pass over the frames, and each frame is binned faster.

2026.10.18 -- Bugfix: frames more often than needed were all analyzed
    * The step analyzed every frame the MD wrote, so positions written every 0.1 ps were
      all used even when nothing needed them: a 1000-atom argon run spent half of a
      12-minute job here. Frames are now skipped down to the interval needed: every 0.1
      ps when residence times are computed, which needs solutes, else about every
      'Longest interval between frames'.

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
