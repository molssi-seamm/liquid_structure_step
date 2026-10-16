2026-10-16 The first version
============================

The Liquid Structure step is phase 6 of the MD analysis campaign of 2026-10-14, whose
design and plan are on the main SEAMM documentation site (developer guide,
campaigns, 2026-10-14: ``PLAN_md_analysis.rst`` and ``PLAN_phase6_structure.rst``).

Decisions (Paul, 2026-10-16):

#. The name: "Structure" was taken by the step for energy-optimized structures.
#. Comparing g(r) with X-ray data at the experiment's resolution is a follow-on.
#. Ion pairing by the step's own classification with cutoffs, which works for
   clusters too.
#. Residence times included, off by default.
#. The solvation shells are the step's own code: SolvationAnalysis, planned at
   first, is GPL-3, while MDAnalysis has been LGPL-3 since 2.8 (so the step requires
   MDAnalysis 2.8 or later).

Checked on 2026-10-16 with 1 molal NaCl in 216 SPC/E waters (50 ps, frames every
0.1 ps): Na-O CN 5.1, Cl-O CN 7.1, the Na shell 5 or 6 waters in 71 % of the frames,
residence times 11 ps (Na) and 4 ps (Cl).
