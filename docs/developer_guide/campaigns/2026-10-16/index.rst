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

The real-data check, 2026-10-16, on mlff's NaCl(aq) run 12197 (one-channel LES
committee, 512 H2O + 9 NaCl, metal units; frames at 1 ps, t = 50-200 ps; on ChemAI),
against mlff's ``analyze_nacl_md.py`` and ``cip_fixed.py``:

- at mlff's minima: Na-O CN 5.50 ± 0.02 (5.52), Cl-O 6.51 ± 0.03 (6.58), Cl-H
  6.01 ± 0.02 (6.05); the first peaks identical (Na-O 2.41 Å, g 5.64; Cl-O 3.21 Å,
  g 3.15; Cl-H 2.21 Å);
- the fraction of Na with a Cl within 3.0 / 3.5 / 4.2 Å: 0.024 / 0.037 / 0.073
  (0.02 / 0.04 / 0.07), and the Na-Cl g(r) at 0.1 Å bins identical;
- it found a weakness, fixed: with nine ion pairs the Na-Cl g(r) is noisy, and the
  automatic first minimum fell on a noise dip at 2.77 Å. The g(r) is now smoothed
  over 0.2 Å and a crossing of 1 must last 0.2 Å; the minimum is then 3.31 Å. Its
  second minimum is not clear, and the step now says so (SSIP not separated)
  rather than reporting none silently.
