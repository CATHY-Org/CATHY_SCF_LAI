# Validation study: float SCF vs SCF in {0,1}, and the role of PMIN

Adds `validation/scf_pmin_study/` only; **no change to `src/`, `input/` or `prepro/`**.

## What it contains
* `README.md` - study write-up (physics read from the Fortran, 17-run matrix, 12-day dry-down, heterogeneous cover, release-rule test, classic-vs-SCF comparison, critical assessment).
* `PLAN.md` - how each step was done.
* `pycathy_module/scf_pmin_study.py` - the study re-expressed through the pyCATHY `CATHY()` class (quick and drydown suites executed; the others defined, not re-run).
* `harness/` - the direct runner (builds both versions, runs the matrix) that produced the numbers in the README; `results*/` - summaries, tables, figures (`.npz` omitted, regenerable).

## Main findings (see README for numbers and caveats)
1. SCF build == classic at SCF = 0 and 1 (E, T within 0.07 %; identical iteration counts in most runs).
2. Float SCF closes the water budget (<= 0.04 % of V0 once root uptake is added back to `mbeconv`) with zero back-steps in every run. A convergence problem of classic with float SCF was **not** reproduced; a difference shows up only as an accounting error (phantom transpiration) for a very thin root zone with PMIN above PCWLT.
3. PMIN must be <= PCWLT when SCF > 0: otherwise the pinned surface node injects water into a stressed root zone (1.5 mm over 12 d at SCF = 1).
4. Binary 0|1 maps and a float SCF are not equivalent at equal mean cover (19-29 % ET difference): float SCF is a sub-grid closure.

## Possible separate issues spotted (not fixed here)
* `src/Makefile` requests `conf_fac.f`; the file is `Conf_Fac.f` (fails on case-sensitive file systems).
* `input/sfbc` ends at 604800 s: runs longer than 7 days stop with "End of file" in `sfvnxt.f`.
* `MAXPRT = 9` in `CATHY.H`: more print times stop the run silently with exit code 0.
* `etran.f` prints per-node debug lines (~44 MB per simulated hour).

## Limits
One closed 10 x 10 m column, uniform soil, 6-12 days, model-vs-model only; no steep slopes, ponding, large meshes, restart tests or field data.
