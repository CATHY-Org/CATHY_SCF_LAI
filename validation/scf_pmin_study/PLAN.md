# Plan: SCF (float) vs SCF ∈ {0,1} (classic), and the role of PMIN, for semi-arid CATHY simulations

**Goal.** Build smoke tests that compare the CATHY build with a *fractional* soil cover fraction
(`SCF_variable`, 0 < SCF < 1, per vegetation type) against the *classic* use (SCF = 0 or 1, one global scalar),
in order to get closer to semi-arid conditions, where CATHY is known to struggle. Understand — from the Fortran —
how `PMIN` (small |PMIN| vs. effectively disabled, e.g. −1e35) controls the evaporation boundary-condition switch,
and when the Feddes-type root-water-uptake stress is (or is not) in play.

Status legend: ✅ done · 🔄 in progress

---

## Step 0 — Set up and verify what is being compared  ✅
1. Clone the pyCATHY wrapper (`BenjMy/pycathy_wrapper`) and read how it selects versions
   (`CATHY(version=...)`, `SCF_variable`, `set_SOIL_defaults`: *"SCF = 1 → Feddes approach, 0 → PMIN approach"*).
2. Get both Fortran code bases:
   * **classic**: `CATHY-Org/CATHY_src_org` at the commit pinned by the wrapper (`ef54607`); `SCF` is one scalar on the `IPEAT SCF` line of `input/soil`.
   * **SCF build**: `CATHY-Org/CATHY_SCF_LAI` (main); `SCF(MAXVEG)` is a per-vegetation-type array on its own trailing line of `input/soil`; `MAXVEG=15`.
3. Diff them: exactly 11 files differ (`CATHY.H, SOILCHAR.H, datin.f, atmnxt.f, atmone.f, atmbak.f, bkstep.f, etran.f, cathy_main.f, inital.f, detout.f`).
   Everything else (solver, `switch.f`, `switch_old.f`, `adrstn.f`) is identical → the PMIN switching logic is the same in both.

## Step 1 — Understand the physics from the Fortran  ✅
Read, in this order, and write the findings in `README.md §2–3`:
| Routine | What to extract |
|---|---|
| `atmnxt.f`, `atmone.f`, `atmbak.f` | `ATMACT = (1−SCF)·ATMPOT` (evaporation share) — only updated while the node is flux-controlled (`IFATM=0`) |
| `etran.f` | `ETP = SCF·|ATMPOT|`, root-profile weights, Feddes stress `GX=min(GX1,GX2)` (PCWLT/PCREF/PCANA), compensation `OMGC`; **PMIN does not enter except in the SCF build's top-layer cap** |
| `switch_old.f`, `adrstn.f`, `switch.f` (called from `conver.f`/`cathy_main.f`/`bkstep.f`) | the `IFATM` state machine: 0 = Neumann, 1 = Dirichlet (pinned at 0 or at PMIN), 2 = ponded; transitions at `PNEW ≤ PMIN` |
| `datin.f`, `cathy_main.f` | input-file layout and which switching routine runs for which `ISIMGR` |
| `detout.f` | what `fort.777` reports (SCF build: 3 columns; classic: transpiration only) |

## Step 2 — Make both builds run reproducibly  ✅
`harness/build.sh` builds `cathy_classic` and `cathy_scf` from scratch copies. Two environment fixes are needed on Linux
and are applied only in the scratch copy: (a) `Makefile` asks for `conf_fac.f` but the file is `Conf_Fac.f`;
(b) the SCF build's `etran.f` prints several lines per surface node per step (~44 MB/hour) — silenced. Neither changes numerics.

## Step 3 — Define a semi-arid test bed  ✅
* Domain: the repo's 10 m × 10 m, 441 surface nodes, 15 layers, 3 m deep (finest layer 6 mm) — uniform properties, so it behaves as 441 identical columns.
* Initial state: uniform dry soil, ψ₀ = −20 m. No-flow bottom/sides: closed water budget (clean mass-balance check).
* Forcing (step function, `IETO=1`): 25 mm rain pulse (day 0–0.25) → 5 mm/d potential ET dry spell → 8 mm pulse (day 4–4.1) → dry to day 6.
  Pulse-driven, evaporation-dominated: the semi-arid regime.
* Vegetation: one class, Feddes PCANA=0, PCREF=−4 m, PCWLT=−150 m, ZROOT=0.5 m, OMGC=1.

## Step 4 — Experiment matrix  ✅
| Axis | Values | Why |
|---|---|---|
| Build | `classic`, `scf` | regression of the new build against the reference |
| SCF | 0, 0.1, 0.3, 1 | 0/1 = classic endpoints; 0.1–0.3 = sparse arid/semi-arid cover (fc = 1 − exp(−k·LAI)) |
| PMIN (m) | −5 (wrapper default), −150 (= PCWLT, file default), −1e4 (≈ physical air-dry), −1e35 (switch disabled) | "small" → "very big" |

17 runs: 13 on the SCF build (SCF {0, 0.3, 1} × 4 PMIN + SCF 0.1 at −150) + 4 on the classic build (SCF 0/1 × PMIN −5/−1e35). SCF 0.1 was dropped outside PMIN −150 to save time (one core).
`harness/study_matrix.py` runs them and stores `results/summary.json` + per-run `.npz` series.

## Step 5 — Metrics  ✅
Cumulative evaporation E (∫ −ATMACT), cumulative transpiration T (∫ QTRAN), E/(E+T), end of "stage 1" (first time evaporation drops below 95 % of its own demand),
E and T at day 3, surface ψ (is it pinned at PMIN? how low does it go?), mass-balance error, back-steps / nonlinear iterations (numerical stress), wall time.

## Step 6 — Smoke tests (`tests/test_scf_pmin_smoke.py`, 1-day runs, several minutes)  ✅ (15/15 pass)
Each test states a physical or code-derived invariant, not a numerical snapshot:
1. **Runs & closes the budget** — both builds exit 0, outputs exist, water-budget residual < 0.1 % of initial storage (`mbeconv` excludes root uptake; add ∫QTRAN dt back).
2. **SCF = 0 ⇒ no transpiration; SCF = 1 ⇒ no evaporation** (PMIN disabled, so the switch cannot interfere).
3. **Partition** — for 0 < SCF < 1 before any limitation: E ≈ (1−SCF)·PET and T ≤ SCF·PET; E + T ≤ PET always (no "double withdrawal").
4. **Regression** — SCF build at SCF = 0 and 1 reproduces the classic build (cumulative E and T within a stated tolerance).
5. **PMIN monotonicity** (SCF = 0) — cumulative E is non-decreasing as PMIN becomes more negative; with finite PMIN the surface head never goes below PMIN.
6. **PMIN-independence of transpiration** (SCF = 1) — PMIN is not in `etran.f`; T must not depend on it (any dependence is reported as a finding, since the switch can still act on the surface node).
7. **Disabled PMIN stress** (SCF < 1, PMIN = −1e35) — record how far the surface head falls and whether back-steps/iterations blow up (reported, with a soft threshold).

## Step 7 — PMIN discussion  ✅
Per version: SCF = 0 (PMIN is the *only* limiter of E), SCF = 1 (Feddes is the *only* limiter of T; PMIN should be inert),
0 < SCF < 1 (both act on the same near-surface node; the SCF build adds a joint cap). Small |PMIN| vs. −1e35. Written up in `README.md §4–5` using the measured numbers.

## Step 8 — Deliverables  ✅
`PLAN.md` (this file), `README.md` (study write-up), `harness/` (build, runner, matrix, plots), `tests/` (pytest smoke tests), `results/` (summary + figures).

## Step 9 — Follow-up: 12-day dry-down  ✅
The 6-day matrix never stressed the root zone. A 12-day dry-down (SCF {0.3, 1} × PMIN {−5, −150}, `harness/study_drydown.py`) shows Feddes stress and
a PMIN-above-PCWLT artefact (phantom infiltration through the pinned surface node). Written up in `README.md §8`; tests in `tests/test_drydown.py` (`-m drydown`).
Not done: targeted test of the `adrstn.f`/`switch_old.f` release rule at SCF=0.3, and `ISIMGR=1`.

## Step 10 — Heterogeneous cover (binary 0|1 vs float SCF)  ✅
`harness/study_extra.py` (block B): clumped patches, 3-cell stripes, cell checkerboard vs uniform float SCF at the nominal and the *realized* cover
(class boundaries are averaged and truncated by `INT()`; the checkerboard collapses to SCF ≈ 0). Results and discussion: `README.md §11`.

## Step 11 — Release rule of the pinned node: explanation and targeted test  ✅
Read `adrstn.f`, `switch_old.f`, `switch.f`; PET staircase for SCF 0/0.3/0.6 from a dry start (block A) and with a shallow water table (block C). `README.md §12`.
Outcome: the inconsistency (comparison with the full `ATMPOT`) exists in the code but the window never opens for PET ≳ 0.6 mm/d at PMIN −150 m.

## Step 12 — Is the "classic + float SCF = numerical problem" claim real?  ✅
`harness/study_critique.py`, `study_critique2.py`: classic vs SCF build at SCF = 0.3 for PMIN −5/−150/−1e35 and thin root zones (worst case for the joint cap).
Outcome: no convergence difference; an *accounting* difference (phantom transpiration, 1.2 % of V0) only with a thin root zone and PMIN above PCWLT. `README.md §10`.

## Step 13 — Why a spatial SCF rather than ZROOT per class; motivation for semi-arid and post-fire  ✅
ZROOT experiment (SCF = 1; 0.006–0.5 m): changes T, never creates E. `README.md §9`.

## Step 14 — Critical assessment of float SCF and a suggested workflow (monthly SCF update)  ✅
`README.md §13–14`. Not tested and flagged: restart continuity between monthly runs, steep/ponding/overland-flow cases, large meshes, field data.

## Step 15 — Replace the pytest smoke tests by a `CATHY()`-based module  ✅ (partly executed)
`pycathy_module/scf_pmin_study.py`: scenarios + invariant checks for all blocks (quick, matrix, drydown, release, hetero, critique), driven through
`pyCATHY.cathy_tools.CATHY` (setters `update_veg_map/update_soil/update_atmbc/update_ic/update_sfbc/update_parm`, `run_processor`).
Executed through the module: `quick` (16/16 checks pass; classic and SCF build identical at the endpoints) and `drydown` (see README §15).
Defined but not re-executed through the module (numbers come from the legacy direct runner, which the module reproduces exactly on every scenario compared): matrix, release, hetero, critique.
pyCATHY pitfalls found while building it are listed in `README.md §15`. The old pytest files are kept for provenance only.

## Known limits (stated up front)
* Closed column with uniform soil: no lateral flow, no groundwater, no runoff; it isolates the ET-partition/PMIN mechanics, not catchment hydrology.
* One soil, one forcing, 6 days: conclusions are about mechanisms and ordering, not calibrated arid-zone predictions.
* SCF is static per vegetation class (no seasonal LAI); `0 < SCF < 1` is flagged in the repo's own notes as having an uncoupled-limiter risk — this study measures it rather than assuming it.
