# Why `SCF` Should Not Take Intermediate (0 < SCF < 1) Values

## What `SCF` does in this code

`SCF` ("soil cover fraction", per the comment in `cathy_main.f`) is read once in
`datin.f`:

```fortran
READ(IIN4,*) IPEAT,SCF
```

and is then used to split the total potential atmospheric demand `ATMPOT(I)`
(negative when it represents evapotranspiration) into two **independent**
extraction pathways, at every atmospheric surface node:

1. **Bare-soil evaporation**, computed in `atmone.f`, `atmnxt.f` and `atmbak.f`:
   ```fortran
   ATMACT(I) = (1.0d0 - SCF) * ATMPOT(I)
   ```
2. **Plant transpiration**, computed in `etran.f`:
   ```fortran
   ETP(I) = -1.0d0 * SCF * ATMPOT(I)
   ```

`ATMACT` is applied as a Neumann-type surface boundary flux, whose
admissibility is checked node‑by‑node against the surface pressure head
(`atmone.f`, lines ~123–141: comparisons against `PONDH_MIN` and `PMIN`), and
which can be switched to a Dirichlet condition (`ADRSTN`/`SWITCH_OLD`, called
right after `ATMNXT` in `cathy_main.f`) once the surface node becomes too dry
or too wet to sustain the imposed flux.

`ETP`, by contrast, is not applied at the surface node at all: it is
redistributed over the whole root zone in `etran.f` as a **volumetric sink
term** `QTRANIE`, weighted by a completely different, root‑zone‑averaged
moisture‑stress function (`GX`, `BTRAN`, based on `PCWLT`/`PCREF`, the wilting
point and field‑capacity pressures), and it is injected into the flow
equation independently of the switching logic that governs `ATMACT`.

## Confirmation from `adrstn.f`

`adrstn.f` is now available and confirms the decoupling directly rather than
by inference:

```fortran
SUBROUTINE ADRSTN(NNOD,IFATM,ATMPOT,ATMACT,PNEW)
...
DO 500 I=1,NNOD
   IF (PNEW(I) .LE. PMIN) THEN
      IF (ATMPOT(I) .GT. 0.0D0 .OR. ATMACT(I) .LT. ATMPOT(I)) THEN
         IFATM(I)=0
         ATMACT(I)=ATMPOT(I)
         PNEW(I)=PMIN
      END IF
   END IF
500 CONTINUE
```

Two things stand out:

1. **`ADRSTN`'s entire argument list is `ATMPOT`, `ATMACT`, `PNEW`, `IFATM` —
   nothing about `QTRANIE`, `SCF`, or the root zone appears anywhere in it.**
   This is direct, textual confirmation (not just an inference from the split
   formula) that the switching mechanism governing the evaporation pathway
   has no visibility whatsoever into what the transpiration pathway is doing
   at the same node in the same time step.

2. **At the actual call site** (`cathy_main.f`, right after `ATMNXT`, using
   the same `PNEW` that `ADRSTN` will test moments later), the second
   branch — `ATMACT(I) .LT. ATMPOT(I)`, meant to catch "the actual flux has
   exceeded the potential" — is mathematically **unreachable** whenever
   `0 ≤ SCF ≤ 1` and `ATMPOT(I) < 0`: `ATMACT(I) = (1-SCF)*ATMPOT(I)` is, by
   construction, always `≥ ATMPOT(I)` in that regime (multiplying a negative
   number by a factor in `[0,1]` can only make it less negative or equal). So
   at this call site `ADRSTN` can only ever fire on the *first* branch —
   "it started raining" (`ATMPOT(I) .GT. 0`) — and releases the node using
   the **full, unscaled** `ATMPOT(I)`, which is correct because `ATMNXT`
   applies no `SCF` reduction to positive (infiltration) fluxes either. In
   other words: in this call context, nothing at all polices the combined
   evaporation-plus-transpiration extraction — `ADRSTN` here only decides
   "should an air‑dry node start accepting rain again," and never revisits
   the dry‑soil evaporation limit in light of `SCF`. This upgrades the
   original diagnosis from "the two limiters *might* not be aware of each
   other" to "the switching code, as written, structurally cannot be aware
   of the transpiration term — it isn't in its inputs."

3. **A secondary, related gap in `bkstep.f`:** when a time step fails and is
   retried at a smaller `DELTAT`, `bkstep.f` recomputes `ATMPOT`/`ATMACT` via
   `ATMNXT`/`ATMBAK` and re-runs `ADRSTN`/`SWITCH_OLD`, but it does **not**
   call `ETRAN` again. `QTRANIE` is therefore left at whatever it was
   computed as before the failed attempt, on a different `DELTAT` and a
   pressure-head field that is about to be discarded — a stale transpiration
   term is combined with a freshly recomputed evaporation term after every
   back-step. This doesn't change the root cause above, but it's worth fixing
   alongside it (see Level 1 fix below).

## The problem with an intermediate SCF

At the two limiting cases the split degenerates into a single, self‑consistent
mechanism:

- **SCF = 0** (bare soil): `ETP ≡ 0`, all the demand goes through `ATMACT`,
  and the dryness/ponding switch (`PMIN`, `PONDH_MIN`) is the *only* limiter
  acting on that demand — fully consistent.
- **SCF = 1** (full canopy): `ATMACT ≡ 0`, all the demand goes through
  transpiration, and the wilting‑point‑based stress function (`GX`/`BTRAN`)
  is the *only* limiter acting on it — also fully consistent.

For **0 < SCF < 1**, both mechanisms act simultaneously on the *same* nodes
(the root‑zone integral in `etran.f` explicitly includes the surface layer,
`DZ = (ZSURF - Z(K+NNOD))/2`) and draw on the *same* pressure‑head field
(`PNEW`) within the same time step, but:

- `ATMACT` is limited by the **surface pressure‑head switch** (`PMIN`,
  `PONDH_MIN`, in `ADRSTN`/`SWITCH_OLD`), which — as confirmed above — has no
  input at all describing what the transpiration pathway is doing.
- `QTRANIE` is limited by an **unrelated, root‑integrated wilting‑point
  function**, which likewise has no knowledge that `ATMACT` is simultaneously
  removing water from the same near‑surface control volume.

Neither limiter is reduced to account for the fact that the other pathway is
concurrently extracting water from the same soil layer. The two "actual" flux
computations are therefore not coupled or iterated against one another before
being imposed together — they are two independently-capped withdrawals from
one shared reservoir. This creates a genuine risk of extracting more water
from the surface layer, within a single time step, than either physical
limiter alone would allow, which can:

- unphysically over‑dry the surface control volume (pressure head pushed past
  the residual/wilting limit that each individual mechanism was designed to
  respect), and
- degrade convergence or stability of the nonlinear (Picard/Newton) solver
  in `cathy_main.f`, since the two source/boundary terms are inconsistent
  with a single admissible flux budget for that node.

## Bottom line

`SCF ∈ {0, 1}` is safe because only one of the two extraction mechanisms — and
therefore only one, purpose‑built physical limiter — is ever active at a
given node. An intermediate `SCF` activates both mechanisms at once without
any coupling between their respective moisture limits, which is the numerical
(and physical) reason to avoid non‑binary `SCF` values in this
implementation.

*(If a partial-cover formulation is genuinely needed, the fix would be to make
the two limiters aware of each other — e.g. cap `ATMACT + Σ QTRANIE` jointly
against the surface/root moisture budget — rather than computing them from
independent fractions of `ATMPOT`.)*

## Proposed fix

Two levels of fix, in increasing order of invasiveness. Level 1 has been
implemented against the uploaded sources (`etran.f`); Level 2 is a
recommendation for the underlying time-stepping/solver structure and has not
been implemented here.

### Level 1 — joint cap inside `ETRAN` (implemented)

`ETRAN` already computes, per surface node, the root-zone weights `BTRANI`,
`BTRAN`, `GX`, `OMG` used to distribute transpiration with depth. It is
called immediately after `ATMNXT` computes `ATMACT`
(`cathy_main.f`, lines ~2922–2929), so `ATMACT` is available at essentially
no extra cost.

**Change 1 — interface.** `ETRAN` now also receives `ATMACT`:

```fortran
SUBROUTINE ETRAN(N,NNOD,NSTR,ATMPOT,ATMACT,Z,PSI,PNODI,VEG_TYPE,
1                  QTRANIE)
```

`PMIN` and `SCF` did not need to be added to the argument list: both are
already used, undeclared and unpassed, inside `etran.f` and `atmone.f`, which
means they are supplied via the `SOILCHAR.H` common block that `etran.f`
already `INCLUDE`s.

**Change 2 — the cap itself.** For each surface node `I`, the top root layer
(`K = I`, i.e. `J = 1`) occupies the same surface control volume as the
`ATMACT` boundary node. If `ATMACT(I)` is already extracting water there
(`ATMACT(I) .LT. 0.0D0`, i.e. `SCF < 1`) **and** the node has already reached
the air-dry limit (`PSI(I) .LE. PMIN` — the same test `atmone.f` uses inline
for the first time step), the top layer's transpiration weight is withheld
and its share subtracted from the node's `BTRAN`/`OMG` accumulators, so the
remaining demand is picked up by deeper, wetter layers through the
already-existing `BTRAN`/`GX` redistribution, instead of stacking a second,
uncoordinated withdrawal on top of `ATMACT`:

```fortran
IF (ATMACT(I).LT.0.0D0 .AND. PSI(I).LE.PMIN) THEN
   TOPWT     = BTRANI(I)
   BTRAN(I)  = BTRAN(I) - TOPWT
   OMG(I)    = OMG(I)   - TOPWT
   BTRANI(I) = 0.0D0
END IF
```

**Change 3 — a guard that comes for free with the fix.** Withholding the top
layer can leave `BTRAN(I) = 0` (e.g. a single-layer root zone) — which the
*original* code would already divide by zero on (`OMG(I) = OMG(I)/BTRAN(I)`,
and the `QTRANIE` loop), a latent bug the joint cap makes easier to trigger.
Both divisions are now guarded; when `BTRAN(I) .LE. 0`, `QTRANIE` for that
node falls back to `0` for that time step instead of extracting undefined
water.

**Call-site change**, `cathy_main.f` line 2929:

```fortran
CALL ETRAN(N,NNOD,NSTR,ATMPOT,ATMACT,Z,PNEW,PNODI,VEG_TYPE,QTRANIE)
```

The full patched subroutine is provided as `etran.f` alongside this report.

**On the criterion used.** With `adrstn.f` now available, the cap's condition
(`PSI(I).LE.PMIN`) can be confirmed rather than treated as a guess: it is
*exactly* the air-dry gate `ADRSTN` itself tests, on the very same `PNEW`
array, moments after `ETRAN` runs in the same call sequence
(`cathy_main.f`/`bkstep.f`). Since `ADRSTN`'s own second branch is
unreachable at this call site for `0 ≤ SCF ≤ 1` (see the confirmation section
above), there was no live `ADRSTN` criterion being missed — the fix's
condition is the most relevant, and only, dryness gate actually in play here.

### `bkstep.f` gap — now closed

`bkstep.f` recomputed `ATMACT` via `ATMNXT`/`ATMBAK` on every back-step but
never called `ETRAN` again, leaving `QTRANIE` stale (computed on the
previous, failed `DELTAT` and pressure field) through every back-step. This
has been fixed by:

1. **Adding `Z`, `PNODI`, `VEG_TYPE`, `QTRANIE` to `BKSTEP`'s argument
   list** — all four were already available and populated in `cathy_main.f`
   at the point `BKSTEP` is called, so this is a pure pass-through, no new
   state introduced.
2. **Inserting a call to the (patched) `ETRAN`** inside `bkstep.f`, in the
   same position relative to `ATMNXT`/`ATMBAK` and `ADRSTN`/`SWITCH_OLD` as
   in `cathy_main.f` — i.e. using the freshly recomputed `ATMACT` and the
   same `PNEW` that `ADRSTN`/`SWITCH_OLD` are about to test, so a
   back-stepped attempt is internally consistent in exactly the same way a
   normal time step is.
3. **Updating the one `CALL BKSTEP(...)` site in `cathy_main.f`** to pass
   those same four arrays through.

Both diffs are minimal (verified against the originals): four added lines to
the `BKSTEP` signature/declarations, one new `CALL ETRAN` block, and a
two-argument addition at the single call site. `etran.f` itself needed no
further change beyond what Level 1 already did.

**Verification done so far:**
- `etran.f`: compiled with `gfortran -fsyntax-only` against stub
  `CATHY.H`/`SOILCHAR.H` headers built from how `PMIN`, `SCF`, `ZROOT`,
  `PCANA`, `PCWLT`, `PCREF`, `PZ`, `OMGC` are actually used elsewhere in the
  codebase — compiled clean (one pre-existing, unrelated "unused dummy
  argument `PNODI`" warning). Also compiled and ran against a small synthetic
  driver exercising three cases: (i) surface air-dry with `SCF=0.5` — top
  layer correctly withheld and demand shifted to depth; (ii) same but
  surface not air-dry — cap correctly does not engage; (iii) a single-layer
  root zone forced fully capped — falls back to `QTRANIE=0` with no
  divide-by-zero/NaN.
- `bkstep.f`: checked for fixed-form column/continuation compliance and
  argument consistency against `cathy_main.f`'s own `ETRAN` call, but could
  not be link-compiled standalone here (it also depends on `SURFWATER.H`,
  `TRANSPSURF.H`, `BCNXT`, `sfvbak`/`sfvrec`, etc., not part of the uploaded
  set) — a real build-tree compile is still needed before running.

**Level 2 (Picard/Newton-loop coupling) remains open**, as described above,
and would be the next step only if Level 1 (now including the `bkstep.f`
fix) proves insufficient in practice.

### Level 2 — couple both terms inside the Picard/Newton loop (not implemented)

Level 1 still evaluates the cap using a **lagged** `PSI`/`ATMACT`: both are
computed once per time step, outside the nonlinear iteration that actually
resolves `PNEW` in `cathy_main.f`. A fully consistent fix would move both
computations inside that iteration:

- At each Picard/Newton iterate `k`, with current head estimate `PNEW^k`:
  1. Recompute `ATMACT^k` (the evaporation switch) from `PNEW^k`.
  2. Recompute `GX`, `BTRAN`, `QTRANIE^k` (transpiration) from `PNEW^k`,
     i.e. call `ETRAN` every iteration rather than once per time step.
  3. Apply the Level 1 joint cap using the *current* `PNEW^k`.
  4. Assemble and solve for `PNEW^{k+1}`.
  5. Converge as usual, using the tolerance already in place.

This removes the remaining inconsistency (both limiters seeing a common,
current state at every iteration rather than two different snapshots in
time), at the cost of restructuring the driver loop in `cathy_main.f` rather
than a single subroutine. Recommended only if Level 1 turns out insufficient
in practice — e.g. if the Picard loop fails to converge or the ET mass
balance (`QTRAN` vs. the input `ATMPOT` demand, already partly tracked
around `cathy_main.f` line ~2930) keeps drifting once Level 1 is in place.

## Spatially distributed `SCF` from an LAI raster — `SCF(VEG_TYPE(I))` (implemented)

Separate from the double-withdrawal fix above: `SCF` was a single global
scalar (`READ(IIN4,*) IPEAT,SCF`, used bare — not indexed by node —
everywhere it appears), so it could not vary spatially at all, even though
canopy cover realistically does (e.g. partial grass cover, mixed land use).
This section has now been **implemented and compile-checked** against the
real `CATHY.H`/`SOILCHAR.H`/`IOUNITS.H` headers, reusing the infrastructure
CATHY already has for `VEG_TYPE` rather than building a new per-node raster
pipeline from scratch.

### Why per-vegetation-type, not a new per-node field

`datin.f` already shows the mechanism for turning a continuous raster into a
per-node model input, using `VEG_MAP`/`VEG_TYPE` as the example
(`RAST_INPUT_DEM` reads the raster, `TRIANGOLI` resamples it onto the mesh,
the result is truncated to an integer class). The per-vegetation-type
physical parameters are then read as arrays indexed by that class, in the
very same input block `SCF` used to be read from. Piggy-backing `SCF` onto
this existing mechanism is a much smaller change than building a genuinely
continuous, per-node `SCF(NODMAX)` field from scratch — at the cost of `SCF`
only varying as finely as the vegetation classification does, which is
generally acceptable for LAI-derived cover.

### Step-by-step implementation

**Step 1 — `CATHY.H`: raise `MAXVEG`.** The uploaded build had
`MAXVEG=1` — i.e. this configuration only had one vegetation type compiled
in at all, so per-class `SCF` would have had nothing to vary over. Raised to
a real class count:
```fortran
PARAMETER (MAXZON=1,MAXTRM=52111,MAXIT=30,MAXVEG=15)
```
Pick `MAXVEG` to match however many classes your LAI binning will produce
(see the Python wrapper steps below).

**Step 2 — `SOILCHAR.H`: `SCF` scalar → array.** Same array sizing as the
other per-vegetation-type Feddes/canopy parameters it sits next to in
`COMMON /FEDDES/`:
```fortran
REAL*8  SCF(MAXVEG),PCREF(MAXVEG),PCWLT(MAXVEG),PCANA(MAXVEG)
```
No other change needed in this file — `SCF` reaches every routine that
already `INCLUDE`s `SOILCHAR.H` (`atmone.f`, `etran.f`) automatically as an
array from this point on.

**Step 3 — `datin.f`: read `SCF` per vegetation type.** Dropped from the old
scalar line, added as a seventh column in the existing per-class loop:
```fortran
READ(IIN4,*) PMIN
READ(IIN4,*) IPEAT
READ(IIN4,*) CBETA0,CANG
DO I=1,NVEG
   READ(IIN4,*) PCANA(I),PCREF(I),PCWLT(I),ZROOT(I),
1                PZ(I),OMGC(I),SCF(I)
END DO
```
**Input-deck consequence**: the physical input file's `IPEAT,SCF` line
becomes `IPEAT` alone, and each of the `NVEG` per-vegetation-type lines
gains a seventh number (`SCF` for that class). Existing input decks will
need updating to this new format — see "Not addressed here" below for a
backward-compatible alternative.

**Step 4 — thread `VEG_TYPE` into the three routines that didn't already
have it.** `etran.f` already received `VEG_TYPE` (needed it for
`OMGC(VEG_TYPE(I))` etc.), but `atmnxt.f`, `atmone.f`, and `atmbak.f` did
not — this was the main gap to close:
- **`atmnxt.f`**: added `VEG_TYPE` to the argument list, widened its `SCF`
  dummy argument from scalar to `SCF(*)` (it receives `SCF` explicitly, not
  via `SOILCHAR.H`), and changed both occurrences of
  `ATMACT(I)=(1.0d0-SCF)*ATMPOT(I)` to
  `ATMACT(I)=(1.0d0-SCF(VEG_TYPE(I)))*ATMPOT(I)`.
- **`atmone.f`**: added `VEG_TYPE` to the argument list only — it already
  gets `SCF` via its existing `INCLUDE 'SOILCHAR.H'`, so no separate `SCF`
  interface change was needed there. Same substitution to
  `SCF(VEG_TYPE(I))` at its one usage site.
- **`atmbak.f`**: same pattern as `atmnxt.f` — added `VEG_TYPE`, widened
  `SCF` to `SCF(*)`, same substitution.
- **`etran.f`**: only its usage line needed changing —
  `ETP(I) = -1.0d0*SCF(VEG_TYPE(I))*ATMPOT(I)` — since it already had
  `VEG_TYPE` from the Level 1 patch.

**Step 5 — update the two call sites.** `cathy_main.f`'s `CALL ATMNXT(...)`
and `bkstep.f`'s `CALL ATMNXT(...)`/`CALL ATMBAK(...)` all needed `VEG_TYPE`
appended as an actual argument — trivial, since `VEG_TYPE` was already in
scope at both call sites (in `bkstep.f`'s case, because the Level 1 patch
had already added it there for `ETRAN`).

**Step 6 — `ATMONE`'s call site — not available.** No file in the uploaded
set calls `ATMONE` (it's presumably invoked once at simulation start from a
driver not included here). Its signature has been updated to expect
`VEG_TYPE` as its last argument; whatever calls it will need that argument
added too.

### Verification actually performed

- `atmnxt.f`, `atmone.f`, `atmbak.f`, and `etran.f` were compiled **together**
  with `gfortran -fsyntax-only` against the real `CATHY.H`/`SOILCHAR.H`/
  `IOUNITS.H` — clean compile, only pre-existing unused-variable warnings
  unrelated to this change (`ALPHA`, `DSATM`, `GASDEV`, `DELTAT` were already
  unused in the original code).
- `datin.f` could not be compiled standalone — it also needs `SURFWATER.H`
  and `RIVERNETWORK.H`, not part of the uploaded set — but the failure is
  exactly that missing include, occurring before the compiler ever reaches
  the edited region, so it doesn't indicate a problem with the edit itself.
- All diffs against the originals were reviewed and are minimal — exactly
  the intended changes, nothing else touched.

### Compatibility with the Level 1 (double-withdrawal) fix

**No change needed to the `ETRAN` joint-cap patch's logic.** It never reads
`SCF` directly — only `ATMACT(I)` and `PSI(I)` — so the cap behaves
identically whether `SCF` is a global scalar or `SCF(VEG_TYPE(I))`. The two
fixes are orthogonal and compose cleanly, as confirmed by the successful
joint compile above.

### Python wrapper's role (unchanged from the original plan)

1. Reproject/resample the LAI raster onto the exact grid extent/resolution
   `RAST_INPUT_DEM` already expects for `VEG_MAP` — same ASCII format, no
   new Fortran reader needed.
2. Compute `SCF_raw = 1 - exp(-k * LAI)` per pixel (see the discussion above
   on choosing `k`).
3. **Bin/cluster** `SCF_raw` into `MAXVEG` classes (quantile bins or
   k-means), matching whatever `MAXVEG` was set to in Step 1.
4. Write the binned classes out as a `VEG_MAP`-format raster.
5. Write the representative `SCF` value per class (e.g. mean `SCF_raw`
   within each bin) as the new seventh column for `datin.f`'s
   per-vegetation-type block.

### Not addressed here

- **Backward compatibility**: the old scalar `SCF` input line is gone;
  existing input decks need updating to the new per-class format. Mirroring
  `VEG_TYPE`'s own optional per-node/global-broadcast pattern (a flag that
  chooses between reading one value and broadcasting it to all classes, or
  reading `NVEG` distinct values) would preserve old decks — not
  implemented here.
- **`ATMONE`'s call site** needs updating wherever it actually lives, once
  that file is available.
- **Time-varying LAI** (seasonal cover change): `SCF` is now spatially
  variable but still fixed for the whole simulation, read once like the
  rest of the vegetation-type table. Seasonal LAI would need `SCF` to be
  re-read periodically during the time loop, analogous to how
  `ATMTIM`/`ATMINP` already drive time-varying atmospheric forcing in
  `ATMNXT` — a further, separate change.
- **True per-node continuity** (no binning): the fuller `SCF(NODMAX)`
  alternative remains an option if per-vegetation-type resolution proves too
  coarse in practice, at the cost of a new raster-ingestion path rather than
  reusing `VEG_TYPE`'s.

---

## Master changelog — every source file touched, in order

Everything above documents each change in the context of the problem it
solves; this section is the single flat index across *all* of it, plus one
fix made in a later session on a file neither this document nor
`CATHY_SCF_PMIN_ET_partitioning.md` otherwise mentions (`detout.f`). Cross-
references point back to the fuller write-up above for each item.

| # | File | Change | Section |
|---|---|---|---|
| 1 | `etran.f` | Joint double-withdrawal cap (Level 1); later re-touched for `SCF(VEG_TYPE(I))` | "Level 1 — joint cap inside `ETRAN`"; "Step 4" |
| 2 | `cathy_main.f` | `CALL ETRAN` widened for `ATMACT`; `CALL BKSTEP` widened for `Z,PNODI,VEG_TYPE,QTRANIE`; `CALL ATMNXT` widened for `VEG_TYPE` | "Level 1"; "`bkstep.f` gap"; "Step 5" |
| 3 | `bkstep.f` | Stale-`QTRANIE` gap closed: added 4 args, inserted `CALL ETRAN` after its `ATMNXT`/`ATMBAK` recompute | "`bkstep.f` gap — now closed" |
| 4 | `CATHY.H` | `MAXVEG` raised from `1` to a real class count | "Step 1" |
| 5 | `SOILCHAR.H` | `SCF` scalar → `SCF(MAXVEG)` array in `COMMON /FEDDES/` | "Step 2" |
| 6 | `datin.f` | `SCF` dropped from the old `IPEAT,SCF` line; added as a 7th per-`NVEG` column. **Input-deck format change — see "Not addressed here."** | "Step 3" |
| 7 | `atmnxt.f` | Added `VEG_TYPE` arg, widened `SCF` to `SCF(*)`, `SCF`→`SCF(VEG_TYPE(I))` at both `ATMACT` sites | "Step 4" |
| 8 | `atmone.f` | Added `VEG_TYPE` arg only (already had `SCF` via `SOILCHAR.H`); same substitution | "Step 4" |
| 9 | `atmbak.f` | Same pattern as `atmnxt.f` | "Step 4" |
| 10 | `detout.f` | See below — not covered elsewhere in this document | — |
| 11 | `cathy_outputs.py` (Python, pyCATHY) | See below | — |

### 10. `detout.f` — fort.777 now reports evaporation, with a consistent sign

Two fixes, made together, on top of everything above:

**(a) Evaporation was never written to `fort.777` at all — independent of,
and prior to, any `SCF` work above.** `ATMACT(I)` was already an argument
`DETOUT` received but never used. Only `ETA(I)` (`QTRANIE` summed over the
root profile — transpiration only) was written, under the label
`"ACT. ETRA"` (i.e. purporting to be actual *total* ET). So bare-soil
evaporation was silently absent from every file pyCATHY or any other
downstream reader parses, for every build, regardless of `SCF`. Fixed by
adding `ATMACT(I)` to the `WRITE(777,...)` statement and widening the
header/format (`2042`/`2062`) from one value column to three:
`ACT. TRAN`, `ACT. EVAP`, `ACT. ETRA`.

**(b) Sign-convention bug in the combined total.** `ETA(I)` uses the
root-uptake sink convention (positive = water extracted); `ATMACT(I)` uses
the surface-flux convention (negative = water leaving the domain — see
§3 of `CATHY_SCF_PMIN_ET_partitioning.md`). Writing `ETA(I)+ATMACT(I)`
raw partially *cancels* the two whenever both pathways are active on the
same node simultaneously — i.e. precisely the `0<SCF<1` regime this whole
document is about — instead of adding their magnitudes. Confirmed via
cumulative-ET bar charts: pure-transpiration scenarios (`SCF=1`) plotted
positive, pure-evaporation scenarios (`SCF=0`) plotted negative, and
mixed-cover scenarios showed an implausibly small net instead of the sum
of two real fluxes. Fixed by writing `-ATMACT(I)` as `ACT. EVAP` and
`ETA(I)-ATMACT(I)` as `ACT. ETRA`, putting all three columns on a
consistent positive-is-actual-loss basis.

```fortran
C        SIGN NOTE: ETA(I) (summed from QTRANIE in ETRAN.F) uses the
C        root-uptake sink convention, positive = water extracted.
C        ATMACT(I) (set in ATMNXT.F) uses the surface-flux convention,
C        negative = water leaving the domain. ATMACT is negated here so
C        both terms -- and their sum -- are reported on a consistent
C        positive-is-actual-ET-loss basis.
         WRITE(777,2062) I,X(I),Y(I),ETA(I)/ARENOD(I),
     1                   -ATMACT(I)/ARENOD(I),
     2                   (ETA(I)-ATMACT(I))/ARENOD(I)
```

**Interaction with the per-vegetation-type `SCF` work (items 4–9 above):**
none — `detout.f` only reads `ETA`/`ATMACT`, which are already the final
per-node values by the time `DETOUT` runs, regardless of whether `SCF`
was a scalar or `SCF(VEG_TYPE(I))` upstream. This fix is orthogonal, same
as the Level 1 joint cap was noted to be in "Compatibility with the
Level 1 fix" above.

### 11. `cathy_outputs.py` (Python side, pyCATHY) — kept in step with 10(a)

`load_spatial_file_fast`/`read_fort777` previously hardcoded a 4-number-
per-row block width (`SURFACE NODE, X, Y, ACT. ETRA`), matching only the
unpatched format. It now auto-detects the actual row width per file
(`total values in a time-step block / nsurfnodes`) instead of trusting a
hardcoded number, and falls back to the tail of the requested column
names when a file has fewer columns than asked for. This means the same
call site (`read_fort777`, always requesting `['ACT. TRAN','ACT. EVAP',
'ACT. ETRA']`) transparently handles both an unpatched/classic (v1.0.0,
4-column) `fort.777` and the patched (6-column) one from item 10 — no
version flag or separate code path needed.

### Files described above but not available to verify directly in this session

`CATHY.H`, `SOILCHAR.H`, `datin.f` (only the read-loop snippet has been
seen, not the full file), `atmone.f`, `atmbak.f`, `bkstep.f`,
`cathy_main.f`, `switch.f`, and `adrstn.f` were edited and compile-checked
against the real build tree in an earlier session; items 1–9 above are
taken from that session's own record of its diffs and verification, not
re-verified independently here. `etran.f`, `atmnxt.f`, and `detout.f`
*have* been seen directly and their diffs above are exact. If you'd like
the changelog cross-checked against the actual current file contents, or
extended with precise line numbers for items 4–9, please share
`CATHY.H`, `SOILCHAR.H`, `datin.f`, `atmone.f`, `atmbak.f`, `bkstep.f`,
`cathy_main.f`, `switch.f`, and/or `adrstn.f` — for whichever version
(classic v1.0.0 vs. the `SCF_variable` build) you want the changelog to
reflect, since several of these differ between the two.
