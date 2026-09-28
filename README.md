# SCF and PMIN in CATHY's ET Partitioning — Mechanics, Arid-Zone Behavior, and LAI Coupling

> **Status note (added after this doc was written):** §1 below describes
> `SCF` as it originally existed — a single global scalar, read once in
> `datin.f`, with `ATMACT`/`ETP` computed from a bare `SCF` (no `VEG_TYPE`
> indexing). That has since been **implemented** as a per-vegetation-type
> array, `SCF(VEG_TYPE(I))`, replacing every scalar usage shown in §1 —
> see `SCF_intermediate_values_issue.md`, section "Spatially distributed
> `SCF` from an LAI raster", for the full per-file diff. §§2–5 below (the
> `PMIN`/switching mechanics, arid-zone tuning guidance, and the
> Beer–Lambert LAI→`fc` conversion) are unaffected by that change and
> still apply as written — only the "global scalar" framing in §1 and the
> `datin.f`/`atmnxt.f`/`etran.f` code snippets there are now historical.
> Separately, `detout.f`'s `fort.777` output has also since been fixed to
> report evaporation (`ATMACT`) at all (previously silently absent) and
> with a sign consistent with transpiration (`ETA`) — see the same
> document's master changelog for detail. Neither change alters the
> physics described here, only what's spatially resolved and what's
> written to disk.

This note documents how the soil cover fraction (`SCF`) and air-dry pressure head (`PMIN`) parameters control evapotranspiration (ET) partitioning in the CATHY model, based on `atmbak.f`, `atmnxt.f`, `datin.f`, `bkstep.f`, `cathy_main.f`, `etran.f`, and `switch.f`.

---

## 1. The two ET pathways

CATHY splits potential evapotranspiration (`ATMPOT`, negative when it represents an evaporative/transpirative demand) into two physically distinct sinks:

| Pathway | Routine | Mechanism |
|---|---|---|
| **Bare-soil evaporation** | `atmnxt.f` / `atmbak.f` | Applied as a surface (Neumann-type) boundary flux |
| **Transpiration** | `etran.f` | Applied as a distributed root-zone volumetric sink |

`SCF` (soil cover fraction) is the static coefficient that decides how much of the demand goes to each pathway:

```fortran
! atmnxt.f / atmbak.f — bare-soil evaporation branch
ATMACT(I) = (1.0d0 - SCF) * ATMPOT(I)

! etran.f — transpiration branch
ETP(I)    = -1.0d0 * SCF * ATMPOT(I)
```

Since `(1-SCF) + SCF = 1`, the two branches are exactly complementary — no demand is lost or double-counted. `SCF` is read once, globally, as a single scalar:

```fortran
! datin.f
READ(IIN4,*) IPEAT,SCF
```

`SCF = 0` → 100% bare-soil evaporation, 0% transpiration.
`SCF = 1` → 0% bare-soil evaporation (even at saturation — the split is moisture-independent), 100% transpiration.
Intermediate values split the demand proportionally and are numerically safe: both equations are linear in `SCF`, with no division or branch logic keyed on its value.

---

## 2. What happens inside the transpiration branch (`etran.f`)

Once `SCF` hands `ETP(I) = SCF × |ATMPOT(I)|` to `ETRAN`, it is distributed over the root profile and reduced by soil-water stress:

- **Root-density weighting** (`BETA`), integrated over depth to `ZROOT(VEG_TYPE)`:
  ```fortran
  BETA = (1 - DEPTH/ZROOT) * EXP(-PZ * DEPTH/ZROOT)
  ```
- **Feddes-type stress function** `GX = min(GX1, GX2)`:
  - `GX1` — dry-side limit: ramps 0→1 between `PCWLT` (wilting point) and `PCREF` (reference head).
  - `GX2` — wet-side (anoxia) limit: a near step-function that collapses to 0 within 0.001 of `PCANA`.
- **Compensation**: uptake is divided by `MAX(OMG, OMGC)`, where `OMG` is the root-weighted average stress and `OMGC(VEG_TYPE)` is a critical threshold below which the plant is assumed to compensate by drawing more from less-stressed layers.

Important: `PMIN` **never appears in `etran.f`**. The transpiration branch has its own independent dryness limit (`PCWLT`/`PCREF`), fully separate from the bare-soil branch's `PMIN`.

**Numerical caveat**: `OMG(I) = OMG(I)/BTRAN(I)` and the final `QTRANIE` both divide by `BTRAN(I)`. If `BTRAN(I) = 0` (e.g. `ZROOT` too small/degenerate for a `VEG_TYPE`), this produces a divide-by-zero regardless of `SCF`. This risk exists as soon as `SCF > 0` anywhere in the domain — every active `VEG_TYPE` needs a sane, non-degenerate `ZROOT`.

---

## 3. What happens inside the bare-soil branch — `PMIN` and `SWITCH`

The `(1-SCF)` fraction of the demand is imposed as a Neumann (flux-controlled) boundary condition at the surface node — but the soil can't sustain an arbitrarily large suction to keep delivering that flux. `PMIN` is the physical floor ("air-dry" pressure head) at which the code gives up on flux control and pins the head instead.

`IFATM` tracks the surface BC state per node:
- `IFATM = 0` → Neumann (flux prescribed, unsaturated)
- `IFATM = 1` → Dirichlet, pinned at a fixed head (saturated-but-not-ponded, or air-dry)
- `IFATM = 2` → Dirichlet, ponded (`PNEW = PL`, the ponding head)

The relevant evaporation case in `switch.f` (case J, `IFATM=0` with `ATMPOT<0`):

```fortran
IF (PNEW(I) .LE. PMIN) THEN
   IFATM(I) = 1
   PNEW(I)  = PMIN
   OVFLNOD(I) = -PONDNOD(I) * ARENOD(I) / DELTAT
   GO TO 500
END IF
```

Once the computed head `PNEW(I)` would drop to or below `PMIN`, the node is switched to Dirichlet and pinned exactly at `PMIN`. **This is a switch in which quantity is prescribed (head instead of flux) — it does not set the flux to zero.**

Once pinned, the actual flux the node delivers is whatever Darcy's law allows at that fixed head:

```
q = -K(ψ) × dH/dx        (Darcy's law, unsaturated form)
```

with `ψ = PMIN`. As the profile just below the surface keeps drying toward equilibrium with `PMIN`, the gradient `dH/dx` shrinks and `K(ψ)` itself collapses (often by orders of magnitude as soil dries), so:

- The flux drops sharply right at the moment of switching (no longer tracking the full potential demand).
- It is **not instantly zero** — there is typically still some gradient from wetter soil below.
- It **decays toward zero** over subsequent time steps as the near-surface profile dries into equilibrium with `PMIN` — governed by conductivity/gradient, not a hard cutoff.

There is also a corresponding wet-dry re-transition once the node is at `IFATM=1` (case H, `switch.f` lines 351–362): if the profile re-wets enough that a Neumann flux at `ATMPOT` could be sustained without violating `PMIN`, the node switches back to `IFATM=0`. This back-and-forth flipping near the threshold is a known source of iteration stiffness in Richards-equation solvers — see §5.

---

## 4. Behavior for semi-arid ecosystems with sparse vegetation

### 4.1 SCF
- Sparse cover means most exposed ground is bare soil → `SCF` should reflect the *actual fractional vegetation cover*, typically **~0.1–0.3** for arid shrub/steppe systems, not values approaching 0.5+.
- Consequence: most of the ET demand is routed through the bare-soil branch, which matches observed arid-zone behavior — evaporation dominates water loss immediately after rain pulses, transpiration is a smaller, more sustained draw from patchy/deep-rooted vegetation.
- **Structural limitation**: `SCF` is a single global scalar, not spatially variable and not indexed by `VEG_TYPE` the way `ZROOT`, `PCANA`, `PCREF`, `PCWLT`, `PZ`, `OMGC` are. It cannot represent within-domain heterogeneity (e.g. bare interspaces vs. vegetated clumps) on its own — see §6 for how to work around this with LAI data.

### 4.2 PMIN
- With `SCF` low, a large share of ET *is* going through the `PMIN`-governed bare-soil branch, so this parameter genuinely matters (unlike an `SCF=1` scenario where it becomes moot).
- Arid soils routinely reach much drier states than a mild threshold like -15 (units consistent with your input file, typically m or cm of head). Real air-dry matric potentials are texture-dependent and often **hundreds to thousands of units more negative** than mesic-soil defaults — clayey soils hold water more tightly and reach far more negative air-dry heads than sandy soils.
- Setting `PMIN` too mild (i.e. not negative enough) triggers the Neumann→Dirichlet switch prematurely, artificially truncating the energy-limited (stage-1) evaporation period and understating how long a real arid soil surface can keep evaporating near its potential rate before becoming supply-limited.
- **Do not simply disable `PMIN` (e.g. `-1e35`) in a low-SCF setup** — unlike the `SCF=1` case where the bare-soil branch is already forced to zero, here disabling `PMIN` removes a real physical brake on the surface node and risks numerically unrealistic drying / solver instability.
- Recommended practice: derive `PMIN` from the actual soil water retention curve (its residual/air-dry point) for the dominant surface texture in the domain, rather than reusing a generic or mesic-soil default.
- **Published `PMIN` (`ψmin`) values and their limits.** `PMIN` is in meters of water head and is a calibrated parameter, not a universal constant. In a paired-catchment CATHY study in temperate southwestern Victoria, Australia, Camporese et al. (2013) assigned -10 m (-0.1 MPa) for pasture and -50 m (-0.5 MPa) for a *Eucalyptus globulus* plantation, the only parameters tuned besides the ponding threshold. Camporese et al. (2014) then showed that a single `ψmin` can reproduce catchment-scale actual ET in a water-limited pasture, that lowering it increases ET and reduces runoff and recharge, and that the resulting ET reduction is analogous to Feddes stress. That validation covers **shallow-rooted vegetation**; the authors conclude the approach is limited for **semi-arid environments and deep-rooted vegetation** [verify this wording against the full text of the conclusions, which was not accessible when this note was written]. A surface-node switch responds only to surface moisture, so deep roots that keep transpiring from wetter layers below a dry surface are better represented by the root-uptake pathway in `etran.f` (`ZROOT`, Feddes, `OMGC`) than by `PMIN`. Do not transfer the -10 m / -50 m values to semi-arid or deep-rooted settings without independent evidence.

### 4.3 Companion parameters worth retuning together
- `PCWLT` / `PCREF` / `PCANA` — xerophytic vegetation typically operates at much more negative wilting-point and reference pressures than mesic vegetation; mesic defaults will choke off transpiration too early.
- `OMGC` — arid-adapted plants are known for strong compensatory root water uptake; consider raising `OMGC` relative to mesic defaults so the compensation term (`MAX(OMG,OMGC)`) can actually engage instead of uniformly suppressing uptake across the profile.
- `ZROOT` — matters strongly here since sparse arid vegetation can be either deep-taprooted (drawing on stored moisture) or shallow-and-lateral (capturing pulse infiltration); this directly controls how much of `ETP` from `etran.f` is actually realized versus stress-limited.

---

## 5. Numerical considerations summary

| Risk | Cause | Mitigation |
|---|---|---|
| Divide-by-zero / NaN in `QTRANIE` | `BTRAN(I)=0` from degenerate `ZROOT` for an active `VEG_TYPE` | Ensure every `VEG_TYPE` with `SCF>0` has a physically sane, non-zero `ZROOT` |
| Iteration oscillation | Node flipping between `IFATM=0` and `IFATM=1` near the `PMIN` threshold (`switch.f` cases G/H/I/J) | Most pronounced when `SCF` is low (more flux forced onto the surface node); avoid an overly mild `PMIN` that sits right in the "active" pressure-head range for your climate |
| Understated stage-1 evaporation | `PMIN` too mild for the actual soil texture | Set `PMIN` from the soil's own retention curve, not a generic default |
| Overstated / unphysical drying | `PMIN` disabled (e.g. `-1e35`) while `SCF<1` still routes real flux through the bare-soil branch | Only disable/relax `PMIN` this aggressively when the bare-soil branch carries negligible or zero demand (e.g. `SCF≈1`) |

---

## 6. Recommendations for coupling with LAI raster maps

`SCF` ("soil cover fraction") is conceptually equivalent to **fractional vegetation cover (`fc`)**, which is routinely derived from LAI (or NDVI) raster products in remote sensing / land-surface modeling. The challenge specific to this codebase is that `SCF` is read as a single global scalar (`datin.f`), while LAI is naturally spatially (and often temporally) distributed.

### 6.1 Converting LAI to an equivalent SCF/fc
The standard relationship is a Beer–Lambert canopy-extinction form:

```
fc = 1 - exp(-k * Ω * LAI)
```

- `k` — light extinction coefficient, typically **0.5–0.6** for randomly distributed canopies, but often lower (**~0.3–0.5**) for the clumped, discontinuous canopies typical of arid shrublands/steppe — using a mesic default `k` will overestimate effective cover for sparse arid vegetation.
- `Ω` — clumping index (≤1). Sparse arid vegetation is frequently strongly clumped (isolated shrubs/tussocks rather than a continuous canopy), and ignoring `Ω` (i.e. assuming `Ω=1`) will bias `fc` high relative to true bare-soil fraction. If a clumping-corrected LAI/NDVI product isn't available, treat this as a source of uncertainty rather than omitting it silently.

### 6.2 Reconciling a spatial `fc` raster with a scalar `SCF`
Since the code as provided only supports one global `SCF`, three practical options, in order of increasing fidelity/effort:

1. **Domain-averaged scalar (quick, lowest fidelity)**: compute an area-weighted mean of the `fc` raster over the catchment/subcatchment and use that as the single `SCF`. Acceptable for a first-order water balance but discards spatial heterogeneity — particularly problematic in arid basins where cover is often concentrated in narrow riparian corridors against otherwise bare hillslopes.
2. **Zonal/sub-basin runs**: partition the domain into zones with contrasting `fc` (e.g. wadi-bottom/riparian vs. bare hillslope), and either (a) run each zone with its own `SCF` if your workflow supports zone-specific parameter files, or (b) at minimum verify the mesh/zone structure already used for `VEG_TYPE`, `ZROOT`, etc. and align an `SCF` estimate to each zone consistently with those.
3. **Code-level extension (highest fidelity — implemented, see status note at
   top)**: since `VEG_TYPE`, `ZROOT`, and the Feddes parameters are already read as rasters/per-zone arrays elsewhere in the input chain (e.g. `IIN3` root map), the more correct long-term fix is to extend `datin.f` to read `SCF` the same way — as a per-node or per-`VEG_TYPE` array rather than a single scalar — and propagate that through `atmnxt.f`, `atmbak.f`, and `etran.f` (replacing the scalar `SCF` with `SCF(VEG_TYPE(I))` or `SCF(I)`). This is a real code change, not a parameter-file change, so treat it as a development task rather than a workaround.

### 6.3 Temporal dynamics
Arid-zone LAI can change sharply with rainfall pulses (green-up after rain, rapid senescence in drydown). A single static `SCF` calibrated from a dry-season LAI composite will understate transpiration (and overstate bare-soil evaporation) during green-up windows, and vice versa if calibrated from a wet-season composite. Two pragmatic approaches:
- If the model/workflow supports periodically re-reading `IIN4`-style parameters (analogous to how atmospheric BCs are time-interpolated via `ATMTIM`/`ATMINP` in `atmnxt.f`), update `SCF` on a similar cadence from a time series of LAI composites (monthly is often a reasonable compromise for semi-arid phenology).
- If not, run bounding scenarios using dry-season and wet-season representative `SCF` endmembers to bracket the sensitivity of your ET partitioning to phenological cover changes, rather than relying on a single static value.

### 6.4 Summary recommendation
- Use `fc = 1 - exp(-k·Ω·LAI)` with a clumping-aware, low-to-moderate `k` (~0.3–0.5) as your best estimate of `SCF` for sparse arid canopies.
- Treat the current global-scalar `SCF` as a domain-average approximation; if spatial heterogeneity in cover is a first-order control on your water balance (likely, in a patchy arid landscape), prioritize the code-level extension to a per-zone/per-node `SCF` over trying to force fidelity out of a single scalar.
- Pair any LAI-derived `SCF` update with a consistency check on `ZROOT`/`PCWLT`/`PCREF`/`PCANA`/`OMGC` for the corresponding `VEG_TYPE`, since these interact directly with how much of the `SCF`-allocated transpiration demand actually gets satisfied.

---

## References

- Camporese, M., Dean, J. F., Dresel, P. E., Webb, J., Daly, E. (2013). Hydrological modelling of paired catchments with competing land uses. *20th International Congress on Modelling and Simulation (MODSIM2013)*, Adelaide, Australia, 1819-1825.
- Camporese, M., Daly, E., Dresel, P. E., Webb, J. A. (2014). Simplified modeling of catchment-scale evapotranspiration via boundary condition switching. *Advances in Water Resources*, 69, 95-105. https://doi.org/10.1016/j.advwatres.2014.04.008
