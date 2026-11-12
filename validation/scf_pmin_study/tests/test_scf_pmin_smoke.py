"""Smoke tests: CATHY with fractional SCF vs classic SCF in {0,1}, and the role of PMIN.

Run:  bash harness/build.sh harness && pytest -v tests            (about 5-8 min on one core)
      pytest -v -m stress tests                                  (PMIN disabled, SCF<1: numerical stress)
Each test asserts a physical or code-derived invariant (see README.md section 3), not a numeric snapshot.
"""
import numpy as np, pytest
from conftest import series, cum, surface_psi_min, budget_residual_pct, PET

# ---------------------------------------------------------------- 1. runs and closes the water budget
@pytest.mark.parametrize("version,scf,pmin", [("scf", 0.0, -150.0), ("scf", 0.3, -150.0), ("scf", 1.0, -1e35),
                                              ("classic", 0.0, -150.0), ("classic", 1.0, -1e35)])
def test_runs_and_closes_budget(run, version, scf, pmin):
    spec, r = run(version, scf, pmin)
    assert r["exit"] == 0 and r["completed"], f"{spec.name} did not complete (exit={r['exit']})"
    for k in ("t", "act", "tr", "psi", "sw"): assert r.get(k) is not None, f"missing output {k}"
    assert np.all(np.isfinite(r["act"])) and np.all(np.isfinite(r["tr"])) and np.all(np.isfinite(r["psi"]))
    assert budget_residual_pct(r) < 0.1, "water-budget residual (mbeconv + root uptake) above 0.1 % of initial storage"

# ---------------------------------------------------------------- 2. endpoints of SCF
def test_scf0_has_no_transpiration(run):
    _, r = run("scf", 0.0, -150.0)
    _, _, _, T = series(r)
    assert T.max() < 1e-3, "SCF=0 must route 100 % of the demand through the surface flux (T = 0)"

def test_scf1_has_no_evaporation_when_pmin_disabled(run):
    _, r = run("scf", 1.0, -1e35)
    t, E, _, _ = series(r)
    dry = t > 0.10
    assert E[dry].max() < 0.05, "SCF=1 and PMIN off: ATMACT=(1-SCF)*ATMPOT=0, no surface evaporation expected"
    assert cum(t, E) < 0.05

# ---------------------------------------------------------------- 3. partition for 0 < SCF < 1
def test_partition_fraction_scf_03(run):
    _, r = run("scf", 0.3, -150.0)
    t, E, tt, T = series(r)
    # window right after the pulse: surface still wet (no PMIN limitation), root zone not yet stressed
    w = (t > 0.06) & (t < 0.20); wt = (tt > 0.06) & (tt < 0.20)
    assert abs(E[w].mean() / PET - 0.7) < 0.05, "E should be (1-SCF)*PET = 70 % of demand while unlimited"
    assert 0.5 * 0.3 < T[wt].mean() / PET <= 0.3 * 1.02, "T should be <= SCF*PET (Feddes can only reduce it)"
    # no double withdrawal: total ET never exceeds the potential demand
    tot = np.interp(t, tt, T) + E
    assert tot[t > 0.06].max() <= PET * 1.02, "E + T exceeded potential ET"

# ---------------------------------------------------------------- 4. regression: SCF build reproduces classic at SCF=0 and 1
@pytest.mark.parametrize("scf,pmin,what", [(0.0, -150.0, "E"), (1.0, -1e35, "T")])
def test_scf_build_matches_classic_at_endpoints(run, scf, pmin, what):
    _, a = run("scf", scf, pmin); _, b = run("classic", scf, pmin)
    ta, Ea, tta, Ta = series(a); tb, Eb, ttb, Tb = series(b)
    x, y = (cum(ta, Ea), cum(tb, Eb)) if what == "E" else (cum(tta, Ta), cum(ttb, Tb))
    assert y > 1.0, "reference quantity unexpectedly small"
    assert abs(x - y) / y < 0.02, f"cumulative {what}: scf-build {x:.3f} mm vs classic {y:.3f} mm"

# ---------------------------------------------------------------- 5. PMIN: the switch is the only evaporation limiter at SCF=0
def test_pmin_monotonic_evaporation_scf0(run):
    E = {}
    for p in (-5.0, -150.0, -1e4):
        _, r = run("scf", 0.0, p); t, e, _, _ = series(r); E[p] = cum(t, e)
    assert E[-5.0] <= E[-150.0] * 1.001 <= E[-1e4] * 1.002, f"E should grow as PMIN gets more negative: {E}"
    assert E[-1e4] - E[-5.0] > 0.2, "PMIN should matter: expected >0.2 mm difference over 1 day"

def test_surface_head_never_below_finite_pmin(run):
    for p in (-5.0, -150.0):
        _, r = run("scf", 0.0, p)
        assert surface_psi_min(r) >= p * (1 + 1e-6), f"surface head went below PMIN={p}"

# ---------------------------------------------------------------- 6. PMIN does not enter etran.f: T at SCF=1 should not depend on it
def test_scf1_transpiration_independent_of_pmin(run):
    _, a = run("scf", 1.0, -1e35); _, b = run("scf", 1.0, -5.0)
    ta, Ea, tta, Ta = series(a); tb, Eb, ttb, Tb = series(b)
    assert cum(tta, Ta) > 0.5
    rel = abs(cum(tta, Ta) - cum(ttb, Tb)) / cum(tta, Ta)
    assert rel < 0.05, f"T depends on PMIN by {100*rel:.1f} % at SCF=1 (hidden coupling through the surface switch)"

# ---------------------------------------------------------------- 7. stress: PMIN disabled with SCF < 1 (slow; report-oriented)
@pytest.mark.stress
@pytest.mark.parametrize("scf", [0.0, 0.3])
def test_pmin_disabled_with_bare_soil_demand(run, scf):
    spec, r = run("scf", scf, -1e35, timeout=900)
    print(f"[{spec.name}] exit={r['exit']} completed={r['completed']} wall={r['wall_s']:.0f}s "
          f"steps={r.get('nsteps')} backsteps={r.get('backsteps')} nl_iter={r.get('nl_iter')} "
          f"psi_surf_min={surface_psi_min(r) if r.get('psi') is not None else None}")
    assert r["exit"] == 0 and r["completed"], "run stalled/failed with PMIN disabled"
    assert np.all(np.isfinite(r["psi"]))
