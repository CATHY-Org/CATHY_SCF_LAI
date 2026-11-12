"""Two targeted follow-up experiments on the SCF build (all PMIN = -150 m = PCWLT, ISIMGR=2, 12 days).
A. Release rule: PET staircase (8 -> 2 -> 1 -> 0.6 -> 0.4 -> 0.3 mm/d) for SCF {0, 0.3, 0.6}.
   Invariant under a *consistent* rule: E(t) <= (1-SCF)*PET(t) and E+T <= PET.
B. Heterogeneous cover at nominal mean 0.5: clumped patches (0|1), 3-cell stripes (0|1), cell checkerboard (0|1),
   vs uniform SCF at the nominal and at the *realized* mean cover (class boundaries are truncated by INT() in datin.f)."""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cathy_runner as c

ROOT = "/home/claude/study/runs_extra"; OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results_extra")
os.makedirs(ROOT, exist_ok=True); os.makedirs(OUT, exist_ok=True)
DAYS = 12.0; DAY = c.DAY
TO = int(os.environ.get("RUN_TIMEOUT", "900"))
sp = os.path.join(OUT, "summary.json"); S = json.load(open(sp)) if os.path.exists(sp) else {}

STAIR = [(0.25, 8.0), (3.0, 2.0), (4.5, 1.0), (6.0, 0.6), (7.5, 0.4), (9.0, 0.3)]      # (start day, PET mm/d)
def stair_forcing(days=DAYS):
    f = [(0.0, 25.0 / 1000 / (0.25 * DAY))]
    for d, p in STAIR: f.append((d * DAY, -p / 1000 / DAY))
    f.append(((days + 10) * DAY, f[-1][1])); return f
def pet_at(t_days):
    p = np.zeros_like(t_days, dtype=float)
    for d, v in STAIR: p = np.where(t_days >= d, v, p)
    return np.where(t_days < 0.25, 0.0, p)

def phantom(a_mm, t):
    m = t > 0.3 * DAY
    return float(np.trapezoid(np.where(m, np.maximum(a_mm, 0), 0), t) / DAY)

def finish(nm, spec, r, extra=None):
    t = r["t"]; a = c.mm_day(r["act"]); ev = np.maximum(-a, 0); tr_t = r["tr_t"]; tr = c.mm_day(r["tr"])
    m = dict(name=nm, scf=spec.scf_list or spec.scf, exit=r["exit"], completed=bool(r["completed"]), wall_s=round(r["wall_s"]),
             nsteps=r["nsteps"], backsteps=r["backsteps"], nl_iter=r["nl_iter"],
             E=float(np.trapezoid(ev, t) / DAY), T=float(np.trapezoid(tr, tr_t) / DAY), phantom_mm=phantom(a, t),
             resid_pct=float(abs(r["mb_cum"][-1] - np.trapezoid(r["tr"], tr_t)) / r["V0"] * 100))
    w = (t > 0.3 * DAY) & (t < 1.0 * DAY)
    m["E_over_PET_d03_1"] = float(ev[w].mean() / (5.0 if "stair" not in nm else 8.0))
    if r["psi"] is not None:
        m["psi_surf_final_min"] = float(r["psi"][-1, 0].min()); m["psi_surf_final_med"] = float(np.median(r["psi"][-1, 0]))
        m["psi_surf_final_max"] = float(r["psi"][-1, 0].max())
    if extra: m.update(extra)
    np.savez_compressed(os.path.join(OUT, nm + ".npz"), t=t, act=r["act"], tr_t=tr_t, tr=r["tr"], psi_t=r["psi_t"],
                        psi=r["psi"][:, :9, :], mb_cum=r["mb_cum"])
    S[nm] = m; json.dump(S, open(sp, "w"), indent=1); print(nm, m["completed"], round(m["E"], 2), round(m["T"], 2), m["wall_s"], flush=True)
    return m

def go(nm, **kw):
    if nm in S and S[nm].get("exit") in (0, -9): return S[nm]
    spec = c.RunSpec(name=nm, version="scf", pmin=-150.0, days=DAYS, psi0=-20.0, dtmax=1800.0, nprt=9, **kw)
    return finish(nm, spec, c.run(spec, ROOT, timeout=TO))

# ---------------- A. release rule
for s in (0.0, 0.3, 0.6):
    go(f"A_stair_scf{s:g}", scf=s, forcing=stair_forcing())

# ---------------- B. heterogeneity
F = c.semi_arid_forcing(((0.0, 0.25, 25.0),), 5.0, DAYS)
cols = np.arange(20)
MAPS = {"clumped": [[1] * 10 + [2] * 10 for _ in range(20)],
        "stripes3": [[1 if (j // 3) % 2 == 0 else 2 for j in cols] for _ in range(20)],
        "checker": [[1 if (i + j) % 2 == 0 else 2 for j in cols] for i in range(20)]}
real = {}
for k, mp in MAPS.items():
    m = go(f"B_{k}", scf_list=[0.0, 1.0], veg_map=mp, forcing=F)
    real[k] = 1.0 - m["E_over_PET_d03_1"]
    S[f"B_{k}"]["realized_scf"] = real[k]
json.dump(S, open(sp, "w"), indent=1)
print("realized mean SCF:", {k: round(v, 3) for k, v in real.items()}, flush=True)
go("B_uniform_0.5", scf=0.5, forcing=F)
for k, v in real.items():
    if 0.02 < v < 0.98 and abs(v - 0.5) > 0.02:
        go(f"B_uniform_real_{k}", scf=round(v, 3), forcing=F)


# ---------------- C. release rule with a shallow water table (capillary supply 1 m below the surface)
STAIR_C = [(0.0, 8.0), (2.0, 4.0), (4.0, 2.5), (6.0, 1.5), (8.0, 1.0), (10.0, 0.6)]
def stairC_forcing(days=DAYS):
    f = [(d * DAY, -p / 1000 / DAY) for d, p in STAIR_C]; f.append(((days + 10) * DAY, f[-1][1])); return f
for s in (0.0, 0.3, 0.6):
    go(f"C_wt_scf{s:g}", scf=s, forcing=stairC_forcing(), wt_height=2.0)
