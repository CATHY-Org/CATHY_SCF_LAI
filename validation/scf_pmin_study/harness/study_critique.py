"""D. Does the CLASSIC build misbehave with a float SCF (colleague's claim) vs the SCF build?  (SCF=0.3)
   E. Can ZROOT (per root-map class) replace a spatial SCF?  (SCF=1, ZROOT varied; PMIN -150)"""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cathy_runner as c
ROOT = "/home/claude/study/runs_crit"; OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results_critique")
os.makedirs(ROOT, exist_ok=True); os.makedirs(OUT, exist_ok=True)
DAY = c.DAY; TO = 900
sp = os.path.join(OUT, "summary.json"); S = json.load(open(sp)) if os.path.exists(sp) else {}
F6 = c.semi_arid_forcing(((0.0, 0.25, 25.0), (4.0, 4.1, 8.0)), 5.0, 6.0)
F12 = c.semi_arid_forcing(((0.0, 0.25, 25.0),), 5.0, 12.0)

def go(nm, days, forcing, **kw):
    if nm in S and S[nm].get("exit") in (0, -9): return
    spec = c.RunSpec(name=nm, days=days, forcing=forcing, psi0=-20.0, dtmax=1800.0, nprt=9, **kw)
    r = c.run(spec, ROOT, timeout=TO)
    t = r["t"]; a = c.mm_day(r["act"]); ev = np.maximum(-a, 0); tt = r["tr_t"]; tr = c.mm_day(r["tr"])
    w = t > 0.3 * DAY; tot = ev + np.interp(t, tt, tr)
    m = dict(name=nm, version=spec.version, scf=spec.scf, pmin=spec.pmin, zroot=spec.zroot, days=days, exit=r["exit"], completed=bool(r["completed"]),
             wall_s=round(r["wall_s"]), nsteps=r["nsteps"], backsteps=r["backsteps"], nl_iter=r["nl_iter"],
             E=float(np.trapezoid(ev, t) / DAY), T=float(np.trapezoid(tr, tt) / DAY),
             resid_pct=float(abs(r["mb_cum"][-1] - np.trapezoid(r["tr"], tt)) / r["V0"] * 100),
             max_ET_over_PET=float(tot[w].max() / 5.0), phantom_mm=float(np.trapezoid(np.where(w, np.maximum(a, 0), 0), t) / DAY),
             psi_surf_min=float(r["psi"][1:, 0].min()), psi_L2_min=float(r["psi"][1:, 1].min()))
    np.savez_compressed(os.path.join(OUT, nm + ".npz"), t=t, act=r["act"], tr_t=tt, tr=r["tr"], psi_t=r["psi_t"], psi=r["psi"][:, :9, 0], mb_cum=r["mb_cum"])
    S[nm] = m; json.dump(S, open(sp, "w"), indent=1)
    print(nm, m["completed"], round(m["E"], 2), round(m["T"], 2), "nl", m["nl_iter"], "bk", m["backsteps"], "resid", round(m["resid_pct"], 4), flush=True)

# E. ZROOT instead of SCF (SCF=1: all demand is transpiration; only the depth of uptake changes)
for z in (0.006, 0.02, 0.1):
    go(f"E_scf1_zroot{z:g}", 6.0, F6, version="scf", scf=1.0, pmin=-150.0, zroot=z)
# D. classic build, float SCF = 0.3
go("D_classic_scf0.3_pmin-5", 6.0, F6, version="classic", scf=0.3, pmin=-5.0)
go("D_classic_scf0.3_pmin-150", 6.0, F6, version="classic", scf=0.3, pmin=-150.0)
go("D_classic_scf0.3_pmin-150_12d", 12.0, F12, version="classic", scf=0.3, pmin=-150.0)
go("D_classic_scf0.3_pmin-5_12d", 12.0, F12, version="classic", scf=0.3, pmin=-5.0)
go("D_classic_scf0.3_pmin-1e35", 6.0, F6, version="classic", scf=0.3, pmin=-1e35)
