"""D2. Worst case for the 'uncoupled limiters' concern: SCF=0.3, thin root zone (top layers only), PMIN above PCWLT,
so the pinned surface node still has a non-zero Feddes factor.  Classic vs SCF build, 6-day forcing."""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cathy_runner as c
ROOT = "/home/claude/study/runs_crit2"; OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results_critique")
os.makedirs(ROOT, exist_ok=True)
sp = os.path.join(OUT, "summary2.json"); S = json.load(open(sp)) if os.path.exists(sp) else {}
DAY = c.DAY; F6 = c.semi_arid_forcing(((0.0, 0.25, 25.0), (4.0, 4.1, 8.0)), 5.0, 6.0)
for z in (0.006, 0.02):
    for v in ("classic", "scf"):
        nm = f"D2_{v}_scf0.3_zroot{z:g}_pmin-5"
        if nm in S: continue
        spec = c.RunSpec(name=nm, version=v, scf=0.3, pmin=-5.0, zroot=z, days=6.0, forcing=F6, psi0=-20.0, dtmax=1800.0, nprt=9)
        r = c.run(spec, ROOT, timeout=1200)
        t = r["t"]; a = c.mm_day(r["act"]); ev = np.maximum(-a, 0); tt = r["tr_t"]; tr = c.mm_day(r["tr"]); w = t > 0.3 * DAY
        tot = ev + np.interp(t, tt, tr)
        S[nm] = dict(version=v, zroot=z, completed=bool(r["completed"]), t_end=float(t[-1] / DAY), wall_s=round(r["wall_s"]), nl_iter=r["nl_iter"], backsteps=r["backsteps"],
                     E=float(np.trapezoid(ev, t) / DAY), T=float(np.trapezoid(tr, tt) / DAY),
                     resid_pct=float(abs(r["mb_cum"][-1] - np.trapezoid(r["tr"], tt)) / r["V0"] * 100),
                     max_ET_over_PET=float(tot[w].max() / 5.0), psi_surf_min=float(r["psi"][1:, 0].min()))
        json.dump(S, open(sp, "w"), indent=1); print(nm, S[nm]["completed"], round(S[nm]["E"], 2), round(S[nm]["T"], 2), "nl", S[nm]["nl_iter"], "bk", S[nm]["backsteps"], "resid", round(S[nm]["resid_pct"], 4), flush=True)
