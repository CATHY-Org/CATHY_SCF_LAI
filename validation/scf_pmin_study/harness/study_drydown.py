"""Follow-up: long dry-down (12 d, one 25 mm pulse, 5 mm/d demand) so Feddes stress actually bites.
Runs: SCF {0.3, 1} x PMIN {-5, -150}.  Tests (i) Feddes-limited transpiration, (ii) whether a finite PMIN
(less negative than PCWLT=-150) lets the surface switch act on a stressed root zone at SCF=1."""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cathy_runner as c
ROOT = "/home/claude/study/runs_dry"; OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results_drydown")
os.makedirs(OUT, exist_ok=True); os.makedirs(ROOT, exist_ok=True)
DAYS = 12.0; F = c.semi_arid_forcing(((0.0, 0.25, 25.0),), 5.0, DAYS)
sp = os.path.join(OUT, "summary.json"); S = json.load(open(sp)) if os.path.exists(sp) else {}
for s, p in [(1.0, -5.0), (1.0, -150.0), (0.3, -5.0), (0.3, -150.0)]:
    nm = f"scf_scf{s:g}_pmin{p:g}"
    if nm in S: continue
    spec = c.RunSpec(name=nm, version="scf", scf=s, pmin=p, days=DAYS, forcing=F, psi0=-20.0, dtmax=1800.0, nprt=9)
    r = c.run(spec, ROOT, timeout=900)
    t = r["t"]; ev = np.maximum(-c.mm_day(r["act"]), 0); tr = c.mm_day(r["tr"])
    S[nm] = dict(scf=s, pmin=p, exit=r["exit"], completed=bool(r["completed"]), t_end=float(t[-1] / c.DAY), wall=round(r["wall_s"]),
                 E=float(np.trapezoid(ev, t) / c.DAY), T=float(np.trapezoid(tr, r["tr_t"]) / c.DAY),
                 T_day={d: float(np.interp(d * c.DAY, r["tr_t"], tr)) for d in (1, 3, 5, 8, 11)},
                 E_day={d: float(np.interp(d * c.DAY, t, ev)) for d in (1, 3, 5, 8, 11)},
                 psi_surf=[float(x) for x in r["psi"][1:, 0, 0]], psi_L5=[float(x) for x in r["psi"][1:, 4, 0]], psi_t=[float(x/c.DAY) for x in r["psi_t"][1:]],
                 resid_pct=float(abs(r["mb_cum"][-1] - np.trapezoid(r["tr"], r["tr_t"])) / r["V0"] * 100), nl=r["nl_iter"])
    np.savez_compressed(os.path.join(OUT, nm + ".npz"), t=t, act=r["act"], tr_t=r["tr_t"], tr=r["tr"], psi_t=r["psi_t"], psi=r["psi"][:, :9, 0])
    json.dump(S, open(sp, "w"), indent=1); print(nm, S[nm]["completed"], round(S[nm]["E"], 2), round(S[nm]["T"], 2), S[nm]["wall"], flush=True)
