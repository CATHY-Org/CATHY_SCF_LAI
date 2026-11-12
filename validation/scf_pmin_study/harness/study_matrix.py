"""Study matrix: SCF (0, 0.1, 0.3, 1) x PMIN (-5, -150, -1e4, -1e35) on the SCF build,
plus classic-build endpoint runs (SCF 0/1 x PMIN -5/-1e35) for the regression comparison.
Semi-arid forcing: 25 mm pulse (day 0-0.25), 5 mm/d dry spell, 8 mm pulse (day 4-4.1), 6 days."""
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cathy_runner as c

ROOT = os.environ.get("STUDY_ROOT", "/home/claude/study/runs_matrix")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
os.makedirs(OUT, exist_ok=True); os.makedirs(ROOT, exist_ok=True)
DAYS = 6.0
PULSES = ((0.0, 0.25, 25.0), (4.0, 4.1, 8.0))
FORC = c.semi_arid_forcing(PULSES, 5.0, DAYS)
PET = 5.0

def tag(v, s, p): return f"{v}_scf{s:g}_pmin{p:g}".replace("+", "")

def metrics(r, spec):
    t = r["t"]; act = c.mm_day(r["act"]); pot = c.mm_day(r["pot"])
    tr_t = r["tr_t"]; tr = c.mm_day(r["tr"])
    ev = np.maximum(-act, 0.0)                       # mm/d evaporation (ATMACT sink)
    cumE = float(np.trapezoid(ev, t) / c.DAY)
    cumT = float(np.trapezoid(tr, tr_t) / c.DAY)
    dry = (t > 0.30 * c.DAY) & (t < 4.0 * c.DAY)     # first dry spell
    first_day = (t > 0.30 * c.DAY) & (t < 1.0 * c.DAY)
    psi = r["psi"]
    m = dict(name=spec.name, version=spec.version, scf=spec.scf, pmin=spec.pmin,
             cumE_mm=cumE, cumT_mm=cumT, cumET_mm=cumE + cumT,
             E_frac=cumE / max(cumE + cumT, 1e-12),
             E_first_dryday_over_demand=float(np.trapezoid(ev[first_day], t[first_day]) / c.DAY / max((1 - spec.scf) * PET * (t[first_day][-1] - t[first_day][0]) / c.DAY, 1e-12)) if (1-spec.scf) > 0 and first_day.sum() > 2 else None,
             E_at_day3_mm_d=float(np.interp(3.0 * c.DAY, t, ev)) if t[-1] >= 3*c.DAY else None,
             T_at_day3_mm_d=float(np.interp(3.0 * c.DAY, tr_t, tr)) if tr_t[-1] >= 3*c.DAY else None,
             ET_max_over_PET=float(max(np.max((ev[dry] if dry.any() else 0)) + 0, 0) / PET),
             psi_surf_min=float(psi[:, 0, :].min()), psi_surf_end=float(psi[-1, 0, 0]),
             psi_L2_end=float(psi[-1, 1, 0]), psi_L9_end=float(psi[-1, 8, 0]),
             sw_surf_end=float(r["sw"][-1, 0, 0]),
             mb_residual_pct_V0=float(abs(r["mb_cum"][-1] - np.trapezoid(r["tr"], r["tr_t"])) / r["V0"] * 100),
             nsteps=r["nsteps"], backsteps=r["backsteps"], nl_iter=r["nl_iter"],
             wall_s=round(r["wall_s"], 1), exit=r["exit"])
    # time (days) when evaporation first drops below 95% of its own demand after day 0.3 (end of stage 1)
    d = (1 - spec.scf) * PET
    if d > 0:
        idx = np.where((t > 0.30 * c.DAY) & (ev < 0.95 * d))[0]
        m["stage1_end_day"] = float(t[idx[0]] / c.DAY) if len(idx) else None
    else: m["stage1_end_day"] = None
    return m

def main():
    runs = []
    for s in (0.0, 1.0):
        for p in (-5.0, -1e35): runs.append(("classic", s, p))
    for p in (-150.0, -1e4, -5.0, -1e35):
        for s in ((0.0, 0.3, 1.0, 0.1) if p == -150.0 else (0.0, 0.3, 1.0)): runs.append(("scf", s, p))
    # de-dup / order: finish endpoints and the three headline PMIN first
    seen, order = set(), []
    for r in runs:
        if r not in seen: seen.add(r); order.append(r)
    summ_path = os.path.join(OUT, "summary.json")
    summ = json.load(open(summ_path)) if os.path.exists(summ_path) else {}
    for (v, s, p) in order:
        nm = tag(v, s, p)
        if nm in summ and summ[nm].get("exit") in (0, -9): continue
        spec = c.RunSpec(name=nm, version=v, scf=s, pmin=p, days=DAYS, forcing=FORC,
                         psi0=-20.0, dtmax=1800.0, nprt=9)
        print(time.strftime("%H:%M:%S"), "RUN", nm, flush=True)
        try:
            r = c.run(spec, ROOT, timeout=int(os.environ.get('RUN_TIMEOUT','660')))
            m = metrics(r, spec) if "t" in r and len(r["t"]) > 3 else dict(name=nm, exit=r["exit"], wall_s=r["wall_s"])
            m["completed"] = bool(r.get("completed")); m["t_end_days"] = float(r["t"][-1] / c.DAY) if "t" in r else None
            if "t" in r and len(r["t"]) > 3 and r.get("psi") is not None:
                np.savez_compressed(os.path.join(OUT, nm + ".npz"), t=r["t"], act=r["act"], pot=r["pot"],
                                    tr_t=r["tr_t"], tr=r["tr"], psi_t=r["psi_t"], psi=r["psi"][:, :9, 0],
                                    sw=r["sw"][:, :9, 0], mb_t=r["mb_t"], mb_cum=r["mb_cum"])
        except Exception as e:
            m = dict(name=nm, exit=-1, error=repr(e))
        summ[nm] = m
        json.dump(summ, open(summ_path, "w"), indent=1)
        print("   ->", {k: m.get(k) for k in ("exit", "wall_s", "cumE_mm", "cumT_mm", "psi_surf_min", "backsteps")}, flush=True)
        # keep disk small: drop heavy stdout/outputs of finished run
        for f in ("stdout.txt",):
            try: os.remove(os.path.join(ROOT, nm, f))
            except OSError: pass

if __name__ == "__main__":
    main()
