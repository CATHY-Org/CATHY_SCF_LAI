"""Figures + summary tables from results/summary.json and results/*.npz."""
import json, os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(HERE, "results")
DAY, AREA = 86400.0, 100.0
mm = lambda q: np.asarray(q) / AREA * DAY * 1000.0
S = json.load(open(os.path.join(R, "summary.json")))
def nm(v, s, p): return f"{v}_scf{s:g}_pmin{p:g}".replace("+", "")
def load(v, s, p):
    f = os.path.join(R, nm(v, s, p) + ".npz")
    if not os.path.exists(f): return None
    z = np.load(f); return {k: z[k] for k in z.files}
def lab(p): return "off (-1e35)" if p <= -1e30 else f"{p:g} m"
COL = {-5.0: "tab:red", -150.0: "tab:orange", -1e4: "tab:green", -1e35: "tab:blue"}
SCFC = {0.0: "tab:brown", 0.1: "tab:olive", 0.3: "tab:green", 1.0: "tab:blue"}

def ev(z): return np.maximum(-mm(z["act"]), 0)

# Fig 1: partition vs SCF at PMIN=-150
fig, ax = plt.subplots(1, 2, figsize=(11, 3.8), sharex=True)
for s in (0.0, 0.1, 0.3, 1.0):
    z = load("scf", s, -150.0)
    if z is None: continue
    ax[0].plot(z["t"] / DAY, ev(z), color=SCFC[s], label=f"SCF={s:g}")
    ax[1].plot(z["tr_t"] / DAY, mm(z["tr"]), color=SCFC[s], label=f"SCF={s:g}")
ax[0].set_title("Surface evaporation  (-ATMACT), PMIN=-150 m"); ax[1].set_title("Transpiration (root uptake)")
for a in ax: a.set_xlabel("day"); a.set_ylabel("mm/d"); a.set_xlim(0.3, 6); a.axhline(5, ls=":", c="k", lw=.7); a.legend(fontsize=8)
ax[0].set_ylim(0, 6); ax[1].set_ylim(0, 6); plt.tight_layout(); plt.savefig(os.path.join(R, "fig1_partition_vs_scf.png"), dpi=130); plt.close()

# Fig 2: PMIN effect at SCF=0 and SCF=0.3 (evaporation)
fig, ax = plt.subplots(1, 2, figsize=(11, 3.8), sharex=True, sharey=True)
for k, s in enumerate((0.0, 0.3)):
    for p in (-5.0, -150.0, -1e4, -1e35):
        z = load("scf", s, p)
        if z is None: continue
        part = "" if S.get(nm("scf", s, p), {}).get("completed", True) else "  [stalled]"
        ax[k].plot(z["t"] / DAY, ev(z), color=COL[p], label=f"PMIN {lab(p)}{part}")
    ax[k].axhline((1 - s) * 5, ls=":", c="k", lw=.7); ax[k].set_title(f"Evaporation, SCF={s:g}  (dotted = its demand)")
    ax[k].set_xlabel("day"); ax[k].set_xlim(0.3, 6); ax[k].set_ylim(0, 5.5); ax[k].legend(fontsize=8)
ax[0].set_ylabel("mm/d"); plt.tight_layout(); plt.savefig(os.path.join(R, "fig2_pmin_effect_evaporation.png"), dpi=130); plt.close()

# Fig 3: cumulative ET components, grid SCF x PMIN
fig, ax = plt.subplots(1, 4, figsize=(13, 3.6), sharey=True)
for k, p in enumerate((-5.0, -150.0, -1e4, -1e35)):
    E, T, X = [], [], []
    for s in (0.0, 0.1, 0.3, 1.0):
        m = S.get(nm("scf", s, p));
        if not m or m.get("cumE_mm") is None: continue
        X.append(f"{s:g}" + ("*" if not m.get("completed", True) else "")); E.append(m["cumE_mm"]); T.append(m["cumT_mm"])
    ax[k].bar(X, E, label="evaporation", color="tab:orange"); ax[k].bar(X, T, bottom=E, label="transpiration", color="tab:green")
    ax[k].set_title(f"PMIN {lab(p)}"); ax[k].set_xlabel("SCF  (* = run stalled, partial)")
ax[0].set_ylabel("cumulative mm (6 d)"); ax[0].legend(fontsize=8); plt.tight_layout()
plt.savefig(os.path.join(R, "fig3_cumulative_E_T.png"), dpi=130); plt.close()

# Fig 4: surface-layer head at print times, SCF=0
fig, ax = plt.subplots(figsize=(6, 3.8))
for p in (-5.0, -150.0, -1e4, -1e35):
    z = load("scf", 0.0, p)
    if z is None: continue
    ax.plot(z["psi_t"] / DAY, np.maximum(z["psi"][:, 0], -1e6), "o-", color=COL[p], label=f"PMIN {lab(p)}")
ax.set_yscale("symlog", linthresh=1); ax.set_xlabel("day"); ax.set_ylabel("surface node head (m, symlog)")
ax.set_title("SCF=0: surface head vs PMIN"); ax.legend(fontsize=8); plt.tight_layout()
plt.savefig(os.path.join(R, "fig4_surface_head.png"), dpi=130); plt.close()

# Table
rows = []
hdr = "| build | SCF | PMIN | done | E mm | T mm | E/(E+T) | stage-1 end (d) | E@3d | T@3d | surf ψ min (m) | back-steps | NL iters | wall s |"
rows += [hdr, "|" + "---|" * (hdr.count("|") - 1)]
def g(m, k, f="{:.2f}"): 
    v = m.get(k); return "–" if v is None else f.format(v)
for n, m in sorted(S.items(), key=lambda kv: (kv[1].get("version", ""), kv[1].get("scf", 0), -kv[1].get("pmin", 0))):
    z = load(m.get("version", "scf"), m.get("scf", 0), m.get("pmin", 0))
    smin = "–" if z is None else f"{z['psi'][1:, 0].min():.4g}"
    rows.append(f"| {m.get('version')} | {m.get('scf'):g} | {lab(m.get('pmin'))} | {'yes' if m.get('completed', m.get('exit')==0) else 'NO ('+g(m,'t_end_days','{:.2f}')+' d)'} | "
                f"{g(m,'cumE_mm')} | {g(m,'cumT_mm')} | {g(m,'E_frac')} | {g(m,'stage1_end_day')} | {g(m,'E_at_day3_mm_d')} | {g(m,'T_at_day3_mm_d')} | {smin} | "
                f"{m.get('backsteps','–')} | {m.get('nl_iter','–')} | {g(m,'wall_s','{:.0f}')} |")
open(os.path.join(R, "summary_table.md"), "w").write("\n".join(rows) + "\n")
print("\n".join(rows))
