import os, json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(H, "results_extra"); DAY, A = 86400., 100.
S = json.load(open(os.path.join(R, "summary.json")))
STAIR = [(0.25, 8.0), (3.0, 2.0), (4.5, 1.0), (6.0, 0.6), (7.5, 0.4), (9.0, 0.3)]
def pet_at(t):
    p = np.zeros_like(t, dtype=float)
    for d, v in STAIR: p = np.where(t >= d, v, p)
    return np.where(t < 0.25, 0.0, p)
mm = lambda q: np.asarray(q) / A * DAY * 1000
def L(k): return np.load(os.path.join(R, k + ".npz"))

# ---- A and C. release rule
STAIR_C = [(0.0, 8.0), (2.0, 4.0), (4.0, 2.5), (6.0, 1.5), (8.0, 1.0), (10.0, 0.6)]
def pet_c(t):
    p = np.zeros_like(t, dtype=float)
    for d, v in STAIR_C: p = np.where(t >= d, v, p)
    return p
def release_block(prefix, petf, stages, t_from, tag):
    lines = ["| SCF | E excess over (1−SCF)·PET after day %g (mm) | max E/((1−SCF)·PET) | E+T excess over PET (mm) | max (E+T)/PET | E / allowed, per PET stage (mm/d) |" % t_from, "|---|---|---|---|---|---|"]
    fig, ax = plt.subplots(1, 3, figsize=(14, 3.8), sharex=True)
    for k, s_ in enumerate((0.0, 0.3, 0.6)):
        nm = f"{prefix}{s_:g}"
        if nm not in S: continue
        z = L(nm); t = z["t"] / DAY; E = np.maximum(-mm(z["act"]), 0); T = np.interp(z["t"], z["tr_t"], mm(z["tr"])); pet = petf(t)
        allow = (1 - s_) * pet; m = t > t_from
        ex = np.trapezoid(np.maximum(E - allow * 1.01, 0)[m], t[m]); ext = np.trapezoid(np.maximum(E + T - pet * 1.01, 0)[m], t[m])
        ratio = np.max(E[m] / np.maximum(allow[m], 1e-9)) if s_ < 1 else np.nan; rt = np.max((E + T)[m] / np.maximum(pet[m], 1e-9))
        stg = []
        for d0, d1 in stages:
            w = (t > d0 + 0.3) & (t < d1); stg.append(f"{E[w].mean():.2f}/{allow[w].mean():.2f}")
        lines.append(f"| {s_:g} | {ex:.3f} | {ratio:.2f} | {ext:.3f} | {rt:.2f} | " + ", ".join(stg) + " |")
        ax[k].step(t, pet, "k:", where="post", label="PET"); ax[k].step(t, allow, "g--", where="post", label="(1−SCF)·PET (E allowed)")
        ax[k].plot(t, E, "r", label="E actual"); ax[k].plot(t, E + T, "b", lw=.8, label="E + T")
        ax[k].set_yscale("log"); ax[k].set_ylim(0.05, 10); ax[k].set_title(f"{tag}, SCF={s_:g}"); ax[k].set_xlabel("day"); ax[k].set_xlim(0.3, 12)
    ax[0].set_ylabel("mm/d (log)"); ax[0].legend(fontsize=7); plt.tight_layout()
    plt.savefig(os.path.join(R, f"fig6_release_rule_{tag}.png"), dpi=130); plt.close()
    open(os.path.join(R, f"release_table_{tag}.md"), "w").write("\n".join(lines) + "\n"); print("\n".join(lines))
release_block("A_stair_scf", pet_at, ((3, 4.5), (4.5, 6), (6, 7.5), (7.5, 9), (9, 12)), 3.0, "A_dry_start")
release_block("C_wt_scf", pet_c, ((2, 4), (4, 6), (6, 8), (8, 10), (10, 12)), 2.0, "C_water_table")

# ---- B. heterogeneity
names = [n for n in ("B_uniform_0.5", "B_uniform_real_clumped", "B_clumped", "B_uniform_real_stripes3", "B_stripes3", "B_checker") if n in S]
fig, ax = plt.subplots(1, 2, figsize=(12, 3.8))
rows = ["| run | realized mean SCF | E (mm) | T (mm) | E+T (mm) | NL iters | back-steps | wall (s) | budget resid % V0 | surf ψ final min / median / max (m) |", "|---|---|---|---|---|---|---|---|---|---|"]
for n in names:
    m = S[n]; sc = m.get("realized_scf", (m["scf"] if isinstance(m["scf"], float) else None))
    rows.append(f"| {n} | {'–' if sc is None else format(sc, '.3f')} | {m['E']:.2f} | {m['T']:.2f} | {m['E']+m['T']:.2f} | {m['nl_iter']} | {m['backsteps']} | {m['wall_s']} | {m['resid_pct']:.3f} | {m.get('psi_surf_final_min',float('nan')):.0f} / {m.get('psi_surf_final_med',float('nan')):.0f} / {m.get('psi_surf_final_max',float('nan')):.0f} |")
    z = L(n); ax[0].plot(z["t"] / DAY, np.maximum(-mm(z["act"]), 0), label=n.replace("B_", "")); ax[1].plot(z["tr_t"] / DAY, mm(z["tr"]), label=n.replace("B_", ""))
ax[0].set_title("Evaporation (domain mean)"); ax[1].set_title("Transpiration (domain mean)")
for a in ax: a.set_xlabel("day"); a.set_ylabel("mm/d"); a.set_xlim(0.3, 12)
ax[0].legend(fontsize=6); plt.tight_layout(); plt.savefig(os.path.join(R, "fig7_heterogeneity.png"), dpi=130); plt.close()
open(os.path.join(R, "hetero_table.md"), "w").write("\n".join(rows) + "\n"); print("\n".join(rows))
