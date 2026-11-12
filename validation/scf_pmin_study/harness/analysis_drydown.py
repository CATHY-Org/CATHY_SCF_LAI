import os, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(H, "results_drydown"); DAY, A = 86400., 100.
def L(k): z = np.load(os.path.join(R, k + ".npz")); return z
fig, ax = plt.subplots(1, 3, figsize=(14, 3.8))
col = {-5.0: "tab:red", -150.0: "tab:blue"}
for s, ls in ((1.0, "-"), (0.3, "--")):
    for p in (-5.0, -150.0):
        z = L(f"scf_scf{s:g}_pmin{p:g}"); a = z["act"] / A * DAY * 1000; t = z["t"]
        ax[0].plot(z["tr_t"] / DAY, z["tr"] / A * DAY * 1000, ls, color=col[p], label=f"SCF={s:g}, PMIN {p:g} m")
        m = t > 0.3 * DAY
        c = np.concatenate([[0], np.cumsum(np.where(m, np.maximum(a, 0), 0)[1:] * np.diff(t)) / DAY])
        ax[1].plot(t / DAY, c, ls, color=col[p], label=f"SCF={s:g}, PMIN {p:g} m")
        if s == 1.0: ax[2].plot(z["psi_t"][1:] / DAY, -z["psi"][1:, 0], "o-", color=col[p], label=f"surface node, PMIN {p:g}")
ax[0].set_title("Transpiration"); ax[0].set_ylabel("mm/d"); ax[0].set_xlim(0.3, 12); ax[0].legend(fontsize=7)
ax[1].set_title("Cumulative infiltration through the surface\nafter the pulse (phantom source)"); ax[1].set_ylabel("mm"); ax[1].legend(fontsize=7)
z = L("scf_scf1_pmin-5"); ax[2].plot(z["psi_t"][1:] / DAY, -z["psi"][1:, 4], "k:^", label="layer 5 (root zone)")
ax[2].set_yscale("log"); ax[2].set_title("SCF=1: head magnitude |ψ| (m)"); ax[2].legend(fontsize=7)
for a in ax: a.set_xlabel("day")
plt.tight_layout(); plt.savefig(os.path.join(R, "fig5_drydown_phantom_infiltration.png"), dpi=130)
