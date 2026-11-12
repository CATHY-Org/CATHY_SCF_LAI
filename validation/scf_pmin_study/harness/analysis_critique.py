import json, os, numpy as np
H = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(H, "results_critique")
S1 = json.load(open(os.path.join(R, "summary.json"))); S2 = json.load(open(os.path.join(R, "summary2.json")))
M = json.load(open(os.path.join(H, "results", "summary.json"))); DD = json.load(open(os.path.join(H, "results_drydown", "summary.json")))
EX = json.load(open(os.path.join(H, "results_extra", "summary.json")))
V0 = 51.14
def resid_from_npz(path):
    z = np.load(path); return abs(z["mb_rel"][-1] - np.trapezoid(z["tr"], z["tr_t"])) / V0 * 100
def scf_matrix(nm):
    m = M[nm]; return m["cumE_mm"], m["cumT_mm"], m["nl_iter"], m["backsteps"], resid_from_npz(os.path.join(H, "results", nm + ".npz"))
def scf_dd(nm):
    m = DD[nm]; return m["E"], m["T"], m["nl"], 0, m["resid_pct"]
rows = ["| scenario (SCF=0.3) | classic E / T (mm) | classic NL iters | classic back-steps | classic budget resid % V0 | SCF-build E / T (mm) | SCF NL iters | SCF back-steps | SCF resid % V0 |", "|---|---|---|---|---|---|---|---|---|"]
def row(label, c, s):
    rows.append(f"| {label} | {c[0]:.2f} / {c[1]:.2f} | {c[2]} | {c[3]} | {c[4]:.3f} | {s[0]:.2f} / {s[1]:.2f} | {s[2]} | {s[3]} | {s[4]:.3f} |")
cl = lambda k: (S1[k]["E"], S1[k]["T"], S1[k]["nl_iter"], S1[k]["backsteps"], S1[k]["resid_pct"])
row("ZROOT 0.5, PMIN −5, 6 d", cl("D_classic_scf0.3_pmin-5"), scf_matrix("scf_scf0.3_pmin-5"))
row("ZROOT 0.5, PMIN −150, 6 d", cl("D_classic_scf0.3_pmin-150"), scf_matrix("scf_scf0.3_pmin-150"))
row("ZROOT 0.5, PMIN off, 6 d", cl("D_classic_scf0.3_pmin-1e35"), scf_matrix("scf_scf0.3_pmin-1e35"))
row("ZROOT 0.5, PMIN −5, 12 d", cl("D_classic_scf0.3_pmin-5_12d"), scf_dd("scf_scf0.3_pmin-5"))
row("ZROOT 0.5, PMIN −150, 12 d", cl("D_classic_scf0.3_pmin-150_12d"), scf_dd("scf_scf0.3_pmin-150"))
for z in (0.02, 0.006):
    c = S2[f"D2_classic_scf0.3_zroot{z:g}_pmin-5"]; s = S2[f"D2_scf_scf0.3_zroot{z:g}_pmin-5"]
    row(f"**ZROOT {z:g}**, PMIN −5, 6 d", (c["E"], c["T"], c["nl_iter"], c["backsteps"], c["resid_pct"]), (s["E"], s["T"], s["nl_iter"], s["backsteps"], s["resid_pct"]))
open(os.path.join(R, "classic_vs_scf_table.md"), "w").write("\n".join(rows) + "\n"); print("\n".join(rows))
zr = ["| ZROOT (m) | E (mm) | T (mm) | NL iters | budget resid % V0 |", "|---|---|---|---|---|"]
for z in (0.006, 0.02, 0.1):
    m = S1[f"E_scf1_zroot{z:g}"]; zr.append(f"| {z:g} | {m['E']:.2f} | {m['T']:.2f} | {m['nl_iter']} | {m['resid_pct']:.3f} |")
m = M["scf_scf1_pmin-150"]; zr.append(f"| 0.5 | {m['cumE_mm']:.2f} | {m['cumT_mm']:.2f} | {m['nl_iter']} | – |")
m = M["scf_scf0_pmin-150"]; zr.append(f"| SCF=0 (bare soil), ZROOT 0.5 | {m['cumE_mm']:.2f} | {m['cumT_mm']:.2f} | {m['nl_iter']} | – |")
open(os.path.join(R, "zroot_table.md"), "w").write("\n".join(zr) + "\n"); print("\n".join(zr))
