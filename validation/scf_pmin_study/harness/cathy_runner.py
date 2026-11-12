"""Minimal driver to run CATHY (classic or SCF build) on the 20x20x15 test column-grid.

Generates input/soil, input/atmbc, input/ic, input/parm for a RunSpec, runs the binary,
and parses the diagnostics used by the smoke tests.  No pyCATHY dependency.
"""
import os, re, shutil, subprocess, time
from dataclasses import dataclass, field
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "bin")
TEMPLATE = os.path.join(HERE, "template")      # input/, prepro/, cathy.fnames
AREA = 100.0                                   # m2 of the catchment (10 m x 10 m)
DAY = 86400.0
NNOD, NLAY = 441, 16                           # surface nodes, node layers (NSTR+1)

def mm_day(q_m3s):  # domain-total flux (m3/s) -> mm/day over the 100 m2 surface
    return np.asarray(q_m3s) / AREA * DAY * 1000.0

@dataclass
class RunSpec:
    name: str
    version: str = "scf"          # 'classic' | 'scf'
    scf: float = 0.0
    pmin: float = -150.0
    days: float = 10.0
    forcing: list = field(default_factory=lambda: semi_arid_forcing())  # [(t_s, m/s)]
    psi0: float = -20.0           # uniform initial pressure head (m)
    isimgr: int = 2               # 2 coupled (SURF -> ADRSTN+SWITCH), 1 subsurface only (SWITCH_OLD)
    dtmax: float = 600.0
    zroot: float = 0.5
    pcwlt: float = -150.0
    pcref: float = -4.0
    nprt: int = 9                 # CATHY.H MAXPRT=9: more print times => silent early stop (exit 0!)
    scf_list: list = None         # SCF build only: one SCF per vegetation class (overrides scf)
    veg_map: list = None          # 20x20 list-of-lists of class ids (1..NVEG) written to input/root_map
    wt_height: float = None       # if set: hydrostatic IC (INDP=3), water table this many m ABOVE the column base (3 m column)

def semi_arid_forcing(pulses=((0.0, 0.25, 25.0),), pet_mm_day=5.0, days=10.0):
    """Step forcing (IETO=1, value held until next breakpoint), m/s, negative = evaporative demand.
    pulses: (start_day, end_day, total_mm) rain events; PET = constant dry-spell demand."""
    pet = -pet_mm_day / 1000.0 / DAY
    ev = sorted(pulses)
    out, t = [], 0.0
    for (d0, d1, mm) in ev:
        if d0 * DAY > t:
            out.append((t, pet)); t = d0 * DAY
        out.append((t, mm / 1000.0 / ((d1 - d0) * DAY))); t = d1 * DAY
    out.append((t, pet)); out.append((days * DAY + 10 * DAY, pet))
    return out

def _read(p): return open(p).read().split("\n")

def make_run(spec: RunSpec, root: str) -> str:
    assert 1 <= spec.nprt <= 9, 'CATHY.H MAXPRT=9'
    d = os.path.join(root, spec.name)
    if os.path.exists(d): shutil.rmtree(d)
    shutil.copytree(TEMPLATE, d)
    for sub in ("output", "vtk"): os.makedirs(os.path.join(d, sub), exist_ok=True)
    shutil.copy(os.path.join(BIN, f"cathy_{spec.version}"), os.path.join(d, "cathy_ft"))
    # --- soil
    soil = _read(os.path.join(TEMPLATE, "input", "soil"))
    veg = f"0.000e+00 {spec.pcref:.6e} {spec.pcwlt:.6e} {spec.zroot:.6e} 1.000e+00 1.000e+00"
    head = [f"{spec.pmin:.6e}\tPMIN"]
    if spec.version == "scf":
        sl = spec.scf_list if spec.scf_list else [spec.scf]
        new = head + ["0\tIPEAT", "0.4 0.225\tCBETA0 CANG"] + [veg] * len(sl) + [" ".join(f"{x:.6e}" for x in sl)] + soil[5:]
    else:  # classic: IPEAT and SCF on one line, no trailing SCF line
        new = head + [f"0 {spec.scf:.6e}\tIPEAT SCF", "0.4 0.225\tCBETA0 CANG", veg] + soil[5:]
    open(os.path.join(d, "input", "soil"), "w").write("\n".join(new))
    if spec.veg_map is not None:
        rm = _read(os.path.join(TEMPLATE, "input", "root_map"))
        body = [" " + " ".join(f"{float(v):.1f}" for v in row) for row in spec.veg_map]
        assert len(body) == 20 and all(len(r) == 20 for r in spec.veg_map)
        open(os.path.join(d, "input", "root_map"), "w").write("\n".join(rm[:6] + body) + "\n")
    # --- atmbc (HSPATM=1 homogeneous, IETO=1 step)
    lines = ["1 1\tHSPATM IETO"]
    for t, v in spec.forcing: lines += [f"{t:.6e}", f"{v:.9e}"]
    open(os.path.join(d, "input", "atmbc"), "w").write("\n".join(lines) + "\n")
    # --- ic: uniform psi0
    if spec.wt_height is None:
        open(os.path.join(d, "input", "ic"), "w").write(f"0 0 INDP IPOND\n{spec.psi0:.6e}\n")
    else:
        open(os.path.join(d, "input", "ic"), "w").write(f"3 0 INDP IPOND\n{spec.wt_height:.6e}\n")
    # --- parm
    p = _read(os.path.join(TEMPLATE, "input", "parm"))
    p[1] = f"    {spec.isimgr}         0.00      0                ISIMGR PONDH_MIN VELREC"
    tmax = spec.days * DAY
    p[11] = f"    0.01    .00001    {spec.dtmax:.1f}   {tmax:.1f}     DELTAT DTMIN  DTMAX  TMAX"
    times = np.linspace(tmax / spec.nprt, tmax, spec.nprt)
    prt = [f"    4        0         {spec.nprt}"] + [f"    {t:.1f}" for t in times]
    tail_i = [i for i, l in enumerate(p) if "NUMVP" in l][0] - 1
    open(os.path.join(d, "input", "parm"), "w").write("\n".join(p[:13] + prt + p[tail_i:]))
    return d

def run(spec: RunSpec, root: str, timeout=3600):
    """Run CATHY. On wall-clock timeout the run is killed and the PARTIAL outputs are parsed
    (exit = -9, 'completed' = False): a stalled solver is itself a result."""
    d = make_run(spec, root)
    t0 = time.time(); timed_out = False
    with open(os.path.join(d, "stdout.txt"), "w") as fo:
        try:
            pr = subprocess.run(["./cathy_ft"], cwd=d, stdout=fo, stderr=subprocess.STDOUT, timeout=timeout)
            rc = pr.returncode
        except subprocess.TimeoutExpired:
            timed_out, rc = True, -9
    res = parse(d, spec)
    res["exit"] = -9 if timed_out else (rc if "t" in res else (rc or -2))
    res["completed"] = (not timed_out) and rc == 0 and "t" in res and res["t"][-1] >= spec.days * DAY * 0.999
    res["wall_s"] = time.time() - t0
    return res

def _table(path, ncol=None, skip_hash=True):
    rows = []
    for l in open(path):
        if skip_hash and l.lstrip().startswith("#"): continue
        try: rows.append([float(x) for x in l.split()])
        except ValueError: continue
    rows = [r for r in rows if r]
    n = ncol or min(len(r) for r in rows)
    return np.array([r[:n] for r in rows])

def parse(d, spec):
    o = {"dir": d}
    try:
        a = _table(os.path.join(d, "output", "hgatmsf"), 10)
        o["t"] = a[:, 2]; o["pot"] = a[:, 3]; o["act"] = a[:, 4]      # m3/s, evaporation sink = negative
    except Exception as e: o["err_atm"] = str(e)
    try:
        tr = _table(os.path.join(d, "fort.666"), 2)
        o["tr_t"] = tr[:, 0]; o["tr"] = tr[:, 1]                       # m3/s, positive = uptake
    except Exception as e: o["err_tr"] = str(e)
    try:
        m = _table(os.path.join(d, "output", "mbeconv"))
        o["mb_t"] = m[:, 2]; o["mb_cum"] = m[:, -2]; o["mb_cum_abs"] = m[:, -1]   # CUM. MBE (m3); NOTE: excludes root uptake
        o["V0"] = float(open(os.path.join(d, "output", "mbeconv")).readline().split("=")[1])
        o["mb_cols"] = m.shape[1]
    except Exception as e: o["err_mb"] = str(e)
    try:
        dc = _table(os.path.join(d, "output", "dtcoupling"))
        o["nsteps"] = int(dc[-1, 0]); o["backsteps"] = int(dc[-1, 3]); o["nl_iter"] = int(dc[:, 5].sum()) if dc.shape[1] > 5 else -1
    except Exception as e: o["err_dt"] = str(e)
    o["psi"], o["psi_t"] = read_nodal(os.path.join(d, "output", "psi"))
    o["sw"], _ = read_nodal(os.path.join(d, "output", "sw"))
    return o

def read_nodal(path):
    """-> array (ntimes, NLAY, NNOD) and times; node k = layer*NNOD + surface node."""
    if not os.path.exists(path): return None, None
    blocks, times, cur, hdr = [], [], [], None
    for l in open(path):
        if "TIME" in l and "NSTEP" in l:
            if cur: blocks.append(np.array(cur))
            cur = []; times.append(float(l.split()[1])); continue
        try: cur += [float(x) for x in l.split()]
        except ValueError: continue
    if cur: blocks.append(np.array(cur))
    n = NNOD * NLAY
    keep = [(t, b) for t, b in zip(times, blocks) if len(b) == n]
    if not keep: return None, None
    return np.array([b.reshape(NLAY, NNOD) for _, b in keep]), np.array([t for t, _ in keep])
