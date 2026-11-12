"""
scf_pmin_study -- reproduce the SCF / PMIN study through the pyCATHY wrapper (``CATHY`` class).

Everything that *sets up and runs* CATHY goes through ``pyCATHY.cathy_tools.CATHY``:
    CATHY(...)            project + sources (version="SCF_variable" or "1.0.0" = classic)
    update_prepo_inputs / update_veg_map / update_soil / update_parm / update_atmbc / update_ic / update_sfbc
    run_processor(recompile=False)
Only two things are done outside the class, both documented in README.md:
    * a one-off silencing of the per-node debug prints in the SCF build's etran.f (44 MB/hour otherwise);
    * each scenario runs in a child process so that a stalled solver can be killed (run_processor has no timeout).

Usage (from the folder containing this file):
    python scf_pmin_study.py build                      # build both master projects (clone + compile, ~1 min)
    python scf_pmin_study.py run quick                  # 1-day checks (replaces the old smoke tests), a few minutes
    python scf_pmin_study.py run matrix|drydown|release|hetero|critique     # the study blocks
    python scf_pmin_study.py run all                    # everything (~2 h on one core)
    python scf_pmin_study.py compare                    # compare with the reference numbers of the direct-runner study
Results: ./results/<suite>.json|csv and ./results/checks.md
"""
import argparse, glob, json, os, re, shutil, subprocess, sys, time, warnings
from dataclasses import dataclass, field, asdict
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.environ.get("SCF_STUDY_WORK", os.path.join(HERE, "work"))
RESULTS = os.path.join(HERE, "results")
DAY, AREA = 86400.0, 100.0           # s/day ; 10 m x 10 m domain
NMAX_PRT = 9                          # CATHY.H MAXPRT of the compiled master (pyCATHY would re-compile to change it)
VERSIONS = {"scf": "SCF_variable", "classic": "1.0.0"}     # pyCATHY version strings


def mm_day(q_m3s):
    return np.asarray(q_m3s) / AREA * DAY * 1000.0


# ----------------------------------------------------------------------------------------------------------------
# forcing builders (net atmospheric flux in m/s, negative = evaporative demand; step function, IETO=1)
# ----------------------------------------------------------------------------------------------------------------
def pulses_forcing(pulses=((0.0, 0.25, 25.0),), pet_mm_day=5.0, days=6.0):
    """pulses = (start_day, end_day, total_mm) rain events; constant PET otherwise. -> (times_s, net_flux)"""
    pet = -pet_mm_day / 1000.0 / DAY
    t, out = 0.0, []
    for d0, d1, mm in sorted(pulses):
        if d0 * DAY > t:
            out.append((t, pet)); t = d0 * DAY
        out.append((t, mm / 1000.0 / ((d1 - d0) * DAY))); t = d1 * DAY
    out.append((t, pet)); out.append(((days + 10) * DAY, pet))
    return out


def stair_forcing(stair, days=12.0, rain_first_mm=None):
    """PET staircase [(start_day, mm/d)...]; optional initial 25 mm pulse over 0.25 d."""
    f = [(0.0, rain_first_mm / 1000.0 / (0.25 * DAY))] if rain_first_mm else []
    f += [(d * DAY, -p / 1000.0 / DAY) for d, p in stair]
    f.append(((days + 10) * DAY, f[-1][1]))
    return f


def pet_of_stair(stair, t_days):
    p = np.zeros_like(t_days, dtype=float)
    for d, v in stair:
        p = np.where(t_days >= d, v, p)
    return p


STAIR_A = [(0.25, 8.0), (3.0, 2.0), (4.5, 1.0), (6.0, 0.6), (7.5, 0.4), (9.0, 0.3)]      # dry start
STAIR_C = [(0.0, 8.0), (2.0, 4.0), (4.0, 2.5), (6.0, 1.5), (8.0, 1.0), (10.0, 0.6)]      # shallow water table

COLS = np.arange(20)
VEG_MAPS = {   # vegetation rasters (class ids) used for the heterogeneity block; 20 x 20 cells
    "clumped": np.array([[1] * 10 + [2] * 10 for _ in range(20)]),
    "stripes3": np.array([[1 if (j // 3) % 2 == 0 else 2 for j in COLS] for _ in range(20)]),
    "checker": np.array([[1 if (i + j) % 2 == 0 else 2 for j in COLS] for i in range(20)]),
}


# ----------------------------------------------------------------------------------------------------------------
@dataclass
class Scenario:
    name: str
    version: str = "scf"                 # 'scf' -> SCF_variable ; 'classic' -> 1.0.0
    scf: float = 0.0                     # uniform SCF (classic: the one global scalar)
    scf_list: list = None                # SCF per vegetation class (SCF build + veg_map)
    veg_map: str = None                  # key of VEG_MAPS (None: one class)
    pmin: float = -150.0
    zroot: float = 0.5
    pcwlt: float = -150.0
    pcref: float = -4.0
    days: float = 6.0
    forcing: list = field(default_factory=list)      # [(t_s, net_flux m/s)]
    psi0: float = -20.0                  # uniform initial head (m) unless wt_height given
    wt_height: float = None              # water table height above the 3 m column base (INDP=3)
    isimgr: int = 2
    dtmax: float = 1800.0
    nprt: int = 9
    timeout: int = 900
    tags: dict = field(default_factory=dict)         # free metadata (suite, expectation...)


# ----------------------------------------------------------------------------------------------------------------
# projects
# ----------------------------------------------------------------------------------------------------------------
def master_dir(version):
    return os.path.join(WORK, f"master_{version}")


def build_master(version):
    """CATHY(version=...) -> clone sources + input deck, silence debug prints (scf only), compile with recompileSrc."""
    from pyCATHY import cathy_tools
    prj = f"master_{version}"
    os.makedirs(WORK, exist_ok=True)
    cwd = os.getcwd()
    shutil.rmtree(os.path.join(WORK, prj), ignore_errors=True)         # always start from a clean clone
    simu = cathy_tools.CATHY(dirName=WORK, prj_name=prj, version=VERSIONS[version], verbose=False)
    if version == "scf":                          # one-off: silence per-node debug WRITE(6,*) in etran.f
        p = os.path.join(WORK, prj, "src", "etran.f")
        out, skip = [], 0
        for l in open(p).read().split("\n"):
            if skip:
                skip -= 1; continue
            if re.match(r"\s*write\(6,\*\)\s*'(I=|ZROOT=)", l, re.I):
                if "BTRAN=" in l: skip = 1
                continue
            out.append(l)
        open(p, "w").write("\n".join(out))
    # COMPILATION. Two routes were tried:
    #  (1) simu.run_processor(recompile=True): calls update_cathyH() first, which rewrites src/CATHY.H with tight
    #      dimensions (MAXVEG=1, NFACEMAX=2205) -> unusable for the 2-class vegetation maps;
    #  (2) simu.recompileSrc() on the upstream header: compiles with "gfortran -O -mcmodel=large -fPIC" and the binary
    #      SEGFAULTS (connection_new_) with gfortran 13.3 on this grid, whereas the upstream Makefile (-O3) runs.
    # => the master is compiled with the upstream Makefile and installed as ./cathy, the name run_processor() expects.
    #    Set SCF_STUDY_COMPILER=pycathy to use route (2) (useful on machines where it works).
    os.chdir(os.path.join(WORK, prj, "src"))
    try:
        if os.environ.get("SCF_STUDY_COMPILER", "make") == "pycathy":
            simu.recompileSrc(verbose=True)
        else:
            if os.path.exists("Conf_Fac.f") and not os.path.exists("conf_fac.f"):
                os.symlink("Conf_Fac.f", "conf_fac.f")              # Makefile asks for the lower-case name
            for f in glob.glob("*.o") + ["cathy_ft"]:
                if os.path.exists(f): os.remove(f)
            subprocess.run(["make", "-j2"], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            shutil.copy("cathy_ft", os.path.join("..", "cathy"))
    finally:
        os.chdir(cwd)
    assert os.path.exists(os.path.join(WORK, prj, "cathy")), "compilation failed"
    return os.path.join(WORK, prj)


def configure(sc: Scenario, run_dir):
    """Write the input deck of one scenario with the CATHY class (no recompile)."""
    import pandas as pd
    from pyCATHY import cathy_tools
    cwd = os.getcwd()
    prj = os.path.basename(run_dir)
    simu = cathy_tools.CATHY(dirName=os.path.dirname(run_dir), prj_name=prj, version=VERSIONS[sc.version], verbose=False)
    simu.update_prepo_inputs()
    # --- vegetation raster ------------------------------------------------------------------------------------
    vm = VEG_MAPS[sc.veg_map] if sc.veg_map else np.ones((20, 20), dtype=int)
    simu.update_veg_map(vm)
    nveg = int(vm.max())
    # --- soil: pass BOTH FP_map and SPP_map explicitly, otherwise update_soil() calls set_SOIL_defaults() and
    #     silently resets PMIN to -5 m (pyCATHY pitfall found while writing this module)
    SPP = simu.set_SOIL_defaults(SPP_map_default=True)
    FP = simu.set_SOIL_defaults(FP_map_default=True, nveg=nveg)
    FP["PCREF"], FP["PCWLT"], FP["ZROOT"] = sc.pcref, sc.pcwlt, sc.zroot
    if sc.version == "scf":
        FP["SCF"] = list(sc.scf_list) if sc.scf_list else [sc.scf] * nveg
    else:
        simu.soil["SCF"] = float(sc.scf)           # classic: one global scalar on the IPEAT line
    simu.soil["PMIN"] = float(sc.pmin)
    simu.update_soil(FP_map=FP, SPP_map=SPP, shellprint_update=False)
    # --- forcing, initial condition, seepage-face file (sfbc ends at 7 d in the shipped deck -> extend) -----------
    t = [int(round(a)) for a, _ in sc.forcing]
    v = [b for _, b in sc.forcing]
    simu.update_atmbc(HSPATM=1, IETO=1, time=t, netValue=v, shellprint_update=False)
    if sc.wt_height is None:
        simu.update_ic(INDP=0, IPOND=0, pressure_head_ini=float(sc.psi0), shellprint_update=False)
    else:
        simu.update_ic(INDP=3, IPOND=0, WTPOSITION=float(sc.wt_height), shellprint_update=False)
    simu.update_sfbc(time=[0, 100000000], no_flow=True)
    # --- parm ---------------------------------------------------------------------------------------------------
    tmax = sc.days * DAY
    prt = [int(round(x)) for x in np.linspace(tmax / sc.nprt, tmax, sc.nprt)]
    # IPRT1=2 is essential: update_veg_map()/update_soil() trigger a mesh-building run (IPRT1=3) and leave parm in that
    # mode, so a plain run_processor() afterwards would only rebuild the mesh in 0 s and report stale outputs.
    simu.update_parm(IPRT1=2, ISIMGR=sc.isimgr, DTMAX=float(sc.dtmax), TMAX=float(tmax), VTKF=0, TIMPRTi=prt,
                     shellprint_update=False)
    os.chdir(cwd)
    return simu


def run_scenario_inprocess(sc: Scenario, out_json=None):
    """Create the run project from the compiled master, configure it with CATHY, run it. Returns the run dir."""
    from pyCATHY import cathy_tools
    assert 1 <= sc.nprt <= NMAX_PRT, "compiled MAXPRT=9"
    run_dir = os.path.join(WORK, "runs", sc.name)
    shutil.rmtree(run_dir, ignore_errors=True)
    os.makedirs(os.path.dirname(run_dir), exist_ok=True)
    shutil.copytree(master_dir(sc.version), run_dir)
    # NOTE: keep grid3d / grid2d.exp (see below) but drop stale RESULT files copied from the master (the repo ships old outputs)
    _outd = os.path.join(run_dir, "output")
    for f in os.listdir(_outd):
        if f not in ("grid3d", "grid2d.exp"):
            p = os.path.join(_outd, f); shutil.rmtree(p) if os.path.isdir(p) else os.remove(p)
    # NOTE: keep output/ as copied from the master. It holds grid3d/grid2d.exp from the mesh-building run; if they are
    # missing CATHY.run_processor() silently re-runs the mesh build (IPRT1=3) and never simulates. Every run starts from
    # a clean copy of the master, and CATHY rewrites all result files.
    for f in ("fort.666", "fort.777"):
        try: os.remove(os.path.join(run_dir, f))
        except OSError: pass
    simu = configure(sc, run_dir)
    cwd = os.getcwd()
    t0 = time.time()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        simu.run_processor(recompile=False, IPRT1=2)
    os.chdir(cwd)
    return run_dir, time.time() - t0, [str(x.message)[:300] for x in w]


# ----------------------------------------------------------------------------------------------------------------
# outputs
# ----------------------------------------------------------------------------------------------------------------
def _table(path, ncol=None):
    rows = []
    for l in open(path):
        if l.lstrip().startswith("#"): continue
        try: rows.append([float(x) for x in l.split()])
        except ValueError: continue
    rows = [r for r in rows if r]
    n = ncol or min(len(r) for r in rows)
    return np.array([r[:n] for r in rows])


def read_nodal(path, nnod=441, nlay=16):
    """-> (ntimes, nlay, nnod) array and times. (pyCATHY's read_outputs('psi'/'sw') returns the same content as DataFrames.)"""
    if not os.path.exists(path): return None, None
    blocks, times, cur = [], [], []
    for l in open(path):
        if "TIME" in l and "NSTEP" in l:
            if cur: blocks.append(np.array(cur))
            cur = []; times.append(float(l.split()[1])); continue
        try: cur += [float(x) for x in l.split()]
        except ValueError: continue
    if cur: blocks.append(np.array(cur))
    keep = [(t, b) for t, b in zip(times, blocks) if len(b) == nnod * nlay]
    if not keep: return None, None
    return np.array([b.reshape(nlay, nnod) for _, b in keep]), np.array([t for t, _ in keep])


def collect(run_dir, sc: Scenario):
    """Time series and metrics of one finished run."""
    o = {}
    a = _table(os.path.join(run_dir, "output", "hgatmsf"), 10)          # pyCATHY: out_CT.read_hgatmsf
    o["t"], o["pot"], o["act"] = a[:, 2], a[:, 3], a[:, 4]                 # m3/s ; evaporation = negative
    tr = _table(os.path.join(run_dir, "fort.666"), 2)
    o["tr_t"], o["tr"] = tr[:, 0], tr[:, 1]                                # m3/s ; uptake = positive
    m = _table(os.path.join(run_dir, "output", "mbeconv"))                 # pyCATHY: simu.read_outputs('mbeconv')
    o["mb_cum"] = m[:, -2]                                                 # CUM. MBE (m3) -- EXCLUDES root uptake
    o["V0"] = float(open(os.path.join(run_dir, "output", "mbeconv")).readline().split("=")[1])
    dc = _table(os.path.join(run_dir, "output", "dtcoupling"))             # pyCATHY: simu.read_outputs('dtcoupling')
    o["nsteps"], o["backsteps"], o["nl_iter"] = int(dc[-1, 0]), int(dc[-1, 3]), int(dc[:, 5].sum())
    o["psi"], o["psi_t"] = read_nodal(os.path.join(run_dir, "output", "psi"))
    o["sw"], _ = read_nodal(os.path.join(run_dir, "output", "sw"))
    return o


def metrics(o, sc: Scenario):
    t, tt = o["t"], o["tr_t"]
    a = mm_day(o["act"]); ev = np.maximum(-a, 0.0); tr = mm_day(o["tr"])
    after = t > 0.3 * DAY
    out = dict(name=sc.name, version=sc.version, scf=sc.scf_list or sc.scf, pmin=sc.pmin, zroot=sc.zroot, days=sc.days,
               E=float(np.trapezoid(ev, t) / DAY), T=float(np.trapezoid(tr, tt) / DAY),
               phantom_mm=float(np.trapezoid(np.where(after, np.maximum(a, 0), 0), t) / DAY),
               budget_resid_pct_V0=float(abs(o["mb_cum"][-1] - np.trapezoid(o["tr"], tt)) / o["V0"] * 100),
               nsteps=o["nsteps"], backsteps=o["backsteps"], nl_iter=o["nl_iter"],
               t_end_days=float(t[-1] / DAY))
    pet = 5.0
    w = (t > 0.3 * DAY) & (t < 1.0 * DAY)
    out["realized_scf"] = float(1.0 - ev[w].mean() / pet) if (w.any() and sc.tags.get("suite") == "hetero") else None   # 5 mm/d forcing, before pinning
    if o["psi"] is not None:
        out["psi_surf_min"] = float(o["psi"][1:, 0].min())       # excludes t = 0 (imposed IC)
        out["psi_surf_final_med"] = float(np.median(o["psi"][-1, 0]))
    tot = ev + np.interp(t, tt, tr)
    out["max_ET_over_PET_after_pulse"] = float(tot[after].max() / pet)
    out["E_day"] = {d: float(np.interp(d * DAY, t, ev)) for d in (1, 3, 5) if t[-1] >= d * DAY}
    out["T_day"] = {d: float(np.interp(d * DAY, tt, tr)) for d in (1, 3, 5) if tt[-1] >= d * DAY}
    end = ev < 0.95 * (1 - (sc.scf if not sc.scf_list else 0)) * pet
    idx = np.where(after & end)[0]
    out["stage1_end_day"] = float(t[idx[0]] / DAY) if len(idx) and (1 - sc.scf) > 0 and not sc.scf_list else None
    return out


def run_one(sc: Scenario):
    """Run in a child process with a wall-clock cap; returns metrics dict (+ completed flag)."""
    os.makedirs(os.path.join(WORK, "specs"), exist_ok=True)
    spec = os.path.join(WORK, "specs", sc.name + ".json")
    json.dump(asdict(sc), open(spec, "w"))
    t0 = time.time(); timed_out = False
    try:
        pr = subprocess.run([sys.executable, os.path.abspath(__file__), "_child", spec], timeout=sc.timeout,
                            capture_output=True, text=True)
        rc, err = pr.returncode, pr.stderr[-600:]
    except subprocess.TimeoutExpired:
        timed_out, rc, err = True, -9, "timeout"
    run_dir = os.path.join(WORK, "runs", sc.name)
    res = dict(name=sc.name, exit=rc, wall_s=round(time.time() - t0), completed=False)
    try:
        o = collect(run_dir, sc)
        res.update(metrics(o, sc))
        res["completed"] = bool((not timed_out) and rc == 0 and o["t"][-1] >= sc.days * DAY * 0.999)
        np.savez_compressed(os.path.join(RESULTS, sc.name + ".npz"), t=o["t"], act=o["act"], tr_t=o["tr_t"], tr=o["tr"],
                            psi_t=o["psi_t"], psi=o["psi"][:, :9, :], mb_cum=o["mb_cum"])
    except Exception as e:                                   # run produced no usable output
        res["error"] = f"{type(e).__name__}: {e}"; res["stderr"] = err
    for junk in ("stdout.txt",): pass
    shutil.rmtree(os.path.join(run_dir, "vtk"), ignore_errors=True)
    for f in ("hgflag",):                                    # can be >100 MB in the phantom-infiltration cases
        try: os.remove(os.path.join(run_dir, "output", f))
        except OSError: pass
    return res


# ----------------------------------------------------------------------------------------------------------------
# suites (each reproduces one block of the study)
# ----------------------------------------------------------------------------------------------------------------
def suite_quick():
    """1-day runs: the invariants of the old smoke tests, now evaluated from CATHY() runs."""
    F = pulses_forcing(((0.0, 0.05, 10.0),), 5.0, 1.0)
    def S(v, s, p): return Scenario(name=f"q_{v}_scf{s:g}_pmin{p:g}".replace("+", ""), version=v, scf=s, pmin=p, days=1.0,
                                    forcing=F, timeout=600, tags=dict(suite="quick"))
    return [S("scf", 0.0, -150.0), S("scf", 0.3, -150.0), S("scf", 1.0, -1e35), S("classic", 0.0, -150.0),
            S("classic", 1.0, -1e35), S("scf", 0.0, -5.0), S("scf", 0.0, -1e4), S("scf", 1.0, -5.0)]


def suite_matrix():
    F = pulses_forcing(((0.0, 0.25, 25.0), (4.0, 4.1, 8.0)), 5.0, 6.0)
    runs = [("classic", s, p) for s in (0.0, 1.0) for p in (-5.0, -1e35)]
    for p in (-150.0, -1e4, -5.0, -1e35):
        for s in ((0.0, 0.3, 1.0, 0.1) if p == -150.0 else (0.0, 0.3, 1.0)):
            runs.append(("scf", s, p))
    return [Scenario(name=f"m_{v}_scf{s:g}_pmin{p:g}".replace("+", ""), version=v, scf=s, pmin=p, days=6.0, forcing=F,
                     timeout=660, tags=dict(suite="matrix")) for v, s, p in runs]


def suite_drydown():
    F = pulses_forcing(((0.0, 0.25, 25.0),), 5.0, 12.0)
    return [Scenario(name=f"d_scf{s:g}_pmin{p:g}", version="scf", scf=s, pmin=p, days=12.0, forcing=F,
                     tags=dict(suite="drydown")) for s, p in [(1.0, -5.0), (1.0, -150.0), (0.3, -5.0), (0.3, -150.0)]]


def suite_release():
    out = []
    for s in (0.0, 0.3, 0.6):
        out.append(Scenario(name=f"r_A_scf{s:g}", scf=s, days=12.0, forcing=stair_forcing(STAIR_A, 12.0, 25.0),
                            tags=dict(suite="release", stair="A")))
        out.append(Scenario(name=f"r_C_scf{s:g}", scf=s, days=12.0, forcing=stair_forcing(STAIR_C, 12.0, None),
                            wt_height=2.0, tags=dict(suite="release", stair="C")))
    return out


def suite_hetero():
    F = pulses_forcing(((0.0, 0.25, 25.0),), 5.0, 12.0)
    out = [Scenario(name=f"h_{k}", scf_list=[0.0, 1.0], veg_map=k, days=12.0, forcing=F, tags=dict(suite="hetero")) for k in VEG_MAPS]
    out += [Scenario(name=f"h_uniform_{s:g}", scf=s, days=12.0, forcing=F, tags=dict(suite="hetero")) for s in (0.5, 0.475, 0.3)]
    return out


def suite_critique():
    F6 = pulses_forcing(((0.0, 0.25, 25.0), (4.0, 4.1, 8.0)), 5.0, 6.0)
    F12 = pulses_forcing(((0.0, 0.25, 25.0),), 5.0, 12.0)
    out = [Scenario(name=f"e_scf1_zroot{z:g}", scf=1.0, zroot=z, days=6.0, forcing=F6, tags=dict(suite="critique")) for z in (0.006, 0.02, 0.1)]
    out += [Scenario(name=f"c_classic_scf0.3_pmin{p:g}".replace("+", ""), version="classic", scf=0.3, pmin=p, days=6.0, forcing=F6,
                     timeout=900, tags=dict(suite="critique")) for p in (-5.0, -150.0, -1e35)]
    out += [Scenario(name=f"c_classic_scf0.3_pmin{p:g}_12d", version="classic", scf=0.3, pmin=p, days=12.0, forcing=F12,
                     tags=dict(suite="critique")) for p in (-5.0, -150.0)]
    return out


SUITES = dict(quick=suite_quick, matrix=suite_matrix, drydown=suite_drydown, release=suite_release,
              hetero=suite_hetero, critique=suite_critique)


# ----------------------------------------------------------------------------------------------------------------
# invariant checks (replace the pytest smoke tests): return (name, passed, detail)
# ----------------------------------------------------------------------------------------------------------------
def _load(name):
    z = np.load(os.path.join(RESULTS, name + ".npz"))
    return {k: z[k] for k in z.files}


def _cum(z, what):
    if what == "E": return float(np.trapezoid(np.maximum(-mm_day(z["act"]), 0), z["t"]) / DAY)
    return float(np.trapezoid(mm_day(z["tr"]), z["tr_t"]) / DAY)


def checks_quick(R):
    g = lambda v, s, p: f"q_{v}_scf{s:g}_pmin{p:g}".replace("+", "")
    ck = []
    for v, s, p in [("scf", 0.0, -150.0), ("scf", 0.3, -150.0), ("scf", 1.0, -1e35), ("classic", 0.0, -150.0), ("classic", 1.0, -1e35)]:
        r = R[g(v, s, p)]; ok = r["completed"] and r["budget_resid_pct_V0"] < 0.1
        ck.append((f"runs & water budget closes: {g(v, s, p)}", ok, f"completed={r['completed']} residual={r['budget_resid_pct_V0']:.4f} % of V0"))
    z0, z1 = _load(g("scf", 0.0, -150.0)), _load(g("scf", 1.0, -1e35))
    ck.append(("SCF=0 -> no transpiration", np.max(mm_day(z0["tr"])) < 1e-3, f"max T = {np.max(mm_day(z0['tr'])):.2e} mm/d"))
    t = z1["t"] / DAY; ev = np.maximum(-mm_day(z1["act"]), 0)
    ck.append(("SCF=1, PMIN off -> no evaporation", ev[t > 0.1].max() < 0.05, f"max E = {ev[t > 0.1].max():.3f} mm/d"))
    z3 = _load(g("scf", 0.3, -150.0)); t = z3["t"] / DAY; tt = z3["tr_t"] / DAY
    E = np.maximum(-mm_day(z3["act"]), 0); T = mm_day(z3["tr"]); w = (t > 0.06) & (t < 0.20); wt = (tt > 0.06) & (tt < 0.20)
    ck.append(("SCF=0.3: E ~ 0.7*PET right after the pulse", abs(E[w].mean() / 5 - 0.7) < 0.05, f"E/PET = {E[w].mean() / 5:.3f}"))
    ck.append(("SCF=0.3: 0.15 PET < T <= 0.3 PET", 0.15 < T[wt].mean() / 5 <= 0.306, f"T/PET = {T[wt].mean() / 5:.3f}"))
    tot = np.interp(t, tt, T) + E
    ck.append(("SCF=0.3: E + T <= PET", tot[t > 0.06].max() <= 5 * 1.02, f"max (E+T)/PET = {tot[t > 0.06].max() / 5:.3f}"))
    for s, p, what in [(0.0, -150.0, "E"), (1.0, -1e35, "T")]:
        a, b = _cum(_load(g("scf", s, p)), what), _cum(_load(g("classic", s, p)), what)
        ck.append((f"SCF build == classic at SCF={s:g} (cumulative {what}, 2 %)", abs(a - b) / b < 0.02, f"scf {a:.3f} mm vs classic {b:.3f} mm"))
    Ec = {p: _cum(_load(g("scf", 0.0, p)), "E") for p in (-5.0, -150.0, -1e4)}
    ck.append(("cumulative E non-decreasing: PMIN -5 -> -150 -> -1e4", Ec[-5.0] <= Ec[-150.0] * 1.001 <= Ec[-1e4] * 1.002, str({k: round(v, 3) for k, v in Ec.items()})))
    for p in (-5.0, -150.0):
        r = R[g("scf", 0.0, p)]
        ck.append((f"surface head never below PMIN={p:g}", r["psi_surf_min"] >= p * (1 + 1e-6), f"min = {r['psi_surf_min']:.4f} m"))
    a, b = _cum(_load(g("scf", 1.0, -1e35)), "T"), _cum(_load(g("scf", 1.0, -5.0)), "T")
    ck.append(("SCF=1: T independent of PMIN (1-day, 5 %)", abs(a - b) / a < 0.05, f"{a:.3f} vs {b:.3f} mm"))
    return ck


def checks_drydown(R):
    ck = []
    for s in (1.0, 0.3):
        r = R[f"d_scf{s:g}_pmin-150"]; ck.append((f"no phantom infiltration at PMIN=PCWLT, SCF={s:g}", r["phantom_mm"] < 0.05, f"{r['phantom_mm']:.3f} mm"))
        r = R[f"d_scf{s:g}_pmin-5"]; ck.append((f"[known behaviour] phantom infiltration at PMIN=-5, SCF={s:g}", r["phantom_mm"] > 0.05, f"{r['phantom_mm']:.3f} mm (expected > 0)"))
    return ck


def checks_release(R):
    ck = []
    for tag, stair, t0 in (("A", STAIR_A, 3.0), ("C", STAIR_C, 2.0)):
        for s in (0.3, 0.6):
            z = _load(f"r_{tag}_scf{s:g}"); t = z["t"] / DAY; E = np.maximum(-mm_day(z["act"]), 0)
            T = np.interp(z["t"], z["tr_t"], mm_day(z["tr"])); pet = pet_of_stair(stair, t); m = t > t0
            ok1 = bool(np.all(E[m] <= (1 - s) * pet[m] * 1.01 + 1e-6)); ok2 = bool(np.all((E + T)[m] <= pet[m] * 1.01 + 1e-6))
            ratio = float(np.max(E[m] / np.maximum((1 - s) * pet[m], 1e-9)))
            ck.append((f"release rule {tag}, SCF={s:g}: E <= (1-SCF)*PET and E+T <= PET", ok1 and ok2, f"max E/allowed = {ratio:.2f}"))
    return ck


def checks_hetero(R):
    ck = []
    for k, (lo, hi) in {"clumped": (0.44, 0.51), "stripes3": (0.27, 0.33), "checker": (-0.01, 0.02)}.items():
        s = R[f"h_{k}"]["realized_scf"]; ck.append((f"raster '{k}' -> realized SCF in [{lo}, {hi}]", lo <= s <= hi, f"{s:.3f}"))
    het, uni = R["h_clumped"], R["h_uniform_0.475"]
    a, b = het["E"] + het["T"], uni["E"] + uni["T"]
    ck.append(("binary clumped ET < 85 % of uniform float ET at equal cover", a < 0.85 * b, f"{a:.1f} vs {b:.1f} mm"))
    return ck


def checks_critique(R):
    ck = []
    for z in (0.006, 0.02, 0.1):
        r = R[f"e_scf1_zroot{z:g}"]; ck.append((f"ZROOT={z:g} m at SCF=1 cannot create bare-soil evaporation", r["E"] < 0.05, f"E = {r['E']:.3f} mm, T = {r['T']:.2f} mm"))
    return ck


CHECKS = dict(quick=checks_quick, drydown=checks_drydown, release=checks_release, hetero=checks_hetero, critique=checks_critique)


# ----------------------------------------------------------------------------------------------------------------
def run_suite(name, only_missing=True):
    os.makedirs(RESULTS, exist_ok=True)
    path = os.path.join(RESULTS, f"{name}.json")
    R = json.load(open(path)) if (only_missing and os.path.exists(path)) else {}
    for sc in SUITES[name]():
        if sc.name in R and R[sc.name].get("exit") in (0, -9) and "error" not in R[sc.name]: continue
        print(time.strftime("%H:%M:%S"), name, sc.name, flush=True)
        R[sc.name] = run_one(sc)
        r = R[sc.name]
        print("   ->", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items() if k in ("completed", "E", "T", "wall_s", "nl_iter", "backsteps", "error")}, flush=True)
        json.dump(R, open(path, "w"), indent=1)
    if name in CHECKS:
        ck = CHECKS[name](R)
        md = [f"## {name}", "", "| check | result | detail |", "|---|---|---|"] + [f"| {n} | {'PASS' if ok else 'FAIL'} | {d} |" for n, ok, d in ck]
        open(os.path.join(RESULTS, f"checks_{name}.md"), "w").write("\n".join(md) + "\n"); print("\n".join(md))
    return R


def compare(ref_dir):
    """Compare the pyCATHY-driven results with the reference numbers of the direct-runner study (relative tol 2 %)."""
    pairs = []   # (our suite, our name, reference json, reference key)
    for nm, r in json.load(open(os.path.join(RESULTS, "matrix.json"))).items() if os.path.exists(os.path.join(RESULTS, "matrix.json")) else []:
        pairs.append((nm, r, "results/summary.json", nm.replace("m_", "", 1)))
    ref = {}
    for f in ("results/summary.json", "results_drydown/summary.json", "results_extra/summary.json", "results_critique/summary.json"):
        p = os.path.join(ref_dir, f)
        if os.path.exists(p): ref[f] = json.load(open(p))
    rows = ["| run | E ours | E ref | T ours | T ref | max rel. diff |", "|---|---|---|---|---|---|"]
    for nm, r, f, key in pairs:
        if f in ref and key in ref[f] and r.get("completed"):
            e0, t0 = ref[f][key]["cumE_mm"], ref[f][key]["cumT_mm"]
            d = max(abs(r["E"] - e0) / max(e0, 1.0), abs(r["T"] - t0) / max(t0, 1.0))
            rows.append(f"| {nm} | {r['E']:.2f} | {e0:.2f} | {r['T']:.2f} | {t0:.2f} | {100 * d:.2f} % |")
    open(os.path.join(RESULTS, "compare_with_reference.md"), "w").write("\n".join(rows) + "\n"); print("\n".join(rows))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["build", "run", "compare", "_child"]); ap.add_argument("arg", nargs="?")
    ap.add_argument("--ref", default=os.path.join(HERE, "..", "harness"), help="folder of the direct-runner results")
    a = ap.parse_args()
    if a.cmd == "build":
        for v in ("scf", "classic"): print("building", v, build_master(v))
    elif a.cmd == "_child":                        # run one scenario (called by run_one)
        sc = Scenario(**json.load(open(a.arg)))
        run_scenario_inprocess(sc)
    elif a.cmd == "run":
        names = list(SUITES) if a.arg == "all" else [a.arg]
        for n in names: run_suite(n)
    elif a.cmd == "compare":
        compare(a.ref)


if __name__ == "__main__":
    main()
