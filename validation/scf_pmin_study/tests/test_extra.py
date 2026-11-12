"""Follow-up tests: (1) how a vegetation raster is turned into SCF classes, heterogeneous vs uniform cover;
(2) the release rule of the pinned (PMIN) surface node when SCF is fractional.

Quick ones run by default (~2 min); slow ones:  pytest -m "release or hetero_long" tests   (~10 min)
"""
import numpy as np, pytest
import cathy_runner as c
from conftest import ROOT

DAY = c.DAY
F_PULSE = lambda days: c.semi_arid_forcing(((0.0, 0.25, 25.0),), 5.0, days)
CLUMPED = [[1] * 10 + [2] * 10 for _ in range(20)]
STRIPES3 = [[1 if (j // 3) % 2 == 0 else 2 for j in range(20)] for _ in range(20)]
CHECKER = [[1 if (i + j) % 2 == 0 else 2 for j in range(20)] for i in range(20)]

def _run(name, days, forcing, **kw):
    spec = c.RunSpec(name=name, version="scf", pmin=-150.0, days=days, forcing=forcing, psi0=-20.0, dtmax=1800.0, nprt=9, **kw)
    return c.run(spec, ROOT, timeout=900)

def realized_scf(r, pet=5.0):
    """Area-weighted mean SCF actually seen by the model: before any PMIN pinning, E = (1-SCF)*PET."""
    t = r["t"]; E = (-c.mm_day(r["act"])).clip(min=0)
    w = (t > 0.3 * DAY) & (t < 1.0 * DAY)
    return 1.0 - float(E[w].mean() / pet)

def total_et(r):
    t = r["t"]; E = (-c.mm_day(r["act"])).clip(min=0)
    return float(np.trapezoid(E, t) / DAY + np.trapezoid(c.mm_day(r["tr"]), r["tr_t"]) / DAY)

# ------------------------------------------------------------------ 1. vegetation raster -> SCF classes (quick)
@pytest.mark.parametrize("name,vmap,lo,hi", [
    ("clumped", CLUMPED, 0.44, 0.51),      # 10 of 21 node columns are pure class 2 (boundary column truncates to class 1)
    ("stripes3", STRIPES3, 0.27, 0.33),
    ("checker", CHECKER, -0.01, 0.02),     # cell-scale mix collapses: node value = mean of neighbours (1.5) -> INT() -> class 1
])
def test_veg_raster_realized_cover(name, vmap, lo, hi):
    r = _run(f"x_{name}_1d", 1.1, F_PULSE(1.1), scf_list=[0.0, 1.0], veg_map=vmap)
    assert r["exit"] == 0 and r["completed"]
    s = realized_scf(r)
    assert lo <= s <= hi, f"{name}: realized mean SCF {s:.3f} not in [{lo}, {hi}]"

def test_uniform_float_scf_is_realized_exactly():
    r = _run("x_uniform05_1d", 1.1, F_PULSE(1.1), scf=0.5)
    assert abs(realized_scf(r) - 0.5) < 0.01

def test_heterogeneous_map_closes_budget():
    r = _run("x_clumped_1d", 1.1, F_PULSE(1.1), scf_list=[0.0, 1.0], veg_map=CLUMPED)
    res = abs(r["mb_cum"][-1] - np.trapezoid(r["tr"], r["tr_t"])) / r["V0"] * 100
    assert res < 0.1

# ------------------------------------------------------------------ 2. slow: binary 0|1 patches vs uniform float SCF at equal realized cover
@pytest.mark.hetero_long
def test_clumped_binary_transpires_less_than_uniform_float_scf():
    het = _run("x_clumped_12d", 12.0, F_PULSE(12.0), scf_list=[0.0, 1.0], veg_map=CLUMPED)
    uni = _run("x_uniform0475_12d", 12.0, F_PULSE(12.0), scf=round(realized_scf(het), 3))
    a, b = total_et(het), total_et(uni)
    assert a < 0.85 * b, f"expected pixel-scale stress to cut total ET by >15 %: binary {a:.1f} mm vs uniform {b:.1f} mm"

# ------------------------------------------------------------------ 3. release rule of the pinned node
STAIR_C = [(0.0, 8.0), (2.0, 4.0), (4.0, 2.5), (6.0, 1.5), (8.0, 1.0), (10.0, 0.6)]
def _stair_forcing(stair, days=12.0, pulse=False):
    f = [(0.0, 25.0 / 1000 / (0.25 * DAY))] if pulse else []
    f += [(d * DAY, -p / 1000 / DAY) for d, p in stair]; f.append(((days + 10) * DAY, f[-1][1])); return f
def _pet(stair, t):
    p = np.zeros_like(t, dtype=float)
    for d, v in stair: p = np.where(t >= d, v, p)
    return p

@pytest.mark.release
@pytest.mark.parametrize("scf", [0.3, 0.6])
def test_pinned_node_never_evaporates_more_than_its_share(scf):
    """SWITCH/ADRSTN release a pinned node only if the diagnosed flux exceeds the FULL ATMPOT (not (1-SCF)*ATMPOT).
    Consequence if the window opened: E would stay above its allotted share (1-SCF)*PET, and E+T could exceed PET.
    Measured: the window does not open for PET >= 0.6 mm/d (soil-limited flux at PMIN=-150 m is <= 0.5 mm/d)."""
    r = _run(f"x_release_wt_scf{scf:g}", 12.0, _stair_forcing(STAIR_C), scf=scf, wt_height=2.0)
    assert r["exit"] == 0 and r["completed"]
    t = r["t"] / DAY; E = (-c.mm_day(r["act"])).clip(min=0); T = np.interp(r["t"], r["tr_t"], c.mm_day(r["tr"])); pet = _pet(STAIR_C, t)
    m = t > 2.0
    assert np.all(E[m] <= (1 - scf) * pet[m] * 1.01 + 1e-6)
    assert np.all((E + T)[m] <= pet[m] * 1.01 + 1e-6)
