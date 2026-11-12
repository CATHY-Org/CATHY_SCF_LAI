import os, sys, pytest
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "harness"))
import cathy_runner as c

ROOT = os.environ.get("SMOKE_ROOT", os.path.join(HERE, "_runs"))
# 1-day smoke forcing: 10 mm pulse over the first 1.2 h, then 5 mm/d potential ET
SMOKE_DAYS = 1.0
SMOKE_FORCING = c.semi_arid_forcing(pulses=((0.0, 0.05, 10.0),), pet_mm_day=5.0, days=SMOKE_DAYS)
PET = 5.0
_cache = {}

def get_run(version, scf, pmin, timeout=600):
    key = (version, scf, pmin)
    if key not in _cache:
        name = f"{version}_scf{scf:g}_pmin{pmin:g}".replace("+", "")
        spec = c.RunSpec(name=name, version=version, scf=scf, pmin=pmin, days=SMOKE_DAYS,
                         forcing=SMOKE_FORCING, psi0=-20.0, dtmax=1800.0, nprt=9)
        _cache[key] = (spec, c.run(spec, ROOT, timeout=timeout))
    return _cache[key]

@pytest.fixture(scope="session")
def run():
    if not (os.path.exists(os.path.join(c.BIN, "cathy_scf")) and os.path.exists(os.path.join(c.BIN, "cathy_classic"))):
        pytest.skip("binaries missing: run harness/build.sh first")
    return get_run

def series(r):
    """evaporation (mm/d, >=0), transpiration (mm/d, >=0), time (d) for the evaporation and transpiration records"""
    ev = (-c.mm_day(r["act"])).clip(min=0.0)
    return r["t"] / c.DAY, ev, r["tr_t"] / c.DAY, c.mm_day(r["tr"])

def cum(t, y): 
    import numpy as np
    return float(np.trapezoid(y, t))   # mm  (t in days, y in mm/d)

def surface_psi_min(r):
    return float(r["psi"][1:, 0, :].min())   # drop t=0 block (the imposed initial condition)

def pytest_configure(config):
    config.addinivalue_line("markers", "stress: slow numerical-stress test (PMIN disabled); run with -m stress")


def budget_residual_pct(r):
    """|CUM.MBE - integral(QTRAN)| / initial storage, in %.  mbeconv does NOT include the root-uptake sink,
    so the raw CUM.MBE equals the transpired volume when SCF>0; closing the budget needs QTRAN added back."""
    import numpy as np
    return float(abs(r["mb_cum"][-1] - np.trapezoid(r["tr"], r["tr_t"])) / r["V0"] * 100.0)
