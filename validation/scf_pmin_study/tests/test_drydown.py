"""12-day dry-down (one 25 mm pulse, 5 mm/d demand). Slow: run with  pytest -m drydown tests.
Invariant: after the rain pulse, the surface must not *add* water to a stressed soil.  CATHY itself flags the
violation in output/hgflag ('INFILTRATION AT AN AIR DRY NODE WITH EVAPORATIVE POTENTIAL FLUX')."""
import numpy as np, pytest
import cathy_runner as c
from conftest import ROOT

def _run(scf, pmin):
    f = c.semi_arid_forcing(((0.0, 0.25, 25.0),), 5.0, 12.0)
    spec = c.RunSpec(name=f"dry_scf{scf:g}_pmin{pmin:g}", version="scf", scf=scf, pmin=pmin, days=12.0, forcing=f,
                     psi0=-20.0, dtmax=1800.0, nprt=9)
    return c.run(spec, ROOT, timeout=900)

def phantom_infiltration_mm(r):
    t, a = r["t"], r["act"] / c.AREA * c.DAY * 1000.0
    m = t > 0.3 * c.DAY
    return float(np.trapezoid(np.where(m, np.maximum(a, 0.0), 0.0), t) / c.DAY)

@pytest.mark.drydown
@pytest.mark.parametrize("scf", [1.0, 0.3])
def test_no_phantom_infiltration_when_pmin_at_or_below_pcwlt(scf):
    r = _run(scf, -150.0)                       # PMIN == PCWLT
    assert r["exit"] == 0 and r["completed"]
    assert phantom_infiltration_mm(r) < 0.05

@pytest.mark.drydown
@pytest.mark.xfail(strict=False, reason="Known behaviour: PMIN (-5 m) less negative than PCWLT (-150 m) pins the surface node wetter "
                   "than the dried root zone -> phantom infiltration (measured 1.50 mm at SCF=1, 0.32 mm at SCF=0.3 over 12 d)")
@pytest.mark.parametrize("scf", [1.0, 0.3])
def test_no_phantom_infiltration_when_pmin_above_pcwlt(scf):
    r = _run(scf, -5.0)
    assert phantom_infiltration_mm(r) < 0.05
