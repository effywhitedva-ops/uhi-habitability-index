r"""Precompute the Mors stellar-activity grid used by xuv_history.py.

Johnstone, Bartel & Gudel (2021, A&A 649, A96) Mors tracks, with Spada et al. (2013) stellar
evolution models from Zenodo record 15729101 (archive md5 in docs/assumptions_v6_es.md), are
sampled on 0.10-1.20 M_sun in 0.05 steps and initial-rotation percentiles 5-95. Ages are
logarithmic from 1 Myr to 0.1 Gyr and linear from 0.1 to 13 Gyr.

Output: data/mors_grid.npz with masses (M_sun), pcts (percentiles), ages (Gyr), Lbol/Lx/Leuv
(erg s^-1), and cumulative CLx/CLeuv (erg). Track values outside their age range are held at their
endpoints, and integration uses rectangles from t=0 with the first sampled luminosity; both are
simple numerical choices. This build script requires ext/Mors and ext/fs255_grid, which are not
included in the release; the generated grid is included so users do not need to rebuild it."""
import os, sys, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "ext", "Mors", "src"))
os.environ.setdefault("STELLARMODELS", os.path.join(os.path.dirname(__file__), "..", "ext", "fs255_grid") + "/")
import Mors as mors
MASSES = np.round(np.arange(0.10, 1.2001, 0.05), 2)
PCTS = np.array([5, 15, 25, 35, 50, 65, 75, 85, 95], float)
AGES = np.concatenate([np.geomspace(0.001, 0.1, 40, endpoint=False), np.linspace(0.1, 13.0, 260)])   # Gyr
out = {k: np.zeros((len(MASSES), len(PCTS), len(AGES))) for k in ("Lbol", "Lx", "Leuv")}
for i, m in enumerate(MASSES):
    for j, p in enumerate(PCTS):
        s = mors.Star(Mstar=float(m), percentile=float(p))
        A = s.Tracks["Age"] / 1e3
        for k in out:
            out[k][i, j] = np.interp(AGES, A, s.Tracks[k], left=s.Tracks[k][0], right=s.Tracks[k][-1])
        print(m, p, flush=True)
# Cumulative energies are later divided by 4 pi a^2 to obtain lifetime fluence in erg cm^-2.
cum = {}
dt = np.diff(np.concatenate([[0.0], AGES])) * 3.15576e16   # s
for k in ("Lx", "Leuv"):
    cum["C" + k] = np.cumsum(out[k] * dt, axis=-1)            # erg (rectangle rule; first step from t=0 using L(1 Myr))
np.savez_compressed(os.path.join(os.path.dirname(__file__), "..", "data", "mors_grid.npz"), masses=MASSES, pcts=PCTS, ages=AGES, **out, **cum)
print("ok")
