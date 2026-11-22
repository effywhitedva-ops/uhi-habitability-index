"""Precalcula pistas de Johnstone+2021 (Mors, con evolucion estelar de Spada+2013, Zenodo 15729101 md5 f76987cf...)
en una malla de masa x percentil de rotacion inicial. Guarda edades [Gyr], Lbol, Lx, Leuv [erg/s] y sus integrales acumuladas."""
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
cum = {}
dt = np.diff(np.concatenate([[0.0], AGES])) * 3.15576e16   # s
for k in ("Lx", "Leuv"):
    cum["C" + k] = np.cumsum(out[k] * dt, axis=-1)            # erg (rectangulos; primer paso desde t=0 con L(1 Myr))
np.savez_compressed(os.path.join(os.path.dirname(__file__), "..", "data", "mors_grid.npz"), masses=MASSES, pcts=PCTS, ages=AGES, **out, **cum)
print("ok")
