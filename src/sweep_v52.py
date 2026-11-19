"""v5.2: barrido 1-D (uno a la vez) del veto hibrido: F_rg y frontera Teff del veto (teff_rg).
Variante secundaria: frontera M-sincronica completa (teff_md = teff_rg). Corre local, 10 000 draws, semilla 42."""
import os, sys, json, time, hashlib
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np, pandas as pd
from multiprocessing import Pool
CAT = "pscomppars.csv"
CONFIGS = [("base", {})] + [(f"Frg_{v}", {"F_rg": float(v)}) for v in (275, 280, 286, 290, 295)] \
        + [(f"Trg_{v}", {"teff_rg": float(v)}) for v in (3700, 3750, 3850, 3900)] \
        + [(f"Tmd_{v}", {"teff_rg": float(v), "teff_md": float(v)}) for v in (3700, 3750, 3850, 3900)]

def job(cfg):
    tag, ov = cfg
    import uhi_code as u, anchors_v5 as A
    assert u.V3["F_rg"] == 282.0 and u.V3["teff_rg"] == 3800.0 and u.V3["teff_md"] == 3800.0, "estado V3 contaminado"
    u.V3.update(ov)   # proceso nuevo por tarea (maxtasksperchild=1): sin fuga de overrides entre configuraciones
    print(tag, "F_rg", u.V3["F_rg"], "teff_rg", u.V3["teff_rg"], "teff_md", u.V3["teff_md"], flush=True)
    t = time.time(); df = pd.read_csv(CAT)
    res = u.run(df, n_draws=10000, seed=42); res.to_csv(f"sweep_v52/{tag}.csv", index=False)
    an = A.run_anchors(f"sweep_v52/{tag}.csv", 10000); an.to_csv(f"sweep_v52/{tag}_anchors.csv", index=False)
    return tag, round(time.time() - t, 1)

if __name__ == "__main__":
    print("md5", hashlib.md5(open(CAT, "rb").read()).hexdigest(), flush=True)
    os.makedirs("sweep_v52", exist_ok=True)
    todo = [c for c in CONFIGS if not os.path.exists(f"sweep_v52/{c[0]}_anchors.csv")]
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 4, maxtasksperchild=1) as p:
        for tag, s in p.imap_unordered(job, todo):
            print(tag, s, flush=True)
