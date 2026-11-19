#!/usr/bin/env python3
"""UHI v5.2 RC: script maestro reproducible. Corre la simulacion completa, las anclas, el barrido del veto,
el benchmark y genera TODAS las tablas, estadisticas y figuras sin pasos manuales.

  python run_uhi_v52.py                 # todo (5σ + 3σ + anclas + barrido 1-D + PHL + tablas + figuras + MRT)
  python run_uhi_v52.py --skip-sweep    # sin el barrido (14 corridas extra)
  python run_uhi_v52.py --no-phl        # sin descargar el Habitable Worlds Catalog (solo ESI calculado)

Usa el snapshot fechado data/pscomppars_2026-11-18.csv (md5 f2f70b24e67a9054940267e8081c81e8) y semilla 42."""
import argparse, hashlib, os, shutil, subprocess, sys, time, urllib.request
ROOT = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(ROOT, "src")
SNAP = os.path.join(ROOT, "data", "pscomppars_2026-11-18.csv")
MD5 = "f2f70b24e67a9054940267e8081c81e8"
PHL_URL = "https://www.hpcf.upr.edu/~abel/phl/hwc/data/hwc.csv"   # PHL @ UPR Arecibo, Habitable Worlds Catalog (no se redistribuye)
ap = argparse.ArgumentParser()
ap.add_argument("--out", default=os.path.join(ROOT, "output"))
ap.add_argument("--draws", type=int, default=10000)
ap.add_argument("--workers", type=int, default=4)
ap.add_argument("--skip-sweep", action="store_true")
ap.add_argument("--no-phl", action="store_true")
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True); os.chdir(a.out)
sys.path.insert(0, SRC); env = dict(os.environ, PYTHONPATH=SRC + os.pathsep + os.environ.get("PYTHONPATH", ""))
md5 = hashlib.md5(open(SNAP, "rb").read()).hexdigest()
assert md5 == MD5, f"snapshot md5 {md5} != {MD5}"
shutil.copy(SNAP, "pscomppars.csv"); print("snapshot md5 OK", md5, flush=True)
if not os.path.exists("sens_out_v5"):
    shutil.copytree(os.path.join(ROOT, "sens_out_v5"), "sens_out_v5")
import warnings; warnings.filterwarnings("ignore")
import pandas as pd, uhi_code as u, anchors_v5 as A
df = pd.read_csv("pscomppars.csv")
for thr, tag in [(5.0, ""), (3.0, "_3sigma")]:
    u.V3["snr_thr"] = thr; t = time.time()
    res = u.run(df, n_draws=a.draws, seed=42); res.to_csv(f"uhi_results_v52{tag}.csv", index=False)
    an = A.run_anchors(f"uhi_results_v52{tag}.csv", a.draws); an.to_csv(f"anchors_v52{tag}.csv", index=False)
    print(f"{thr:g}σ: puntuados {len(res)}, vivos {(res.P_alive > 0).sum()}, {time.time() - t:.0f} s", flush=True)
u.V3["snr_thr"] = 5.0
if not a.skip_sweep:
    subprocess.run([sys.executable, "-W", "ignore", os.path.join(SRC, "sweep_v52.py"), str(a.workers)], check=True, env=env)
if not a.no_phl:
    try:
        os.makedirs("phl", exist_ok=True); urllib.request.urlretrieve(PHL_URL, "phl/hwc.csv")
        print("PHL md5", hashlib.md5(open("phl/hwc.csv", "rb").read()).hexdigest(), "(certificado v5.2: e971b8dcc305a2fb142d2840aaa8ca37)", flush=True)
    except Exception as ex:
        print("PHL no disponible; se usa solo el ESI calculado:", ex, flush=True)
subprocess.run([sys.executable, os.path.join(SRC, "make_tables.py"), "."], check=True, env=env)
subprocess.run([sys.executable, os.path.join(SRC, "make_figures.py"), "."], check=True, env=env)
subprocess.run([sys.executable, os.path.join(SRC, "make_mrt.py"), "uhi_results_v52.csv", "uhi_results_v52.mrt"], check=True, env=env)
print("listo:", a.out)
