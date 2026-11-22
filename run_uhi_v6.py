#!/workspace/uhi_v6/.venv/bin/python
"""UHI v6 Fase 1: corrida Monte Carlo completa (semilla 42, snapshot md5 f2f70b24e67a9054940267e8081c81e8).
  python run_uhi_v6.py --thr 5 --variant main          # Rocky, Hycean y combinado
Variantes: main | compat | teff (D1 logistica en Teff) | bt3d (costa 3D Berta-Thompson+2026) | jackson (D2) |
           madhu (D4 Hycean optimista) | only_clim | only_xuv | only_d3 | single (D1+D2+D3 sin clases)"""
import argparse, hashlib, os, sys, time, json, warnings
warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(ROOT, "src"))
import pandas as pd, numpy as np
import uhi_code as u, uhi_code_v6 as v6, anchors_v5 as A5
SNAP = os.path.join(ROOT, "data", "pscomppars_2026-11-18.csv")
VAR = {"main": {}, "compat": dict(v6.CFG_COMPAT), "teff": dict(climate_mode="teff"), "bt3d": dict(shore_model="bt3d"),
       "jackson": dict(xuv_star="jackson"), "madhu": dict(hycean_climate="madhu"),
       "only_clim": dict(v6.CFG_COMPAT, climate_mode="rot"), "only_xuv": dict(v6.CFG_COMPAT, xuv=True),
       "only_d3": dict(v6.CFG_COMPAT, d3=True), "single": dict(v6.CFG_V6, branches=False),
       # --- pruebas de sensibilidad al prior de clase (2026-11-20) ---
       "pA1_rocoso": dict(class_prior=(2.0, 1.0, 1.0)), "pA2_volatil": dict(class_prior=(0.5, 1.0, 1.0)),
       "pA3_ck17": {}, "pA4_empM": {}, "pA5_empM_P5": {},
       "pB1_msini_min": dict(msini="min"), "pB2_msini_howard": dict(msini="howard"),
       "pC1_mr_ck17": dict(mr="ck"), "pD1_linea_earthlike": dict(rock_line="earthlike"), "pD2_linea_h2o100": dict(rock_line="h2o100")}
# sobreescrituras de parametros por variante (ruta de masa P_R(M) = sigmoide((ln M_v - ln M)/s))
VARP = {"pA1_rocoso": dict(mv_lo=5.0, mv_hi=10.0), "pA2_volatil": dict(mv_lo=1.5, mv_hi=3.0),
        "pA3_ck17": dict(mv_lo=2.04, mv_hi=2.04),                       # Chen & Kipping 2017: transicion Terran->Neptunian 2.04 M_E
        "pA4_empM": dict(mv_lo=4.10, mv_hi=4.10, mv_s=0.40),            # ajuste logistico al archivo (dens_ok, M*<0.6, N=101): src/fit_empirical_prior.py
        "pA5_empM_P5": dict(mv_lo=2.83, mv_hi=2.83, mv_s=0.70)}         # idem con P>5 d (N=44)
ap = argparse.ArgumentParser(); ap.add_argument("--thr", type=float, default=5.0); ap.add_argument("--variant", default="main")
ap.add_argument("--draws", type=int, default=10000); ap.add_argument("--outdir", default="results")
ap.add_argument("--legacy", action="store_true", help="fixes v6.1 apagados (regresion contra v6.0)"); a = ap.parse_args()
assert hashlib.md5(open(SNAP, "rb").read()).hexdigest() == "f2f70b24e67a9054940267e8081c81e8"
u.V3["snr_thr"] = a.thr
P = {**u.V3, **v6.V6, **VARP.get(a.variant, {})}
cfg = dict(v6.CFG_V6); cfg.update(VAR[a.variant])
if a.legacy:
    cfg["rs_floor"] = False; import xuv_history as _X; _X.INNES_S10[0] = 0.0852
OUT = a.outdir if os.path.isabs(a.outdir) else os.path.join(ROOT, a.outdir)
tag = f"{a.variant}_{int(a.thr)}s"
t = time.time(); df = pd.read_csv(SNAP)
res = v6.run_v6(df, a.draws, seed=42, cfg=cfg, P=P, progress=lambda i, N: print(tag, i, N, f"{time.time()-t:.0f}s", flush=True))
os.makedirs(OUT, exist_ok=True)
summ = dict(tag=tag, cfg=cfg, thr=a.thr, draws=a.draws, seconds=round(time.time() - t), version=v6.VERSION)
for b, r in res.items():
    r.to_csv(os.path.join(OUT, f"uhi_v6_{b}_{tag}.csv"), index=False)
    summ[b] = dict(scored=len(r), alive=int((r.P_alive > 0).sum()))
an = v6.run_v6(A5.synth_rows(), a.draws, seed=42, cfg=cfg, P=P, joint=False)
for b, r in an.items():
    r.to_csv(os.path.join(OUT, f"anchors_v6_{b}_{tag}.csv"), index=False)
summ["alive_P001"] = {b: int((r.P_alive >= 0.01).sum()) for b, r in res.items()}
json.dump(summ, open(os.path.join(OUT, f"summary_{tag}.json"), "w"), indent=1, default=str)
print(json.dumps(summ, default=str), flush=True)
