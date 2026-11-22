#!/usr/bin/env python3
r"""Run the full UHI v6.1 Monte Carlo on the frozen NASA Exoplanet Archive snapshot.

Example: python run_uhi_v6.py --variant main --thr 5 --outdir results/v6.1

Arguments: --thr is the SNR detection threshold in sigma (5 production, 3 sensitivity scenario);
--variant selects a model configuration; --draws sets draws per planet (10,000 for release);
--outdir is relative to the repository root unless absolute; --legacy disables v6.1 fixes for
regression against v6.0 (R* floor off and the Innes M-star 10-bar preprint value 0.0852).

Input: data/pscomppars_2026-11-18.csv, whose md5 must be f2f70b24e67a9054940267e8081c81e8.
The seed is fixed to 42. Output files are uhi_v6_<branch>_<variant>_<thr>s.csv for rocky, Hycean and
combined; anchors_v6_<branch>_<tag>.csv for synthetic Earth/Venus/Mars controls; and summary_<tag>.json
with configuration, runtime, number scored, number with P_alive>0, and number with P_alive>=0.01.
Branch CSV ranks are pooled across facilities; filter transiting==1 for JWST and 0 for ELT/PCS lists.

Variants: main; compat (v5.2 physics); teff (Teff-based synchronous climate transition); bt3d
(Berta-Thompson et al. 2026 bolometric shoreline); jackson (Jackson et al. 2012 XUV); madhu
(Madhusudhan et al. 2021 Hycean climate); only_clim, only_xuv, only_d3 (single D1/D2/D3 module);
single (D1+D2+D3 without classes); plus class-prior tests pA*, msini tests pB*, M-R test pC1, and
composition-line tests pD* (assumptions, prior tests (a)-(d)). Non-main variants are sensitivity checks.

Note: the first line is a local interpreter path from the development machine. Run with
"python run_uhi_v6.py" after installing requirements.txt, rather than executing the file directly."""
import argparse, hashlib, os, sys, time, json, warnings
warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(ROOT, "src"))
import pandas as pd, numpy as np
import uhi_code as u, uhi_code_v6 as v6, anchors_v5 as A5
SNAP = os.path.join(ROOT, "data", "pscomppars_2026-11-18.csv")
# Model-switch overrides for each variant; see uhi_code_v6.CFG_V6 for keys and defaults.
VAR = {"main": {}, "compat": dict(v6.CFG_COMPAT), "teff": dict(climate_mode="teff"), "bt3d": dict(shore_model="bt3d"),
       "jackson": dict(xuv_star="jackson"), "madhu": dict(hycean_climate="madhu"),
       "only_clim": dict(v6.CFG_COMPAT, climate_mode="rot"), "only_xuv": dict(v6.CFG_COMPAT, xuv=True),
       "only_d3": dict(v6.CFG_COMPAT, d3=True), "single": dict(v6.CFG_V6, branches=False),
       # --- class-prior sensitivity tests (2026-11-20) ---
       "pA1_rocoso": dict(class_prior=(2.0, 1.0, 1.0)), "pA2_volatil": dict(class_prior=(0.5, 1.0, 1.0)),
       "pA3_ck17": {}, "pA4_empM": {}, "pA5_empM_P5": {},
       "pB1_msini_min": dict(msini="min"), "pB2_msini_howard": dict(msini="howard"),
       "pC1_mr_ck17": dict(mr="ck"), "pD1_linea_earthlike": dict(rock_line="earthlike"), "pD2_linea_h2o100": dict(rock_line="h2o100")}
# per-variant parameter overrides (mass route P_R(M) = sigmoid((ln M_v - ln M)/s))
VARP = {"pA1_rocoso": dict(mv_lo=5.0, mv_hi=10.0), "pA2_volatil": dict(mv_lo=1.5, mv_hi=3.0),
        "pA3_ck17": dict(mv_lo=2.04, mv_hi=2.04),                       # Chen & Kipping 2017: Terran->Neptunian transition at 2.04 M_E
        "pA4_empM": dict(mv_lo=4.10, mv_hi=4.10, mv_s=0.40),            # logistic fit to the archive (dens_ok, M*<0.6, N=101): src/fit_empirical_prior.py
        "pA5_empM_P5": dict(mv_lo=2.83, mv_hi=2.83, mv_s=0.70)}         # same fit restricted to P>5 d (N=44)
ap = argparse.ArgumentParser(); ap.add_argument("--thr", type=float, default=5.0); ap.add_argument("--variant", default="main")
ap.add_argument("--draws", type=int, default=10000); ap.add_argument("--outdir", default="results")
ap.add_argument("--legacy", action="store_true", help="fixes v6.1 apagados (regresion contra v6.0)"); a = ap.parse_args()
# Refuse to run if the input catalogue differs from the frozen snapshot used for the release.
assert hashlib.md5(open(SNAP, "rb").read()).hexdigest() == "f2f70b24e67a9054940267e8081c81e8"
u.V3["snr_thr"] = a.thr
P = {**u.V3, **v6.V6, **VARP.get(a.variant, {})}
cfg = dict(v6.CFG_V6); cfg.update(VAR[a.variant])
# Legacy mode modifies an imported module constant only within this process, to reproduce v6.0 for the
# regression check in results/v6.1_legacy_check.
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
# Synthetic controls are rerun with a fresh seed-42 generator. joint=False because they do not participate in
# catalogue ranking.
an = v6.run_v6(A5.synth_rows(), a.draws, seed=42, cfg=cfg, P=P, joint=False)
for b, r in an.items():
    r.to_csv(os.path.join(OUT, f"anchors_v6_{b}_{tag}.csv"), index=False)
summ["alive_P001"] = {b: int((r.P_alive >= 0.01).sum()) for b, r in res.items()}
json.dump(summ, open(os.path.join(OUT, f"summary_{tag}.json"), "w"), indent=1, default=str)
print(json.dumps(summ, default=str), flush=True)
