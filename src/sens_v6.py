"""Sobol (Saltelli, SALib) de UHI v6 Fase 1 a umbral SNR FIJO (escenario; D5): python sens_v6.py THR N WORKERS OUTDIR
Evalua el MISMO codigo de produccion (simulate_planet_v6) con 100 draws por planeta sobre los 1213 planetas de la
muestra de sensibilidad de v5 (sens_planets_v5.txt), con numeros aleatorios comunes (semilla por planeta = crc32 del nombre).
Salidas por rama (rocky, hycean, combined): solapamiento top-10 con la corrida base, Spearman vs base, numero de vivos y
UHI_rank de los planetas clave. S1/ST con calc_second_order=False (como v5)."""
import os, sys, time, zlib, json, warnings
warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.stats import spearmanr
import uhi_code as u, uhi_code_v6 as v6

# (nombre, base, lo, hi, descripcion)
PARAMS = [
 # --- v5.2 heredados con mayor S_T medio en el Sobol de v5 (fR_d se excluye: Rocky ya no tiene borde superior de radio) ---
 ("floor_ppm", 20.0, 10.0, 50.0, "S_obs: piso sistematico [ppm]"),
 ("n_H", 2.0, 1.0, 5.0, "S_obs: escalas de altura"),
 ("iwa_hi", 4.0, 3.0, 6.0, "S_obs nt: IWA max [lambda/D]"),
 ("iwa_lo", 2.0, 1.5, 3.0, "S_obs nt: IWA min [lambda/D]"),
 ("n_bins", 10.0, 1.0, 30.0, "S_obs: canales del piso"),
 ("xuv_mid", 10.0, 3.0, 30.0, "S_estrella: F_XUV actual con 50% de estres"),
 ("tsat_M", 1.0, 0.5, 3.0, "S_estrella: t_sat M<0.6 [Gyr]"),
 ("lpco2_hi", 2.0, 0.0, 2.0, "Clima M: log10 P_CO2 max"),
 ("olr_Amin", 150.0, 120.0, 190.0, "Clima: A minimo de OLR [W/m2]"),
 ("flare_late", 0.5, 0.2, 0.8, "S_estrella: flares M<0.35"),
 ("cloud_min", 0.5, 0.1, 1.0, "S_obs: nubes en transmision"),
 ("sig_model", 0.5, 0.2, 1.0, "S_obs: incertidumbre de SNR (ln)"),
 ("yang_eff_lo", 0.5, 0.0, 1.0, "Clima 3D: eficacia minima (rotacion rapida)"),
 # --- D1 ---
 ("lock_m_c", 0.60, 0.50, 0.70, "D1: centro de M_lock (U +-0.05) [Msun]"),
 ("lock_a_c", 0.20, 0.12, 0.30, "D1: centro de a_lock (U +-0.05) [AU]"),
 ("prot_lo_c", 5.0, 3.0, 7.0, "D1: P_lo de la eficacia 3D (U +-1) [d]"),
 ("prot_hi_c", 20.0, 12.0, 30.0, "D1: P_hi de la eficacia 3D (U +-5) [d]"),
 # --- D2 ---
 ("log_rsat", -3.05, -3.20, -2.90, "D2: log10 R_sat (Wright+2018 via Pass)"),
 ("tsat_fac", 1.0, 0.5, 2.0, "D2: factor sobre t_sat de Pass+2025"),
 ("fflare_c", 1.25, 1.0, 1.5, "D2: f_flare (fijo en el draw)"),
 ("post_drop_log", -2.0, -3.0, -1.5, "D2: log10 caida post-saturacion (FC)"),
 ("pct_c", 50.0, 25.0, 75.0, "D2: centro del percentil de rotacion Mors (U +-20)"),
 ("shore_logc", 0.0, -1.0, 1.0, "D2: desplazamiento log10 de la costa XUV (ARM + c)"),
 ("shore_sig", 1.0, 0.5, 2.0, "D2: ancho ln de la costa"),
 # --- D4 ---
 ("wv_c", 0.14, 0.05, 0.25, "D4: ancho del valle de radios [R_E] (U +-0.03)"),
 ("mv_c", 4.5, 2.0, 8.0, "D4: M_v de la ruta de masa [M_E] (U +-1.5)"),
 ("log2_piH", 0.0, -1.0, 1.0, "D4: log2 prior relativo de la clase Hycean"),
 ("lph2_c", 1.0, 0.3, 1.7, "D4: centro de log10 P_H2 [bar] (U +-0.3)"),
 ("fcl_c", 1.5, 0.8, 2.5, "D4: tolerancia sobre S_crit de Innes (U +-0.25)"),
 ("eta_log", 0.0, -1.0, 1.0, "D4: log10 factor sobre eta de Owen & Wu"),
 ("tloss_ref_log", 0.0, -1.0, 0.65, "D4: log10 t_loss de referencia [Gyr]"),
 ("hyc_alb", 0.12, 0.05, 0.60, "D4: albedo de Bond Hycean (Innes 0.12; Leconte+2024 >0.5)"),
]
NAMES = [p[0] for p in PARAMS]; BOUNDS = [[p[2], p[3]] for p in PARAMS]
KEY = ["TOI-7166 b", "GJ 3378 b", "GJ 273 b", "Proxima Cen b", "Proxima Cen d", "TRAPPIST-1 d", "TRAPPIST-1 e",
       "TRAPPIST-1 f", "TOI-700 d", "TOI-700 e", "K2-18 b", "LHS 1140 b", "Kepler-186 f", "TOI-4336 A b"]
BR = ["rocky", "hycean", "combined"]
ROOT = os.path.dirname(HERE)
_D = None; _THR = None


def to_P(x):
    P = {**u.V3, **v6.V6}; P["snr_thr"] = _THR; cfg = dict(v6.CFG_V6)
    s = dict(zip(NAMES, x))
    for k in ["floor_ppm", "n_H", "iwa_hi", "iwa_lo", "n_bins", "xuv_mid", "tsat_M", "lpco2_hi", "olr_Amin", "flare_late",
              "cloud_min", "sig_model", "yang_eff_lo", "log_rsat", "tsat_fac", "shore_logc", "shore_sig", "hyc_alb"]:
        P[k] = float(s[k])
    P["iwa_hi"] = max(P["iwa_hi"], P["iwa_lo"])
    P["lock_m_lo"], P["lock_m_hi"] = s["lock_m_c"] - 0.05, s["lock_m_c"] + 0.05
    P["lock_a_lo"], P["lock_a_hi"] = s["lock_a_c"] - 0.05, s["lock_a_c"] + 0.05
    P["prot_lo_lo"], P["prot_lo_hi"] = s["prot_lo_c"] - 1, s["prot_lo_c"] + 1
    P["prot_hi_lo"], P["prot_hi_hi"] = s["prot_hi_c"] - 5, s["prot_hi_c"] + 5
    P["fflare_lo"] = P["fflare_hi"] = s["fflare_c"]
    P["post_drop"] = 10 ** s["post_drop_log"]
    P["pct_lo"], P["pct_hi"] = max(5.0, s["pct_c"] - 20), min(95.0, s["pct_c"] + 20)
    P["wv_lo"], P["wv_hi"] = max(0.02, s["wv_c"] - 0.03), s["wv_c"] + 0.03
    P["mv_lo"], P["mv_hi"] = s["mv_c"] - 1.5, s["mv_c"] + 1.5
    P["lph2_lo"], P["lph2_hi"] = s["lph2_c"] - 0.3, s["lph2_c"] + 0.3
    P["fcl_lo"], P["fcl_hi"] = s["fcl_c"] - 0.25, s["fcl_c"] + 0.25
    P["eta_fac"] = 10 ** s["eta_log"]; P["tloss_ref"] = 10 ** s["tloss_ref_log"]
    cfg["class_prior"] = (1.0, 2.0 ** s["log2_piH"], 1.0)
    return P, cfg


def init(thr):
    global _D, _THR
    _THR = thr; u.V3["snr_thr"] = thr
    names = [l.rstrip("\n") for l in open(os.path.join(ROOT, "data", "sens_planets_v5.txt")) if l.strip()]
    df = pd.read_csv(os.path.join(ROOT, "data", "pscomppars_2026-11-18.csv"))
    _D = v6._tiers(u.prepare(df)); _D = _D[_D.pl_name.isin(names)].reset_index(drop=True)


def evaluate(x, n=100):
    P, cfg = to_P(x)
    out = {b: np.zeros(len(_D)) for b in BR}
    for i, row in enumerate(_D.itertuples()):
        rng = np.random.default_rng(zlib.crc32(str(row.pl_name).encode()))
        res = v6.simulate_planet_v6(row, n, rng, cfg=cfg, P=P)
        for b in BR:
            out[b][i] = res[b]["UHI_rank"]
    return np.vstack([out[b] for b in BR])


def work(args):
    j, x = args
    return j, evaluate(x)


def outputs(Y, base):
    """Y: (3, planetas) -> vector de salidas."""
    names = list(_D.pl_name); vals = []
    for k, b in enumerate(BR):
        y, y0 = Y[k], base[k]
        t0 = set(np.argsort(-y0, kind="stable")[:10]); t = set(np.argsort(-y, kind="stable")[:10])
        vals += [len(t0 & t) / 10.0, spearmanr(y0, y).statistic if np.std(y) > 0 else 0.0, float((y > 0).sum())]
        vals += [y[names.index(p)] if p in names else np.nan for p in KEY]
    return np.array(vals)


def out_names():
    return [f"{b}:{o}" for b in BR for o in ["top10_overlap", "spearman_vs_base", "count_alive"] + [f"UHI_{p}" for p in KEY]]


if __name__ == "__main__":
    from SALib.sample import sobol as sobol_sample
    from SALib.analyze import sobol
    thr, N, W, outdir = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    os.makedirs(outdir, exist_ok=True)
    init(thr)
    prob = dict(num_vars=len(NAMES), names=NAMES, bounds=BOUNDS)
    X0 = sobol_sample.sample(prob, N, calc_second_order=False, seed=11)   # Saltelli/Sobol aleatorizado (SALib)
    base_x = np.array([p[1] for p in PARAMS])
    t = time.time(); base = evaluate(base_x); print("base eval", round(time.time() - t, 1), "s", flush=True)
    np.save(os.path.join(outdir, "base_raw.npy"), base)
    raw_path = os.path.join(outdir, "Y_raw.npy")
    R = np.full((len(X0), 3, len(_D)), np.nan, np.float32)
    # reanudacion: usa el checkpoint valido con mas filas completas (Y_raw.npy o copias rotativas Y_raw_bk0/1.npy)
    best = -1
    for cand in [raw_path, os.path.join(outdir, "Y_raw_bk0.npy"), os.path.join(outdir, "Y_raw_bk1.npy")]:
        if not os.path.exists(cand):
            continue
        try:
            R0 = np.load(cand)
        except Exception as e:
            print("checkpoint ilegible", cand, repr(e), flush=True); continue
        if R0.shape == R.shape and int(np.isfinite(R0).all(axis=(1, 2)).sum()) > best:
            best = int(np.isfinite(R0).all(axis=(1, 2)).sum()); R = R0; print("reanuda desde", cand, "filas_completas", best, flush=True)
    todo = [(j, X0[j]) for j in range(len(X0)) if not np.isfinite(R[j]).all()]
    print(f"thr={thr} N={N} D={len(NAMES)} evals={len(X0)} pendientes={len(todo)}", flush=True)
    # v6.1 ejecucion (20 nov 2026): checkpoint atomico frecuente + copia rotativa; no cambia la ciencia ni las muestras
    CK_ROWS, CK_SECS, BK_SECS = int(os.environ.get("UHI_CK_ROWS", 50)), float(os.environ.get("UHI_CK_SECS", 120)), float(os.environ.get("UHI_BK_SECS", 1800))
    def save_atomic(path, arr):
        tmp = path + ".tmp"
        with open(tmp, "wb") as fh:
            np.save(fh, arr); fh.flush(); os.fsync(fh.fileno())
        os.replace(tmp, path)
        dfd = os.open(os.path.dirname(os.path.abspath(path)), os.O_RDONLY); os.fsync(dfd); os.close(dfd)
    def checkpoint(c, force_bk=False):
        nonlocal_state["ck"] = time.time(); save_atomic(raw_path, R)
        done = int(np.isfinite(R).all(axis=(1, 2)).sum())
        print(c, len(todo), f"{time.time() - t:.0f}s", f"filas_completas={done}/{len(X0)}", time.strftime("%H:%M:%S"), flush=True)
        if force_bk or time.time() - nonlocal_state["bk"] >= BK_SECS:
            nonlocal_state["bk"] = time.time()
            save_atomic(os.path.join(outdir, "Y_raw_bk%d.npy" % nonlocal_state["k"]), R); nonlocal_state["k"] ^= 1
    nonlocal_state = dict(ck=time.time(), bk=time.time(), k=0)
    with Pool(W, initializer=init, initargs=(thr,), maxtasksperchild=200) as pool:
        for c, (j, y) in enumerate(pool.imap_unordered(work, todo, chunksize=1)):
            R[j] = y
            if (c + 1) % CK_ROWS == 0 or time.time() - nonlocal_state["ck"] >= CK_SECS:
                checkpoint(c + 1)
    checkpoint(len(todo), force_bk=True)
    Y = np.vstack([outputs(R[j], base) for j in range(len(X0))]); np.save(os.path.join(outdir, "Y_outputs.npy"), Y)
    rows = []
    for k, on in enumerate(out_names()):
        y = Y[:, k]
        if not np.isfinite(y).all() or np.std(y) == 0:
            rows += [dict(output=on, param=p, S1=np.nan, S1_conf=np.nan, ST=np.nan, ST_conf=np.nan, constant=True) for p in NAMES]
            continue
        Si = sobol.analyze(prob, y, calc_second_order=False, seed=11)
        rows += [dict(output=on, param=p, S1=Si["S1"][i], S1_conf=Si["S1_conf"][i], ST=Si["ST"][i], ST_conf=Si["ST_conf"][i], constant=False)
                 for i, p in enumerate(NAMES)]
    pd.DataFrame(rows).to_csv(os.path.join(outdir, "sobol_results.csv"), index=False)
    json.dump(dict(thr=thr, N=N, D=len(NAMES), evals=len(X0), draws=100, planets=len(_D), seconds=round(time.time() - t)),
              open(os.path.join(outdir, "meta.json"), "w"))
    print("listo", round(time.time() - t), "s", flush=True)
