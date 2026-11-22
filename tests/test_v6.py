"""Tests de UHI v6 Fase 1.  python3 tests/test_v6.py [n_compat]
 1) compat: v6 con todos los modulos apagados reproduce v5.2 (catalogo completo, 5 sigma y 3 sigma, y anclas Tierra/Venus/Marte);
 2) ARM manual (Pass+2025 Ec. 10-11 con valores de catalogo) para 13 planetas, vs los valores calculados a mano en el diseno;
 3) planetas ancla: corren en v6 y se reporta clase, ARM, P_ret, UHI por rama."""
import os, sys, time, hashlib, warnings
warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
import numpy as np, pandas as pd
import uhi_code as u, uhi_code_v6 as v6, xuv_history as X, anchors_v5 as A5

DF = pd.read_csv(os.path.join(ROOT, "data", "pscomppars_2026-11-18.csv"))
ANCHORS = ["TRAPPIST-1 d", "TRAPPIST-1 e", "Proxima Cen b", "Proxima Cen d", "Ross 128 b", "K2-125 b", "TOI-4336 A b",
           "TOI-7166 b", "LHS 1140 b", "K2-18 b", "TOI-700 d", "TOI-700 e"]
MANUAL_ARM = {"TRAPPIST-1 d": -1.96, "TRAPPIST-1 e": -1.36, "TRAPPIST-1 f": -0.88, "Proxima Cen b": -0.71, "Proxima Cen d": -2.02,
              "Ross 128 b": -0.77, "GJ 1061 d": -0.47, "GJ 273 b": -0.01, "Wolf 1061 c": -0.01, "GJ 3378 b": -0.10,
              "LHS 1140 b": 0.51, "TOI-7166 b": -0.16, "TOI-4336 A b": -0.30, "Teegarden's Star b": -0.90}
PASS_TABLE2 = {"TRAPPIST-1 d": -1.86, "TRAPPIST-1 e": -1.29}
COLS = ["UHI_rank", "P_alive", "UHI_cond_med", "UHI_mean", "P_top10", "rank_med"] + [s + "_med" for s in u.SUB]


def compat(thr, ref_csv, ref_anchor):
    u.V3["snr_thr"] = thr; P = {**v6.V3, **v6.V6}
    t = time.time(); r6 = v6.run_v6(DF, 10000, seed=42, cfg=v6.CFG_COMPAT, P=P)["combined"]
    ref = pd.read_csv(ref_csv)
    m = ref.merge(r6, on="planet", suffixes=("_5", "_6"))
    assert len(m) == len(ref) == len(r6), (len(m), len(ref), len(r6))
    assert (ref.planet.values == r6.planet.values).all(), "orden distinto"
    mx = max(float(np.nanmax(np.abs(m[c + "_5"] - m[c + "_6"]))) for c in COLS)
    an = A5.synth_rows()
    a6 = v6.run_v6(an, 10000, seed=42, cfg=v6.CFG_COMPAT, P=P, joint=False)["combined"]
    a5 = pd.read_csv(ref_anchor)
    ma = a5.merge(a6, on="planet", suffixes=("_5", "_6"))
    mxa = max(float(np.nanmax(np.abs(ma[c + "_5"] - ma[c + "_6"]))) for c in ["UHI_rank", "P_alive", "UHI_cond_med"])
    print(f"compat {thr:g}sigma: {len(r6)} planetas, max|dif| catalogo = {mx:.3g}, anclas = {mxa:.3g}, vivos {(r6.P_alive>0).sum()}  ({time.time()-t:.0f} s)")
    assert mx < 1e-12 and mxa < 1e-12   # diferencias solo por el redondeo de la escritura CSV de v5.2 (repr float)
    u.V3["snr_thr"] = 5.0


def compat_bitwise(nrows=800):
    """v5.2 en memoria (uhi_code.run) vs v6 con CFG_COMPAT sobre las primeras filas del snapshot: igualdad EXACTA."""
    sub = DF.head(nrows)
    a = u.run(sub, 2000, seed=42); b = v6.run_v6(sub, 2000, seed=42, cfg=v6.CFG_COMPAT, P={**u.V3, **v6.V6})["combined"]
    for c in COLS + ["UHI_p16", "UHI_p84", "ESI_med"]:
        assert np.array_equal(a[c].values, b[c].values, equal_nan=True), c
    assert (a.planet.values == b.planet.values).all()
    print(f"compat bit a bit (v5.2 en memoria vs v6 compat, {len(a)} planetas, 2000 draws): IDENTICO")


def manual_arm():
    v52 = pd.read_csv(os.path.join(ROOT, "data", "uhi_results_v52.csv")).set_index("planet")
    d = DF.set_index("pl_name"); out = []
    for nm, ref in MANUAL_ARM.items():
        r = d.loc[nm]; Ms = r.st_mass; L = 10 ** r.st_lum
        tsat = 7.6 if r.hostname == "TRAPPIST-1" else max(5.9 - 15.4 * Ms, 0.01)
        tt = min(tsat, r.st_age) if np.isfinite(r.st_age) else tsat
        I = 1.25 * 10 ** -3.05 * X.fl_pass(Ms) * L * X.LSUN_CGS * tt * X.GYR_S / (4 * np.pi * (r.pl_orbsmax * X.AU_CM) ** 2) / X.I_XUV_EARTH
        vesc = 11.186 * np.sqrt(v52.loc[nm, "Mp_med"] / v52.loc[nm, "Rp_med"])
        a = X.arm(vesc, I); out.append((nm, ref, round(a, 2)))
    print("ARM (Pass+2025 Ec.10-11, valores de catalogo): planeta, manual-diseno, test")
    for o in out:
        print("  ", o)
    dev = max(abs(o[1] - o[2]) for o in out)
    print("  max |dif| =", round(dev, 3), "| Pass Tabla 2: TRAPPIST-1 d -1.86, e -1.29")
    assert dev <= 0.06, dev


def anchors_v6(n=10000):
    res = v6.run_v6(DF, n, seed=42, cfg=v6.CFG_V6, names=ANCHORS, joint=False)
    R = res["rocky"].set_index("planet"); H = res["hycean"].set_index("planet"); Cb = res["combined"].set_index("planet")
    v52 = pd.read_csv(os.path.join(ROOT, "data", "uhi_results_v52.csv")).set_index("planet")
    tab = pd.DataFrame({"evid": R.evidence, "P_R": R.P_R, "P_H": R.P_H, "P_G": R.P_G, "f_sync": R.f_sync, "ARM_med": R.ARM_med, "ARM_all": R.ARM_all_med,
                        "manual_ARM": pd.Series(MANUAL_ARM), "Pret_R": R.P_retain_med, "UHIrank_R": R.UHI_rank,
                        "UHIrank_H": H.UHI_rank, "UHIrank_C": Cb.UHI_rank, "UHIrank_v52": v52.UHI_rank}).loc[ANCHORS]
    print(tab.round(3).to_string())
    for nm in ANCHORS:
        assert np.isfinite(R.loc[nm, "UHI_rank"]) and np.isfinite(H.loc[nm, "UHI_rank"])
        assert abs(R.loc[nm, "P_R"] + R.loc[nm, "P_H"] + R.loc[nm, "P_G"] - 1) < 1e-9
    for nm in ["TRAPPIST-1 d", "TRAPPIST-1 e", "Proxima Cen b", "Proxima Cen d", "Ross 128 b", "LHS 1140 b", "TOI-4336 A b"]:   # masa medida (TOI-7166 b: masa M-R; v6 usa Otegi, no CK)
        if nm in MANUAL_ARM and np.isfinite(R.loc[nm, "ARM_med"]):
            assert abs(R.loc[nm, "ARM_med"] - MANUAL_ARM[nm]) < 0.35, (nm, R.loc[nm, "ARM_med"])
    return tab


if __name__ == "__main__":
    manual_arm()
    t = anchors_v6(); t.to_csv(os.path.join(ROOT, "tables", "anchors_v6_check.csv"))
    compat_bitwise()
    if "--fast" not in sys.argv:
        compat(5.0, os.path.join(ROOT, "data", "uhi_results_v52.csv"), os.path.join(ROOT, "data", "anchors_v52.csv"))
        compat(3.0, os.path.join(ROOT, "data", "uhi_results_v52_3sigma.csv"), os.path.join(ROOT, "data", "anchors_v52_3sigma.csv"))
    print("TODOS LOS TESTS OK")
