"""Controles sinteticos v5: Tierra, Venus y Marte alrededor de un analogo solar a 10 pc, en geometria
transitante y no transitante. Magnitudes del Sol a 10 pc: M_J = 3.67, M_Ks = 3.27 (Willmer 2018, ApJS 236, 47).
Parametros planetarios: valores estandar del Sistema Solar; errores pequenos (SUPUESTO: 1% en R y M)."""
import numpy as np, pandas as pd, uhi_code as u

PL = {"Earth": (365.256, 1.0, 1.0, 0.0167), "Venus": (224.701, 0.949, 0.815, 0.0068), "Mars": (686.980, 0.532, 0.107, 0.0934)}


def synth_rows():
    rows = []
    for geo, tf in [("tr", 1), ("nt", 0)]:
        for nm, (per, rp, mp, e) in PL.items():
            rows.append(dict(pl_name=f"SYN {nm} ({geo})", hostname="SYN Sun", discoverymethod="synthetic", tran_flag=tf,
                             rv_flag=1, pl_controv_flag=0, pl_orbper=per, pl_orbpererr1=1e-3, pl_orbpererr2=-1e-3,
                             pl_orbsmax=np.nan, pl_orbsmaxerr1=np.nan, pl_orbsmaxerr2=np.nan,
                             pl_rade=rp, pl_radeerr1=0.01 * rp, pl_radeerr2=-0.01 * rp,
                             pl_bmasse=mp, pl_bmasseerr1=0.01 * mp, pl_bmasseerr2=-0.01 * mp, pl_bmassprov="Mass",
                             pl_orbeccen=e, pl_orbeccenerr1=1e-4, pl_orbeccenerr2=-1e-4, pl_insol=np.nan, pl_eqt=np.nan,
                             st_teff=5772.0, st_tefferr1=10.0, st_tefferr2=-10.0, st_rad=1.0, st_raderr1=0.01, st_raderr2=-0.01,
                             st_mass=1.0, st_masserr1=0.01, st_masserr2=-0.01, st_lum=0.0, st_lumerr1=np.nan, st_lumerr2=np.nan,
                             st_age=4.6, st_ageerr1=0.1, st_ageerr2=-0.1, st_met=0.0, sy_dist=10.0, sy_jmag=3.67, sy_kmag=3.27,
                             pl_rade_reflink="", pl_bmasse_reflink=""))
    return pd.DataFrame(rows)


def run_anchors(catalog_csv, n=10000, seed=42):
    res = u.run(synth_rows(), n_draws=n, seed=seed, joint=False)
    cat = pd.read_csv(catalog_csv)
    sc = cat.UHI_rank.values
    res["pct_vs_catalog"] = [100.0 * np.mean(sc < x) for x in res.UHI_rank]
    res["pos_if_inserted"] = [int(np.sum(sc > x)) + 1 for x in res.UHI_rank]
    return res


if __name__ == "__main__":
    import sys
    r = run_anchors(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 10000)
    r.to_csv(sys.argv[3] if len(sys.argv) > 3 else "anchors_v5.csv", index=False)
    print(r[["planet", "UHI_rank", "P_alive", "UHI_cond_med", "pct_vs_catalog", "pos_if_inserted"] + u.SUB_MED if hasattr(u, "SUB_MED") else ["planet", "UHI_rank", "P_alive", "UHI_cond_med", "pct_vs_catalog", "pos_if_inserted"]].round(4).to_string())
