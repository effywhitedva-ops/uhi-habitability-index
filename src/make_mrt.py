"""Escribe la tabla completa de resultados (N = 6112) en formato AAS Machine-Readable Table (MRT) con astropy.
Uso: python make_mrt.py uhi_results_v52.csv uhi_results_v52.mrt"""
import sys, numpy as np, pandas as pd
from astropy.table import Table
COLS = {"rank_pos": ("", "UHI v5.2 rank position (5 sigma, by UHI_rank)"), "planet": ("", "Planet name (NASA Exoplanet Archive pl_name)"),
        "UHI_rank": ("", "UHI_rank = P_alive x conditional median UHI"), "P_alive": ("", "Probability UHI > 0 (hurdle)"),
        "UHI_cond_med": ("", "Median UHI conditional on UHI > 0"), "UHI_cond_q25": ("", "25th percentile of UHI | UHI > 0"),
        "UHI_cond_q75": ("", "75th percentile of UHI | UHI > 0"), "P_top10": ("", "Probability of being in the joint top 10"),
        "rank_med": ("", "Median joint rank across draws"), "rank_ci5": ("", "5th percentile of joint rank"), "rank_ci95": ("", "95th percentile of joint rank"),
        "UHI_med": ("", "Unconditional median UHI"), "UHI_mean": ("", "Posterior mean UHI"), "tier": ("", "Observability tier T0-T3"),
        "ESI_med": ("", "Earth Similarity Index (Schulze-Makuch+2011, 4 properties), MC median"),
        "S_fis_med": ("", "Median S_fis"), "S_clima_med": ("", "Median S_clima"), "S_atm_med": ("", "Median S_atm"), "S_estrella_med": ("", "Median S_star"),
        "S_geo_med": ("", "Median S_geo"), "S_mag_med": ("", "Median S_mag"), "S_obs_med": ("", "Median S_obs"),
        "Rp_med": ("", "Median planet radius (Earth radii)"), "Mp_med": ("", "Median planet mass (Earth masses)"), "S_med": ("", "Median instellation (Earth units)"),
        "Teq_med": ("K", "Median equilibrium temperature"), "albedo_med": ("", "Median Bond albedo"), "transiting": ("", "1 = transiting"),
        "st_teff": ("K", "Host effective temperature (catalog)"), "sy_dist": ("pc", "System distance (catalog)")}
r = pd.read_csv(sys.argv[1]).sort_values("rank_pos")[list(COLS)]
t = Table.from_pandas(r)
for c, (u, d) in COLS.items():
    t[c].description = d
    if u: t[c].unit = u
    if t[c].dtype.kind == "f": t[c].format = "%.6g"
t.write(sys.argv[2], format="ascii.mrt", overwrite=True)
print("MRT", sys.argv[2], len(t), "rows")
