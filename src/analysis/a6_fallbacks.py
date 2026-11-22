r"""LEGACY AUDIT SCRIPT (v6.0-to-v6.1 review). Kept as a record; it is not rerunnable from the release
as shipped because its input paths or inputs are not included (see README, "Legacy scripts").

Numerical checks of documented fallbacks (docs/assumptions_v6_es.md, fallback list).

(1) Gladman et al. (1996) tidal-locking time for an Earth-like planet as a function of stellar mass and
orbital distance, compared with the D1 p_sync transition; Q=100, k2=0.3, initial 24 h rotation and
alpha=0.33 are illustrative values, with Q=10 and Q=1000 as brackets. (3) Erkaev et al. (2007)
Roche-lobe enhancement K for surviving Hycean planets, to estimate how much ignoring it shifts S_atm.

Run from the repository root. Reads results/uhi_v6_hycean_{main,madhu}_5s.csv (v6.0-era location;
results/v6.1 in the release) and writes results/analysis_phase1/fallback_checks.json."""
import numpy as np, pandas as pd, json
G=6.674e-11; MSUN=1.98847e30; AU=1.495978707e11; ME=5.9722e24; RE=6.371e6; GYR=3.15576e16
def tlock(a_au, mstar, mp=1.0, rp=1.0, Q=100., k2=0.3, P0_h=24., alpha=0.33):
    r"""Return tidal-locking time in Gyr (Gladman et al. 1996 form).

    a_au AU, mstar M_sun, mp M_E, rp R_E, Q and k2 dimensionless, P0_h initial spin period in hours,
    alpha the moment-of-inertia factor. Uses SI internally."""
    w=2*np.pi/(P0_h*3600); I=alpha*mp*ME*(rp*RE)**2
    return w*(a_au*AU)**6*I*Q/(3*G*(mstar*MSUN)**2*k2*(rp*RE)**5)/GYR
out={"tlock_Gyr_Earthlike_Q100_k2_0.3_P0_24h":{}}
for m in (0.1,0.3,0.5,0.6,0.8,1.0):
    for a in (0.05,0.1,0.2,0.3,0.4):
        out["tlock_Gyr_Earthlike_Q100_k2_0.3_P0_24h"][f"M{m}_a{a}"]=float(f"{tlock(a,m):.3g}")
# semimajor axis where t_lock = 1 Gyr for each M* (Q=100, plus Q=10 and Q=1000)
out["a_lock_1Gyr_AU"]={m: round(float((1/tlock(1,m))**(1/6)),3) for m in (0.1,0.3,0.5,0.6,0.8,1.0)}
out["a_lock_1Gyr_AU_Q10"]={m: round(float((1/tlock(1,m,Q=10))**(1/6)),3) for m in (0.1,0.3,0.5,0.6,0.8,1.0)}
out["a_lock_1Gyr_AU_Q1000"]={m: round(float((1/tlock(1,m,Q=1000))**(1/6)),3) for m in (0.1,0.3,0.5,0.6,0.8,1.0)}
# Erkaev K for surviving Hycean planets (Innes and Madhusudhan), using class medians
for v in ("main","madhu"):
    H=pd.read_csv(f"results/uhi_v6_hycean_{v}_5s.csv"); H=H[H.P_alive>0]
    xi=(H.Mp_cls_med*ME/(3*H.st_mass*MSUN))**(1/3)*H.a_med*AU/(H.Rp_cls_med*RE)
    K=1-1.5/xi+0.5/xi**3
    out[f"erkaev_K_{v}"]=dict(n=int(len(H)),K_min=round(float(K.min()),3),K_p5=round(float(K.quantile(.05)),3),K_med=round(float(K.median()),3),
        dlnS_atm_shift_med=round(float(np.log(K).median()),3), n_K_lt_0_9=int((K<0.9).sum()),
        top10_K={p:round(float(k),3) for p,k in zip(H.sort_values("UHI_rank",ascending=False).planet[:10],K.loc[H.sort_values("UHI_rank",ascending=False).index[:10]])})
json.dump(out,open("results/analysis_phase1/fallback_checks.json","w"),indent=1); print(json.dumps(out,indent=1))
