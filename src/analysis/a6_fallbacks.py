"""Chequeos de los 5 fallbacks: (1) t_lock de Gladman+1996 en la frontera p_sync; (3) K de Erkaev+2007 en la rama Hycean."""
import numpy as np, pandas as pd, json
G=6.674e-11; MSUN=1.98847e30; AU=1.495978707e11; ME=5.9722e24; RE=6.371e6; GYR=3.15576e16
def tlock(a_au, mstar, mp=1.0, rp=1.0, Q=100., k2=0.3, P0_h=24., alpha=0.33):
    w=2*np.pi/(P0_h*3600); I=alpha*mp*ME*(rp*RE)**2
    return w*(a_au*AU)**6*I*Q/(3*G*(mstar*MSUN)**2*k2*(rp*RE)**5)/GYR
out={"tlock_Gyr_Earthlike_Q100_k2_0.3_P0_24h":{}}
for m in (0.1,0.3,0.5,0.6,0.8,1.0):
    for a in (0.05,0.1,0.2,0.3,0.4):
        out["tlock_Gyr_Earthlike_Q100_k2_0.3_P0_24h"][f"M{m}_a{a}"]=float(f"{tlock(a,m):.3g}")
# a donde t_lock = 1 Gyr y 10 Gyr para cada M*
out["a_lock_1Gyr_AU"]={m: round(float((1/tlock(1,m))**(1/6)),3) for m in (0.1,0.3,0.5,0.6,0.8,1.0)}
out["a_lock_1Gyr_AU_Q10"]={m: round(float((1/tlock(1,m,Q=10))**(1/6)),3) for m in (0.1,0.3,0.5,0.6,0.8,1.0)}
out["a_lock_1Gyr_AU_Q1000"]={m: round(float((1/tlock(1,m,Q=1000))**(1/6)),3) for m in (0.1,0.3,0.5,0.6,0.8,1.0)}
# Erkaev K en planetas Hycean vivos (Innes y Madhu), con medianas de clase
for v in ("main","madhu"):
    H=pd.read_csv(f"results/uhi_v6_hycean_{v}_5s.csv"); H=H[H.P_alive>0]
    xi=(H.Mp_cls_med*ME/(3*H.st_mass*MSUN))**(1/3)*H.a_med*AU/(H.Rp_cls_med*RE)
    K=1-1.5/xi+0.5/xi**3
    out[f"erkaev_K_{v}"]=dict(n=int(len(H)),K_min=round(float(K.min()),3),K_p5=round(float(K.quantile(.05)),3),K_med=round(float(K.median()),3),
        dlnS_atm_shift_med=round(float(np.log(K).median()),3), n_K_lt_0_9=int((K<0.9).sum()),
        top10_K={p:round(float(k),3) for p,k in zip(H.sort_values("UHI_rank",ascending=False).planet[:10],K.loc[H.sort_values("UHI_rank",ascending=False).index[:10]])})
json.dump(out,open("results/analysis_phase1/fallback_checks.json","w"),indent=1); print(json.dumps(out,indent=1))
