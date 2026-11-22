"""Robustez al prior (independiente de prior_sensitivity_v6.py): rocky, hycean, combinado; 5 y 3 sigma."""
import numpy as np, pandas as pd, json
from scipy.stats import spearmanr, kendalltau
V = ["pA1_rocoso","pA2_volatil","pA3_ck17","pA4_empM","pA5_empM_P5","pB1_msini_min","pB2_msini_howard","pC1_mr_ck17","pD1_linea_earthlike","pD2_linea_h2o100"]
L = lambda b, v, t: pd.read_csv(f"results/uhi_v6_{b}_{v}_{t}s.csv").set_index("planet")
rows = []
for t in (5, 3):
    for b in ("rocky", "hycean", "combined"):
        m = L(b, "main", t)
        for v in V:
            d = L(b, v, t)
            t10m, t10d = list(m.index[:10]), list(d.index[:10])
            u20 = list(dict.fromkeys(list(m.index[:20]) + list(d.index[:20])))
            alive = sorted(set(m.index[m.P_alive > 0]) | set(d.index[d.P_alive > 0]))
            # P_top10 suma ~10 si hay >=10 vivos por draw; overlap ponderado
            rows.append(dict(thr=t, branch=b, variant=v, top10_overlap=len(set(t10m) & set(t10d)),
                             top20_overlap=len(set(m.index[:20]) & set(d.index[:20])),
                             spearman_alive_union=round(spearmanr(m.loc[alive, "UHI_rank"], d.reindex(alive).UHI_rank.fillna(0)).statistic, 3),
                             kendall_top20_union=round(kendalltau(m.loc[u20, "UHI_rank"], d.reindex(u20).UHI_rank.fillna(0)).statistic, 3),
                             n_alive_main=int((m.P_alive > 0).sum()), n_alive_var=int((d.P_alive > 0).sum()),
                             n1_main=m.index[0], n1_var=d.index[0],
                             entran=";".join(p for p in t10d if p not in t10m), salen=";".join(p for p in t10m if p not in t10d),
                             max_abs_dUHI_top10main=round(float((d.reindex(t10m).UHI_rank - m.loc[t10m].UHI_rank).abs().max()), 4)))
df = pd.DataFrame(rows)
df["material"] = (df.top10_overlap <= 6) | (df.spearman_alive_union < 0.8)
df.to_csv("results/analysis_phase1/prior_robustness_all_branches.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60)
print(df.drop(columns=["entran","salen"]).to_string())
print(df[df.material][["thr","branch","variant","top10_overlap","spearman_alive_union","entran","salen"]].to_string())
