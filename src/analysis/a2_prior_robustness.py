r"""LEGACY AUDIT SCRIPT (v6.0-to-v6.1 review). Kept as a record; it is not rerunnable from the release
as shipped because its input paths or inputs are not included (see README, "Legacy scripts").

Independent class-prior robustness check for rocky, Hycean and combined branches at 5 and 3 sigma.

Run from the repository root. It reads results/uhi_v6_<branch>_<variant>_<thr>s.csv (the v6.1 release
stores them in results/v6.1) for main and ten prior-test variants. It reports top-10/top-20 overlap,
Spearman correlation among planets alive in either run, Kendall tau on the top-20 union, survivor
counts, top-ranked planets, entries/exits and the maximum UHI_rank change. A scenario is flagged as
material if top-10 overlap is <=6 or Spearman <0.8 (reporting thresholds). Writes
results/analysis_phase1/prior_robustness_all_branches.csv."""
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
            # P_top10 sums to ~10 when >=10 planets survive per draw; weighted overlap
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
