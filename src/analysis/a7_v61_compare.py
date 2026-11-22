"""v6.1 vs v6.0: vivos, top 10, robustez al prior (3 ramas), Innes vs Madhusudhan, regresion legacy. Escribe results/analysis_phase1/v61_*."""
import numpy as np, pandas as pd, json, glob, os
from scipy.stats import spearmanr, kendalltau
O = "results/analysis_phase1"; out = {}
V = ["pA1_rocoso","pA2_volatil","pA3_ck17","pA4_empM","pA5_empM_P5","pB1_msini_min","pB2_msini_howard","pC1_mr_ck17","pD1_linea_earthlike","pD2_linea_h2o100"]
L = lambda d, b, v, t: pd.read_csv(f"{d}/uhi_v6_{b}_{v}_{t}s.csv").set_index("planet")
# 0) regresion legacy (fixes apagados) vs v6.0
reg = {}
for f in sorted(glob.glob("results/v6.1_legacy_check/uhi_v6_*_main_5s.csv")):
    a = pd.read_csv("results/" + os.path.basename(f)); b = pd.read_csv(f)
    cols = ["UHI_rank", "P_alive", "UHI_cond_med", "P_top10", "rank_med", "P_R", "P_H", "P_G"]
    reg[os.path.basename(f)] = dict(order_equal=bool((a.planet.values == b.planet.values).all()), max_abs_diff=max(float(np.nanmax(np.abs(a[c] - b[c]))) for c in cols))
out["legacy_regression"] = reg
# 1) vivos y top 10
alive, tops = {}, {}
for t in (5, 3):
    for b in ("rocky", "hycean", "combined"):
        for v in ("main", "madhu"):
            o, n = L("results", b, v, t), L("results/v6.1", b, v, t)
            alive[f"{b}_{v}_{t}s"] = dict(v60_gt0=int((o.P_alive > 0).sum()), v61_gt0=int((n.P_alive > 0).sum()), v60_ge001=int((o.P_alive >= .01).sum()), v61_ge001=int((n.P_alive >= .01).sum()))
            if v == "main":
                to, tn = list(o.index[:10]), list(n.index[:10]); al = sorted(set(o.index[o.P_alive >= .01]) | set(n.index[n.P_alive >= .01]))
                tops[f"{b}_{t}s"] = dict(overlap=len(set(to) & set(tn)), entran=[p for p in tn if p not in to], salen=[p for p in to if p not in tn],
                    spearman_alive001=round(float(spearmanr(o.loc[al, "UHI_rank"], n.reindex(al).UHI_rank.fillna(0)).statistic), 3),
                    top10_v61=[(p, round(float(n.UHI_rank[p]), 4), int(o.rank_pos[p])) for p in tn], max_abs_dUHI=round(float((n.UHI_rank - o.UHI_rank.reindex(n.index)).abs().max()), 4))
out["alive"] = alive; out["top10_shift"] = tops
# 2) robustez al prior v6.1 vs v6.0 (mismo calculo que a2)
rows = []
for d, lab in (("results", "v6.0"), ("results/v6.1", "v6.1")):
    for t in (5, 3):
        for b in ("rocky", "hycean", "combined"):
            m = L(d, b, "main", t)
            for v in V:
                x = L(d, b, v, t); al = sorted(set(m.index[m.P_alive > 0]) | set(x.index[x.P_alive > 0]))
                rows.append(dict(version=lab, thr=t, branch=b, variant=v, top10_overlap=len(set(m.index[:10]) & set(x.index[:10])),
                                 spearman=round(float(spearmanr(m.loc[al, "UHI_rank"], x.reindex(al).UHI_rank.fillna(0)).statistic), 3),
                                 salen=";".join(p for p in m.index[:10] if p not in set(x.index[:10]))))
pr = pd.DataFrame(rows); pr["material"] = (pr.top10_overlap <= 6) | (pr.spearman < 0.8); pr.to_csv(f"{O}/v61_prior_robustness.csv", index=False)
piv = pr.pivot_table(index=["thr", "branch", "variant"], columns="version", values="top10_overlap").reset_index()
out["prior_overlap_table"] = piv.values.tolist()
# 3) Innes vs Madhu v6.1
im = {}
for t in (5, 3):
    I, M = L("results/v6.1", "hycean", "main", t), L("results/v6.1", "hycean", "madhu", t)
    M = M.reindex(I.index); aI, aM = I.P_alive > 0, M.P_alive > 0; both = aI & aM
    tI = list(I.sort_values("UHI_rank", ascending=False).index[:10]); tM = list(M.sort_values("UHI_rank", ascending=False).index[:10])
    Cb = L("results/v6.1", "combined", "main", t); Cm = L("results/v6.1", "combined", "madhu", t)
    im[f"{t}s"] = dict(alive_innes=int(aI.sum()), alive_madhu=int(aM.sum()), alive_innes_001=int((I.P_alive >= .01).sum()), alive_madhu_001=int((M.P_alive >= .01).sum()),
        top10_overlap=len(set(tI) & set(tM)), spearman_both=round(float(spearmanr(I.loc[both, "UHI_rank"], M.loc[both, "UHI_rank"]).statistic), 3),
        UHI_max_innes=round(float(I.UHI_rank.max()), 4), UHI_max_madhu=round(float(M.UHI_rank.max()), 4), top10_innes=tI, top10_madhu=tM,
        combined_top10_overlap=len(set(Cb.index[:10]) & set(Cm.index[:10])))
out["innes_vs_madhu_v61"] = im
json.dump(out, open(f"{O}/v61_vs_v60.json", "w"), indent=1, default=str)
print(json.dumps({k: out[k] for k in ("legacy_regression", "alive", "innes_vs_madhu_v61")}, indent=0, default=str))
for k, v in tops.items(): print(k, {kk: vv for kk, vv in v.items() if kk != "top10_v61"}); print("   ", v["top10_v61"])
print(pr[pr.version == "v6.1"].to_string())
