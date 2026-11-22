"""Resumen y figura del Sobol v6 (5 sigma y 3 sigma): S_T medio por parametro sobre salidas no constantes, por rama."""
import os, json, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
res = {}
fig, axs = plt.subplots(2, 3, figsize=(17, 11))
for i, t in enumerate((5, 3)):
    d = pd.read_csv(os.path.join(ROOT, f"sens_out_v6_{t}s", "sobol_results.csv")); d = d[~d.constant.astype(bool)]
    d["branch"] = d.output.str.split(":").str[0]
    for j, b in enumerate(["rocky", "hycean", "combined"]):
        s = d[d.branch == b].groupby("param").agg(ST_mean=("ST", "mean"), ST_sd=("ST", "std"), S1_mean=("S1", "mean"), n=("ST", "size")).sort_values("ST_mean", ascending=False)
        s.to_csv(os.path.join(ROOT, "tables", f"sobol_v6_{b}_{t}s_meanST.csv"))
        res[f"{b}_{t}s"] = s.head(10).round(4).reset_index().values.tolist()
        top = s.head(15)[::-1]
        axs[i, j].barh(top.index, top.ST_mean, xerr=top.ST_sd, color="#3182bd", alpha=0.85)
        axs[i, j].barh(top.index, top.S1_mean.clip(lower=0), color="#08306b", alpha=0.9, height=0.4)
        axs[i, j].set_title(f"{b} - {t}σ (S_T medio ± sd sobre {int(s.n.max())} salidas; S1 oscuro)", fontsize=10)
        axs[i, j].tick_params(axis="y", labelsize=8)
fig.suptitle("UHI v6 Fase 1 - Sobol (Saltelli, N=256, D=32, 100 draws x 1214 planetas)"); fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig(os.path.join(ROOT, "figures", "fig_v6_sobol_ST.png"), dpi=140)
json.dump(res, open(os.path.join(ROOT, "tables", "sobol_v6_top.json"), "w"), indent=1)
for k, v in res.items():
    print(k, v[:6])
