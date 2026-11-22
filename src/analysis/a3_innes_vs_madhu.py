"""Comparacion formal Innes+2023 (principal) vs Madhusudhan+2021 (variante) en la rama Hycean. Independiente del pipeline."""
import numpy as np, pandas as pd, json, sys
sys.path.insert(0, "src"); import xuv_history as X
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import spearmanr
OUT = "results/analysis_phase1"; FIG = "figures/analysis_phase1"; import os; os.makedirs(FIG, exist_ok=True)
summ = {}
sub = ["S_fis_med", "S_atm_med", "S_estrella_med", "S_geo_med", "S_mag_med", "S_obs_med", "P_H", "tloss_all_med"]
fig, axs = plt.subplots(1, 3, figsize=(19, 6.2))
for t in (5, 3):
    I = pd.read_csv(f"results/uhi_v6_hycean_main_{t}s.csv").set_index("planet")
    M = pd.read_csv(f"results/uhi_v6_hycean_madhu_{t}s.csv").set_index("planet").reindex(I.index)
    # 1) que cambia: solo S_clima_H? (todas las demas medianas iguales)
    diffs = {c: float(np.nanmax(np.abs(I[c] - M[c]))) for c in sub}
    aI, aM = I.P_alive > 0, M.P_alive > 0
    # 2) criterio climatico: S/S_crit Innes en el borde mas permisivo de la principal (P=1 bar, f_cl=2)
    teff = I.st_teff.fillna(3457.0); smax = 2.0 * X.s_crit_innes(teff, 1.0)
    teq0 = 278.6 * I.S_med ** 0.25            # T_eq(A=0) para S en S_E (T_eq,Tierra(A=0) = 278.6 K)
    madhu_ok = teq0 <= 500
    innes_never = I.S_med > smax
    both = aI & aM
    tI, tM = list(I.sort_values("UHI_rank", ascending=False).index[:10]), list(M.sort_values("UHI_rank", ascending=False).index[:10])
    summ[f"{t}s"] = dict(max_abs_diff_other_subindices=diffs, alive_innes=int(aI.sum()), alive_madhu=int(aM.sum()), alive_both=int(both.sum()),
        alive_innes_only=int((aI & ~aM).sum()), alive_madhu_only=int((aM & ~aI).sum()),
        madhu_alive_with_S_above_innes_max=int((aM & innes_never).sum()),
        P_clima_H_innes_mean_on_madhu_alive=round(float(I.loc[aM, "P_clima_H"].mean()), 4),
        S_med_quantiles_madhu_alive=[round(float(x), 3) for x in I.loc[aM, "S_med"].quantile([.1, .5, .9])],
        S_med_quantiles_innes_alive=[round(float(x), 3) for x in I.loc[aI, "S_med"].quantile([.1, .5, .9])],
        spearman_on_both_alive=round(float(spearmanr(I.loc[both, "UHI_rank"], M.loc[both, "UHI_rank"]).statistic), 3),
        top10_innes=tI, top10_madhu=tM, overlap=sorted(set(tI) & set(tM)),
        innes_top10_rank_in_madhu={p: int(M.UHI_rank.rank(ascending=False, method="min")[p]) for p in tI},
        madhu_top10_rank_in_innes={p: int(I.UHI_rank.rank(ascending=False, method="min")[p]) for p in tM},
        madhu_top10_S_med={p: round(float(I.S_med[p]), 3) for p in tM}, innes_top10_S_med={p: round(float(I.S_med[p]), 3) for p in tI},
        ratio_UHI_cond_on_both=round(float((M.loc[both, "UHI_cond_med"] / I.loc[both, "UHI_cond_med"]).median()), 3))
    if t == 5:
        ax = axs[0]
        ax.scatter(I.loc[aM & ~aI, "st_teff"], I.loc[aM & ~aI, "S_med"], s=6, c="#fdae6b", label=f"vivo solo con Madhusudhan ({int((aM & ~aI).sum())})")
        ax.scatter(I.loc[both, "st_teff"], I.loc[both, "S_med"], s=8, c="#3182bd", label=f"vivo con ambos ({int(both.sum())})")
        tt = np.linspace(2300, 7000, 200)
        for p_, ls in ((1, "-"), (10, "--"), (100, ":")):
            ax.plot(tt, X.s_crit_innes(tt, p_), "k" + ls, lw=1, label=f"Innes S_crit {p_} bar")
        ax.plot(tt, 2 * X.s_crit_innes(tt, 1), color="gray", lw=1, label="Innes borde más permisivo (1 bar, f_cl=2)")
        ax.axhline((500 / 278.6) ** 4, color="#e6550d", lw=1.5, label="Madhusudhan T_eq(A=0)=500 K")
        for p in ["K2-18 b", "TOI-270 d", "LHS 1140 b", "TOI-7166 b", "Kepler-1630 b", "Wolf 1061 d", "TOI-2257 b"]:
            if p in I.index: ax.annotate(p, (I.st_teff[p], I.S_med[p]), fontsize=7)
        ax.set_yscale("log"); ax.set_xlabel("T_eff [K]"); ax.set_ylabel("S [S_E] (mediana)"); ax.set_ylim(1e-3, 50)
        ax.set_title("Criterio climático Hycean: Innes vs Madhusudhan (5σ)", fontsize=10); ax.legend(fontsize=7, loc="lower right")
        ax = axs[1]
        r = lambda d: d.UHI_rank.rank(ascending=False, method="min")
        x, y = r(I)[aI | aM], r(M)[aI | aM]
        ax.scatter(x, y, s=6, c=np.log10(I.S_med[aI | aM]), cmap="viridis"); ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("puesto Innes (principal)"); ax.set_ylabel("puesto Madhusudhan"); ax.plot([1, 5e3], [1, 5e3], "k:", lw=0.8)
        ax.set_title(f"Puestos Hycean 5σ (color = log S); Spearman vivos-ambos {summ['5s']['spearman_on_both_alive']}", fontsize=10)
        ax = axs[2]
        ax.bar(["Innes", "Madhusudhan"], [aI.sum(), aM.sum()], color=["#3182bd", "#e6550d"])
        ax.set_title("Vivos en la rama Hycean (5σ y 3σ iguales)", fontsize=10)
        for i, v in enumerate([aI.sum(), aM.sum()]): ax.text(i, v, str(int(v)), ha="center", va="bottom")
fig.tight_layout(); fig.savefig(f"{FIG}/fig_a3_innes_vs_madhu.png", dpi=130)
json.dump(summ, open(f"{OUT}/innes_vs_madhu_formal.json", "w"), indent=1)
print(json.dumps(summ, indent=1)[:6000])
