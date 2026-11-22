import pandas as pd, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, axs = plt.subplots(2, 3, figsize=(18, 11))
for i, t in enumerate((5, 3)):
    for j, b in enumerate(("rocky", "hycean", "combined")):
        d = pd.read_csv(f"results/analysis_phase1/sobol_robust_{b}_{t}s.csv").head(15)[::-1]
        ax = axs[i, j]; col = np.where(d.ci_width < 0.05, "#31a354", "#de2d26")
        ax.barh(d.param, d.ST_mean, xerr=[d.ST_mean - d.ci_lo, d.ci_hi - d.ST_mean], color=col, alpha=.85, capsize=2)
        ax.barh(d.param, d.S1_mean.clip(lower=0), color="#08306b", height=.35)
        ax.set_title(f"{b} {t}σ: S_T medio robusto (IC95 bootstrap); verde = ancho IC < 0.05", fontsize=9); ax.tick_params(axis="y", labelsize=8)
fig.suptitle("UHI v6 Fase 1: Sobol corregido (sin salidas degeneradas con UHI casi siempre 0); S1 en azul oscuro")
fig.tight_layout(rect=(0, 0, 1, .97)); fig.savefig("figures/analysis_phase1/fig_a1_sobol_ST_robusto.png", dpi=130)
