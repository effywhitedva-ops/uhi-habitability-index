"""UHI v5.2 RC: figuras vectoriales (PDF) + vistas previas PNG, generadas desde CSV. Uso: python make_figures.py [dir_base]"""
import os, sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
matplotlib.rcParams.update({"pdf.fonttype": 42, "font.size": 9})
B = sys.argv[1] if len(sys.argv) > 1 else "."
P = lambda *a: os.path.join(B, *a)
OUT = P("figures_v52"); os.makedirs(OUT, exist_ok=True)


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".pdf")); fig.savefig(os.path.join(OUT, name + ".png"), dpi=160); plt.close(fig)


# (a) Sobol post-v5
g = pd.read_csv(P("sens_out_v5", "importance_meanST.csv")).set_index("param").sort_values("mean_ST", ascending=False).head(15).iloc[::-1]
fig, ax = plt.subplots(figsize=(7.2, 6.0)); y = np.arange(len(g))
ax.barh(y + 0.2, g.mean_ST, 0.4, xerr=g.mean_STconf, color="#4c72b0", capsize=2, label="$S_{Ti}$ (total) ± IC 95 %")
ax.barh(y - 0.2, g.mean_S1.clip(lower=0), 0.4, color="#dd8452", label="$S_i$ (primer orden)")
ax.set_yticks(y); ax.set_yticklabels(g.index); ax.set_xlabel("Índice de Sobol, promedio sobre las salidas no constantes")
ax.set_title("Sensibilidad de Sobol post-v5 (UHI_rank): 15 supuestos dominantes de 51\n"
             "Saltelli N=128 (6784 corridas) + 512 LHS, 100 draws/planeta.\nCorrido sobre v5: F_rg (veto v5.1) y rampa w_eff (v5.2) NO incluidos", fontsize=8.5)
ax.legend(loc="lower right", fontsize=8); fig.tight_layout(); save(fig, "fig_a_sobol_post_v5")

# (b) P_alive vs mediana condicional
r = pd.read_csv(P("uhi_results_v52.csv")).sort_values("rank_pos"); a = pd.read_csv(P("anchors_v52.csv"))
lv = r[r.P_alive > 0]; t20 = r.head(20)
fig, ax = plt.subplots(figsize=(10.5, 6.4)); fig.subplots_adjust(right=0.74)
ax.scatter(lv.UHI_cond_med, lv.P_alive, s=6, c="#bbbbbb", label=f"vivos (P_alive>0), N={len(lv)}")
for c in [0.1, 0.2, 0.3, 0.4, 0.5]:
    xx = np.linspace(c, 1, 100); ax.plot(xx, c / xx, c="#dddddd", lw=0.6, zorder=0)
    ax.text(min(1.0, c / 0.985), 0.985, f"{c}", fontsize=6, color="#999999", ha="center")
ax.scatter(t20.UHI_cond_med, t20.P_alive, s=22, c=np.where(t20.transiting == 1, "#1b9e77", "#d95f02"), edgecolor="k", lw=0.4, zorder=3)
for _, q in t20.iterrows():
    ax.annotate(f"{int(q.rank_pos)}", (q.UHI_cond_med, q.P_alive), fontsize=6.5, xytext=(3, 2), textcoords="offset points")
fig.text(0.755, 0.95, "Top 20 (5σ)\n" + "\n".join(f"{int(q.rank_pos)}. {q.planet}" for _, q in t20.iterrows()), fontsize=7, va="top")
mk = {"Earth": ("*", "#1f77b4"), "Venus": ("X", "#d62728"), "Mars": ("P", "#ff7f0e")}
for _, q in a.iterrows():
    nm = q.planet.split()[1]; m, c = mk[nm]
    ax.scatter(q.UHI_cond_med, q.P_alive, marker=m, s=90 if nm == "Earth" else 60, c=c, edgecolor="k", lw=0.5, zorder=4,
               label=f"{q.planet}: P_alive={q.P_alive:.2f}, med.cond.={q.UHI_cond_med:.3f}")
ax.scatter([], [], s=22, c="#1b9e77", edgecolor="k", lw=.4, label="top 20, transita"); ax.scatter([], [], s=22, c="#d95f02", edgecolor="k", lw=.4, label="top 20, no transita")
ax.set_xlabel("mediana condicional UHI | UHI>0"); ax.set_ylabel("P_alive = P(UHI > 0)")
ax.set_xlim(-0.02, 0.75); ax.set_ylim(-0.04, 1.04)
ax.set_title("UHI v5.2 (5σ): P_alive vs mediana condicional (curvas grises: UHI_rank constante)\nVenus y Marte sintéticos a 10 pc en P_alive = 0", fontsize=9)
ax.legend(fontsize=6, loc="center left", bbox_to_anchor=(1.01, 0.2)); save(fig, "fig_b_palive_vs_condmed")

# (c) UHI (5σ, 3σ) vs ESI / PHL
bt = pd.read_csv(P("tables_v52", "uhi_results_v52_with_benchmark.csv"))
r3 = pd.read_csv(P("uhi_results_v52_3sigma.csv")).set_index("planet")
fig, axs = plt.subplots(1, 2, figsize=(12, 6.2), gridspec_kw={"width_ratios": [1, 1.25]})
ax = axs[0]; l5 = bt[bt.P_alive > 0]
ax.scatter(l5.ESI_med, l5.UHI_rank, s=7, c="#4c72b0", alpha=.6, label=f"5σ vivos (N={len(l5)})")
l3 = r3[r3.P_alive > 0]; ax.scatter(l3.ESI_med, l3.UHI_rank, s=7, c="#dd8452", alpha=.45, marker="^", label=f"3σ vivos (N={len(l3)})")
ax.axvline(0.8, c="gray", ls=":", lw=1); ax.text(0.805, ax.get_ylim()[1] * 0.97 if ax.get_ylim()[1] > 0 else 0.5, "ESI = 0.8", fontsize=7, color="gray")
ax.set_xlabel("ESI (Schulze-Makuch+2011, 4 propiedades; mediana MC)"); ax.set_ylabel("UHI_rank")
ax.set_title("UHI_rank vs ESI (planetas vivos)", fontsize=9); ax.legend(fontsize=7)
ax = axs[1]; t = bt.sort_values("rank_pos").head(20).iloc[::-1]; y = np.arange(len(t))
ax.scatter(t.rank_pos, y, c="#4c72b0", s=30, label="UHI 5σ", zorder=3)
ax.scatter(t.planet.map(r3.rank_pos), y, c="#dd8452", marker="^", s=30, label="UHI 3σ", zorder=3)
ax.scatter(t.ESI_rank_SM11, y, c="#55a868", marker="s", s=26, label="ESI (SM11, este catálogo)", zorder=3)
if "PHL_ESI_rank" in t:
    ax.scatter(t.PHL_ESI_rank, y, c="#c44e52", marker="D", s=24, label="ESI PHL/HWC (emparejados)", zorder=3)
for yi, (_, q) in zip(y, t.iterrows()):
    xs = [q.rank_pos, r3.rank_pos.get(q.planet, np.nan), q.ESI_rank_SM11, q.get("PHL_ESI_rank", np.nan)]
    xs = [x for x in xs if pd.notna(x)]; ax.plot([min(xs), max(xs)], [yi, yi], c="#cccccc", lw=0.8, zorder=1)
ax.set_xscale("log"); ax.set_yticks(y); ax.set_yticklabels([f"{int(p)}. {n}" for p, n in zip(t.rank_pos, t.planet)], fontsize=7)
ax.set_xlabel("posición en cada ranking (escala log; 1 = mejor)"); ax.set_title("Top 20 UHI v5.2 (5σ): posición en UHI 3σ, ESI y PHL", fontsize=9)
ax.legend(fontsize=7, loc="lower right"); fig.tight_layout(); save(fig, "fig_c_uhi_vs_esi_phl")
print("figuras:", sorted(os.listdir(OUT)))
