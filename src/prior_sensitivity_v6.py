r"""Class-prior sensitivity of rocky/combined rankings and formal Innes vs Madhusudhan comparison.

Reads uhi_v6_{rocky,combined,hycean}_<variant>_<thr>s.csv from UHI_RES (default results/) and writes
prior_sens_* and innes_vs_madhu_* tables to UHI_TAB (default tables/) plus figures to UHI_FIG
(default figures/). For the v6.1 release layout, set UHI_RES=results/v6.1. These environment
variables only redirect files; they do not change computations.

Scenarios correspond to docs/assumptions_v6_es.md prior tests (a)-(d): rocky/volatile class priors,
alternative mass-route transitions (Chen & Kipping 2017 and empirical fits), true-mass treatment for
M sin i (minimum mass or a Howard et al. 2010 prior), an alternative M-R relation, and alternative
composition boundaries. A planet is robust if it stays in the top 10 in every scenario, stable if it
stays in the top 20 with spread <=10 places, and prior-dependent otherwise; these labels are
reporting conventions. The Hycean comparison treats Innes et al. (2023) as production and
Madhusudhan et al. (2021) as an optimistic variant."""
import os, json, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import spearmanr, kendalltau
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RD, TD, FD = (os.environ.get(k, os.path.join(ROOT, d)) for k, d in (("UHI_RES", "results"), ("UHI_TAB", "tables"), ("UHI_FIG", "figures")))
for _d in (TD, FD): os.makedirs(_d, exist_ok=True)
R = lambda f: os.path.join(RD, f); T = lambda f: os.path.join(TD, f); F = lambda f: os.path.join(FD, f)
# Scenario labels used in tables and figures (keys are run_uhi_v6.py variant names).
SC = {"pA1_rocoso": "A1 prior rocoso (π=2:1:1, M_v 5-10)", "pA2_volatil": "A2 prior volátil (π=0.5:1:1, M_v 1.5-3)",
      "pA3_ck17": "A3 transición CK17 (M_v=2.04)", "pA4_empM": "A4 empírico archivo M (M_v=4.10, s=0.40)",
      "pA5_empM_P5": "A5 empírico M, P>5 d (M_v=2.83, s=0.70)", "pB1_msini_min": "B1 M = M sin i", "pB2_msini_howard": "B2 M sin i + prior Howard+2010",
      "pC1_mr_ck17": "C1 M-R Chen & Kipping 2017", "pD1_linea_earthlike": "D1 frontera línea Earth-like", "pD2_linea_h2o100": "D2 frontera 100% H2O"}
def md_table(d, index=False):
    r"""Render a DataFrame as a Markdown table, formatting floating values to three decimals."""
    d = d.reset_index() if index else d
    fmt = lambda v: f"{v:.3f}" if isinstance(v, (float, np.floating)) else str(v)
    return "| " + " | ".join(map(str, d.columns)) + " |\n|" + "---|" * len(d.columns) + "\n" + "\n".join("| " + " | ".join(fmt(v) for v in r) + " |" for r in d.values)
load = lambda b, v, t: pd.read_csv(R(f"uhi_v6_{b}_{v}_{t}s.csv")).set_index("planet")
out = {}
for t in [5, 3]:
    avail = [v for v in SC if os.path.exists(R(f"uhi_v6_rocky_{v}_{t}s.csv"))]
    if not avail:
        continue
    for b in ["rocky", "combined"]:
        m = load(b, "main", t); runs = {v: load(b, v, t) for v in avail}
        # global metrics per scenario
        g = []
        for v, d in runs.items():
            u = list(dict.fromkeys(list(m.index[:20]) + list(d.index[:20])))
            alive = list(set(m.index[m.P_alive > 0]) | set(d.index[d.P_alive > 0]))
            g.append(dict(escenario=SC[v], top10_overlap=len(set(m.index[:10]) & set(d.index[:10])), top20_overlap=len(set(m.index[:20]) & set(d.index[:20])),
                          kendall_top20_union=kendalltau(m.loc[u, "UHI_rank"], d.loc[u, "UHI_rank"]).statistic,
                          spearman_top20_union=spearmanr(m.loc[u, "UHI_rank"], d.loc[u, "UHI_rank"]).statistic,
                          spearman_vivos=spearmanr(m.loc[alive, "UHI_rank"], d.loc[alive, "UHI_rank"]).statistic, n_vivos=int((d.P_alive > 0).sum()),
                          nuevo_n1=d.index[0]))
        g = pd.DataFrame(g); g.to_csv(T(f"prior_sens_global_{b}_{t}s.csv"), index=False)
        # per planet: main top 20
        rows = []
        for p in m.index[:20]:
            rk = {v: int(runs[v].loc[p, "rank_pos"]) for v in avail}; ur = {v: runs[v].loc[p, "UHI_rank"] for v in avail}
            pt = {v: runs[v].loc[p, "P_top10"] for v in avail}
            worst = max(rk, key=rk.get); best = min(rk, key=rk.get)
            rr = list(rk.values()) + [int(m.loc[p, "rank_pos"])]
            spread = max(rr) - min(rr)
            if max(rr) <= 10: flag = "robusto (top 10 en todos)"
            elif max(rr) <= 20 and spread <= 10: flag = "estable (top 20 en todos)"
            else: flag = "DEPENDE DEL PRIOR"
            rows.append(dict(planet=p, evidencia=m.loc[p, "evidence"], mass_source=m.loc[p, "mass_source"], rank_main=int(m.loc[p, "rank_pos"]),
                             rank_min=min(rr), rank_max=max(rr), peor_escenario=SC[worst].split(" ")[0], mejor_escenario=SC[best].split(" ")[0],
                             UHI_main=m.loc[p, "UHI_rank"], UHI_min=min(list(ur.values()) + [m.loc[p, "UHI_rank"]]), UHI_max=max(list(ur.values()) + [m.loc[p, "UHI_rank"]]),
                             P_R_main=m.loc[p, "P_R"], P_R_min=min(runs[v].loc[p, "P_R"] for v in avail), P_R_max=max(runs[v].loc[p, "P_R"] for v in avail),
                             P_top10_main=m.loc[p, "P_top10"], P_top10_min=min(list(pt.values()) + [m.loc[p, "P_top10"]]), P_top10_max=max(list(pt.values()) + [m.loc[p, "P_top10"]]),
                             veredicto=flag, **{f"rank_{v.split('_')[0]}": rk[v] for v in avail}))
        pp = pd.DataFrame(rows); pp.to_csv(T(f"prior_sens_planets_{b}_{t}s.csv"), index=False)
        open(T(f"prior_sens_{b}_{t}s.md"), "w").write(f"# Sensibilidad al prior de clase - {b}, {t} sigma\n\n## Global (vs principal)\n\n" + md_table(g) +
                                                      "\n\n## Top 20 principal\n\n" + md_table(pp.drop(columns=[c for c in pp.columns if c.startswith("rank_p")])) + "\n")
        out[f"{b}_{t}s"] = dict(global_=g.round(3).to_dict("records"), planets=pp[["planet", "evidencia", "rank_main", "rank_min", "rank_max", "peor_escenario",
                                "UHI_main", "UHI_min", "UHI_max", "P_R_min", "P_R_max", "P_top10_min", "P_top10_max", "veredicto"]].round(3).to_dict("records"))
        if b == "rocky":
            M = np.array([[int(m.loc[p, "rank_pos"])] + [int(runs[v].loc[p, "rank_pos"]) for v in avail] for p in pp.planet])
            fig, ax = plt.subplots(figsize=(13, 9))
            im = ax.imshow(np.minimum(M, 40), cmap="RdYlGn_r", vmin=1, vmax=40, aspect="auto")
            for i in range(M.shape[0]):
                for j in range(M.shape[1]):
                    ax.text(j, i, str(M[i, j]), ha="center", va="center", fontsize=7)
            ax.set_xticks(range(M.shape[1])); ax.set_xticklabels(["principal"] + [SC[v] for v in avail], rotation=40, ha="right", fontsize=8)
            ax.set_yticks(range(M.shape[0])); ax.set_yticklabels([f"{p} [{e}] - {f.split(' ')[0]}" for p, e, f in zip(pp.planet, pp.evidencia, pp.veredicto)], fontsize=8)
            plt.colorbar(im, ax=ax, label="posición (tope 40)"); ax.set_title(f"Rocky {t}σ: posición del top 20 principal bajo cada prior de clase")
            fig.tight_layout(); fig.savefig(F(f"fig_v6_prior_heatmap_{t}s.png"), dpi=140); plt.close(fig)
            fig, ax = plt.subplots(figsize=(9, 7)); y = np.arange(len(pp))[::-1]
            ax.hlines(y, pp.P_top10_min, pp.P_top10_max, color="grey", lw=3); ax.plot(pp.P_top10_main, y, "ko", label="principal")
            ax.set_yticks(y); ax.set_yticklabels(pp.planet, fontsize=8); ax.set_xlabel("P_top10 (rango sobre los 10 escenarios de prior)")
            ax.set_title(f"Rocky {t}σ: dispersión de P_top10 por prior de clase"); ax.legend(); fig.tight_layout(); fig.savefig(F(f"fig_v6_prior_ptop10_{t}s.png"), dpi=140); plt.close(fig)

# ---- Innes vs Madhusudhan (Hycean), formal comparison ----
rows = []; md = ["# Rama Hycean: Innes+2023 (principal) vs Madhusudhan+2021 (variante optimista)\n"]
# Hycean candidates and benchmarks reported individually in the Innes/Madhusudhan comparison.
KEYP = ["K2-18 b", "TOI-7166 b", "K2-332 b", "K2-125 b", "TOI-2094 b", "TOI-4336 A b", "LHS 1140 b", "TOI-2257 b", "TOI-270 d", "Kepler-1630 b", "Wolf 1061 d", "Kapteyn c"]
for t in [5, 3]:
    i_, m_ = load("hycean", "main", t), load("hycean", "madhu", t)
    ci, cm = load("combined", "main", t), load("combined", "madhu", t)
    alive = list(set(i_.index[i_.P_alive > 0]) | set(m_.index[m_.P_alive > 0]))
    rows.append(dict(thr=t, vivos_Innes=int((i_.P_alive > 0).sum()), vivos_Madhu=int((m_.P_alive > 0).sum()), top10_overlap=len(set(i_.index[:10]) & set(m_.index[:10])),
                     spearman_vivos=spearmanr(i_.loc[alive, "UHI_rank"], m_.loc[alive, "UHI_rank"]).statistic,
                     UHI_max_Innes=i_.UHI_rank.iloc[0], UHI_max_Madhu=m_.UHI_rank.iloc[0], n_Hycean_en_top10_combinado_Innes=int(sum(ci.loc[p, "P_H"] > ci.loc[p, "P_R"] for p in ci.index[:10])),
                     n_Hycean_en_top10_combinado_Madhu=int(sum(cm.loc[p, "P_H"] > cm.loc[p, "P_R"] for p in cm.index[:10])),
                     combinado_top10_overlap=len(set(ci.index[:10]) & set(cm.index[:10]))))
    side = pd.DataFrame({"pos": range(1, 11), "Innes: planeta": i_.index[:10], "Innes: UHI_rank": i_.UHI_rank.values[:10], "Innes: P_alive": i_.P_alive.values[:10],
                         "Madhu: planeta": m_.index[:10], "Madhu: UHI_rank": m_.UHI_rank.values[:10], "Madhu: P_alive": m_.P_alive.values[:10]})
    key = pd.DataFrame([dict(planeta=p, P_H=i_.loc[p, "P_H"], S_med=i_.loc[p, "S_med"], Innes_rank=int(i_.loc[p, "rank_pos"]), Innes_UHI=i_.loc[p, "UHI_rank"],
                             Innes_P_clima_H=i_.loc[p, "P_clima_H"], Madhu_rank=int(m_.loc[p, "rank_pos"]), Madhu_UHI=m_.loc[p, "UHI_rank"],
                             Madhu_P_clima_H=m_.loc[p, "P_clima_H"], S_atm_H_med=m_.loc[p, "tloss_all_med"]) for p in KEYP if p in i_.index])
    side.to_csv(T(f"innes_vs_madhu_top10_{t}s.csv"), index=False); key.to_csv(T(f"innes_vs_madhu_key_{t}s.csv"), index=False)
    md += [f"\n## {t} sigma - top 10 lado a lado\n\n", md_table(side), f"\n\n## {t} sigma - planetas clave (t_loss en Gyr, mediana)\n\n", md_table(key.rename(columns={"S_atm_H_med": "t_loss_med_Gyr"})), "\n"]
    out[f"innes_madhu_{t}s"] = dict(side=side.round(3).values.tolist(), key=key.round(3).to_dict("records"))
summ = pd.DataFrame(rows); summ.to_csv(T("innes_vs_madhu_summary.csv"), index=False)
md.insert(1, "\n## Resumen\n\n" + md_table(summ) + "\n")
open(T("innes_vs_madhu.md"), "w").write("".join(md)); out["innes_madhu_summary"] = summ.round(3).to_dict("records")
fig, axs = plt.subplots(1, 2, figsize=(14, 5.5))
for ax, (nm, d) in zip(axs, [("Innes+2023 (principal)", load("hycean", "main", 5)), ("Madhusudhan+2021 (variante)", load("hycean", "madhu", 5))]):
    dd = d.head(10)[::-1]; ax.barh(dd.index, dd.UHI_rank, color="#3182bd"); ax.set_title(f"Hycean 5σ - {nm}"); ax.set_xlabel("UHI_rank")
    for yy, (v, pa) in enumerate(zip(dd.UHI_rank, dd.P_alive)): ax.text(v, yy, f" {v:.3f} (P_alive {pa:.2f})", va="center", fontsize=7)
    ax.set_xlim(0, max(dd.UHI_rank.max() * 1.7, 0.01)); ax.tick_params(axis="y", labelsize=8)
fig.tight_layout(); fig.savefig(F("fig_v6_innes_vs_madhu_5s.png"), dpi=140); plt.close(fig)
json.dump(out, open(T("prior_sens_summary.json"), "w"), indent=1, default=float)
print(json.dumps(out.get("innes_madhu_summary"), default=float))
