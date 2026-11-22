r"""LEGACY AUDIT SCRIPT (v6.0-to-v6.1 review). Kept as a record; it is not rerunnable from the release
as shipped because its input paths or inputs are not included (see README, "Legacy scripts").

Build the Phase-1 tables and figures from saved Monte Carlo results; it does not rerun the model.

Inputs are read from results/ in the repository root: uhi_v6_<branch>_<variant>_<thr>s.csv for main
and variant runs, summary_<variant>_5s.json, and the v5.2 reference tables in data/. The v6.1
release stores main/madhu/prior-test CSVs under results/v6.1 and does not ship the teff, bt3d, jackson,
only_clim, only_xuv, only_d3 and single variant runs. Regenerating these tables therefore requires
running those variants with --outdir results (5 and 3 sigma as used below); the shipped tables are
already generated.

Outputs (tables/ and figures/): top-10 lists per branch, v5.2-to-v6 rank shifts with D1/D2/D3
ablation attribution, Spearman and top-N overlap comparisons, class probabilities by evidence route,
summary_v6.json, top-10, rank-shift and class-probability figures, and the D2 Pass/Mors continuity
check. Branch top-10 lists here are pooled across facilities, not split into JWST and ELT lists."""
import os, sys, json, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import spearmanr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = lambda f: os.path.join(ROOT, "results", f); T = lambda f: os.path.join(ROOT, "tables", f); F = lambda f: os.path.join(ROOT, "figures", f)
# Released v5.2 catalogue rankings at 5 and 3 sigma, kept only as comparison baselines.
# Released v5.2 catalogue rankings at 5 and 3 sigma, kept only as comparison baselines.
V52 = {5: pd.read_csv(os.path.join(ROOT, "data", "uhi_results_v52.csv")), 3: pd.read_csv(os.path.join(ROOT, "data", "uhi_results_v52_3sigma.csv"))}
BR = ["rocky", "hycean", "combined"]
# Well-known targets always included in the class-probability table.
# Well-known targets always included in the class-probability table.
ANCH = ["TRAPPIST-1 d", "TRAPPIST-1 e", "Proxima Cen b", "Proxima Cen d", "Ross 128 b", "K2-125 b", "TOI-4336 A b", "TOI-7166 b",
        "LHS 1140 b", "K2-18 b", "TOI-700 d", "TOI-700 e"]
load = lambda b, v, t: pd.read_csv(R(f"uhi_v6_{b}_{v}_{t}s.csv"))
COLS = ["rank_pos", "planet", "UHI_rank", "P_alive", "UHI_cond_med", "P_top10", "rank_med", "P_R", "P_H", "P_G", "evidence", "tier", "st_teff", "S_med"]
EVC = {"densidad": "#1b9e77", "radio": "#7570b3", "prior": "#d95f02"}
out = {}


def md_table(d, index=False):
    r"""Render a DataFrame as a Markdown table, formatting floating values to three decimals."""
    d = d.reset_index() if index else d
    fmt = lambda v: f"{v:.3f}" if isinstance(v, (float, np.floating)) else str(v)
    return "| " + " | ".join(map(str, d.columns)) + " |\n|" + "---|" * len(d.columns) + "\n" + "\n".join("| " + " | ".join(fmt(v) for v in r) + " |" for r in d.values)


# ---- T1: top 10 per branch and threshold ----
md = ["# UHI v6 Fase 1 - Top 10 por rama (10,000 draws, semilla 42)\n"]
for t in (5, 3):
    for b in BR:
        d = load(b, "main", t).head(10)[COLS]
        d.to_csv(T(f"top10_{b}_{t}s.csv"), index=False)
        md.append(f"\n## {b} - {t} sigma\n\n" + md_table(d) + "\n")
        out[f"top10_{b}_{t}s"] = d[["planet", "UHI_rank", "P_alive", "UHI_cond_med", "P_top10", "P_R", "P_H", "evidence"]].round(3).values.tolist()
open(T("top10_v6.md"), "w").write("".join(md))

# ---- T2: combined rank changes vs v5.2, with attribution by ablation (ablations at 5 sigma only) ----
abl = {"solo_D1_clima": "only_clim", "solo_D2_XUV": "only_xuv", "solo_D3": "only_d3", "D1+D2+D3_sin_clases": "single"}
rows = []
for t in (5, 3):
    c6 = load("combined", "main", t).set_index("planet"); v5 = V52[t].set_index("planet")
    top = list(dict.fromkeys(list(V52[t].planet.head(20)) + list(c6.index[:20])))
    for p in top:
        r = dict(thr=t, planet=p, rank_v52=int(v5.loc[p, "rank_pos"]), UHI_rank_v52=v5.loc[p, "UHI_rank"], rank_v6=int(c6.loc[p, "rank_pos"]),
                 UHI_rank_v6=c6.loc[p, "UHI_rank"], P_R=c6.loc[p, "P_R"], P_H=c6.loc[p, "P_H"], P_G=c6.loc[p, "P_G"], evidence=c6.loc[p, "evidence"],
                 Rp_med=v5.loc[p, "Rp_med"], ARM_med=load("rocky", "main", t).set_index("planet").loc[p, "ARM_all_med"])
        if t == 5:
            for k, v in abl.items():
                a = load("combined", v, 5).set_index("planet"); r["rank_" + k] = int(a.loc[p, "rank_pos"]); r["UHI_" + k] = a.loc[p, "UHI_rank"]
        r["delta"] = r["rank_v52"] - r["rank_v6"]
        rows.append(r)
sh = pd.DataFrame(rows); sh.to_csv(T("rank_shifts_v52_v6.csv"), index=False)
open(T("rank_shifts_v52_v6.md"), "w").write(md_table(sh))

# ---- T3: Spearman correlations and top-N overlaps ----
def comp(a, b, col="UHI_rank"):
    r"""Compare two ranked result tables on a shared score column.

    Returns Spearman rho for all planets and for planets alive in either run, the number in that union,
    and top-10/top-20 name overlaps. Inputs must contain planet and the score column (default UHI_rank)."""
    m = a[["planet", col]].merge(b[["planet", col]], on="planet")
    alive = (m[col + "_x"] > 0) | (m[col + "_y"] > 0)
    return dict(spearman_all=spearmanr(m[col + "_x"], m[col + "_y"]).statistic,
                spearman_alive_union=spearmanr(m.loc[alive, col + "_x"], m.loc[alive, col + "_y"]).statistic,
                n_alive_union=int(alive.sum()),
                top10_overlap=len(set(a.planet.head(10)) & set(b.planet.head(10))),
                top20_overlap=len(set(a.planet.head(20)) & set(b.planet.head(20))))
sp = []
for t in (5, 3):
    for b in BR:
        sp.append(dict(comp=f"v6 {b} vs v5.2", thr=t, **comp(load(b, "main", t), V52[t])))
for b in BR:
    sp.append(dict(comp=f"v6 {b}: 5 sigma vs 3 sigma", thr=0, **comp(load(b, "main", 5), load(b, "main", 3))))
sp.append(dict(comp="v5.2: 5 sigma vs 3 sigma", thr=0, **comp(V52[5], V52[3])))
for v in ["teff", "bt3d", "jackson", "madhu", "only_clim", "only_xuv", "only_d3", "single"]:
    for b in (BR if v not in ("only_clim", "only_xuv", "only_d3", "single") else ["combined"]):
        sp.append(dict(comp=f"variante {v} vs main ({b})", thr=5, **comp(load(b, v, 5), load(b, "main", 5))))
sp.append(dict(comp="variante madhu vs main (hycean)", thr=3, **comp(load("hycean", "madhu", 3), load("hycean", "main", 3))))
sp = pd.DataFrame(sp); sp.to_csv(T("spearman_v6.csv"), index=False); open(T("spearman_v6.md"), "w").write(md_table(sp))

# ---- T4: class probabilities with an evidence column ----
c5 = load("combined", "main", 5).set_index("planet"); r5 = load("rocky", "main", 5).set_index("planet"); h5 = load("hycean", "main", 5).set_index("planet")
hm = load("hycean", "madhu", 5).set_index("planet")
pl = list(dict.fromkeys(ANCH + list(r5.index[:10]) + list(h5.index[:10]) + list(hm.index[:10]) + list(V52[5].planet.head(10))))
cp = pd.DataFrame({"evidencia": c5.loc[pl, "evidence"], "P_R": c5.loc[pl, "P_R"], "P_H": c5.loc[pl, "P_H"], "P_G": c5.loc[pl, "P_G"],
                   "Rp_cat": c5.loc[pl, "Rp_med"], "Mp_cat": c5.loc[pl, "Mp_med"], "mass_source": c5.loc[pl, "mass_source"],
                   "dens_ok": c5.loc[pl, "dens_ok"], "rank_rocky_5s": r5.loc[pl, "rank_pos"], "rank_hycean_5s": h5.loc[pl, "rank_pos"],
                   "rank_hycean_madhu_5s": hm.loc[pl, "rank_pos"], "rank_comb_5s": c5.loc[pl, "rank_pos"]}).rename_axis("planet")
cp.to_csv(T("class_probabilities_v6.csv")); open(T("class_probabilities_v6.md"), "w").write(md_table(cp, index=True))
ev = c5.groupby("evidence")[["P_R", "P_H", "P_G"]].agg(["mean", "count"]); ev.to_csv(T("class_by_evidence_v6.csv"))

# ---- summary ----
summ = {}
for t in (5, 3):
    for b in BR:
        d = load(b, "main", t); summ[f"{b}_{t}s"] = dict(scored=len(d), alive=int((d.P_alive > 0).sum()))
    summ[f"v52_{t}s"] = dict(scored=len(V52[t]), alive=int((V52[t].P_alive > 0).sum()))
for v in ["teff", "bt3d", "jackson", "madhu", "only_clim", "only_xuv", "only_d3", "single"]:
    s = json.load(open(R(f"summary_{v}_5s.json"))); summ[f"{v}_5s"] = {k: s[k] for k in s if k in BR}
summ["global_class_mean"] = c5[["P_R", "P_H", "P_G"]].mean().round(4).to_dict()
summ["class_by_evidence_counts"] = c5.evidence.value_counts().to_dict()
json.dump(dict(summary=summ, tops=out), open(T("summary_v6.json"), "w"), indent=1, default=float)

# ---- figures ----
def topfig(t):
    r"""Draw the three-branch top-10 bar figure for threshold t (5 or 3 sigma).

    Bars are UHI_rank coloured by the class-evidence route; labels show P_R, P_H and P_top10."""
    fig, ax = plt.subplots(1, 3, figsize=(16, 5.2))
    for k, b in enumerate(BR):
        d = load(b, "main", t).head(10)[::-1]
        y = np.arange(len(d))
        ax[k].barh(y, d.UHI_rank, color=[EVC.get(e, "grey") for e in d.evidence])
        ax[k].errorbar(d.P_alive * d.UHI_cond_q25, y, fmt="none") if False else None
        ax[k].set_yticks(y); ax[k].set_yticklabels([f"{p}  (P_R={r:.2f}, P_H={h:.2f})" for p, r, h in zip(d.planet, d.P_R, d.P_H)], fontsize=8)
        for yy, v, pt in zip(y, d.UHI_rank, d.P_top10):
            ax[k].text(v, yy, f" {v:.3f} | P10={pt:.2f}", va="center", fontsize=7)
        ax[k].set_title(f"{b} ({t}σ)"); ax[k].set_xlabel("UHI_rank = P_alive × mediana condicional")
        ax[k].set_xlim(0, max(d.UHI_rank.max() * 1.6, 0.01))
    from matplotlib.patches import Patch
    fig.legend(handles=[Patch(color=c, label=f"evidencia: {e}") for e, c in EVC.items()], loc="lower center", ncol=3, fontsize=9)
    fig.suptitle(f"UHI v6 Fase 1 - top 10 por rama, umbral {t}σ (Innes+2023 para Hycean)"); fig.tight_layout(rect=(0, 0.05, 1, 0.96))
    fig.savefig(F(f"fig_v6_top10_{t}s.png"), dpi=150); plt.close(fig)
for t in (5, 3):
    topfig(t)

for t in (5, 3):
    c6 = load("combined", "main", t).set_index("planet"); v5 = V52[t].set_index("planet")
    top = list(dict.fromkeys(list(V52[t].planet.head(15)) + list(c6.index[:15])))
    fig, ax = plt.subplots(figsize=(8, 10))
    far = []
    for p in top:
        a, b = v5.loc[p, "rank_pos"], c6.loc[p, "rank_pos"]
        col = "#c0392b" if b > a + 5 else ("#27ae60" if b < a - 5 else "#7f8c8d")
        ax.plot([0, 1], [a, min(b, 60)], "-o", color=col, lw=1.5)
        ax.text(-0.03, a, f"{p} #{a}", ha="right", va="center", fontsize=7.5)
        if b > 60:
            far.append(f"{p}: #{b} (P_R={c6.loc[p, 'P_R']:.2f})")
        else:
            ax.text(1.03, b, f"#{b} {p}", ha="left", va="center", fontsize=7.5)
    if far:
        ax.text(1.03, 59, "fuera del top 60 en v6:\n" + "\n".join(far), ha="left", va="bottom", fontsize=7, color="#c0392b")
    ax.set_ylim(62, 0); ax.set_xlim(-0.6, 1.6); ax.set_xticks([0, 1]); ax.set_xticklabels(["v5.2", "v6 combinado"]); ax.set_ylabel("posicion")
    ax.set_title(f"Cambios de rango v5.2 → v6 (combinado, {t}σ)"); fig.tight_layout(); fig.savefig(F(f"fig_v6_rank_shifts_{t}s.png"), dpi=150); plt.close(fig)

cpl = cp.sort_values(["evidencia", "P_R"])
fig, ax = plt.subplots(figsize=(9, max(6, 0.28 * len(cpl))))
y = np.arange(len(cpl))
ax.barh(y, cpl.P_R, color="#8c6d31", label="R (rocoso / mundo de agua sin H2)")
ax.barh(y, cpl.P_H, left=cpl.P_R, color="#3182bd", label="H (Hycean)")
ax.barh(y, cpl.P_G, left=cpl.P_R + cpl.P_H, color="#bdbdbd", label="G (mini-Neptuno rico en gas)")
ax.set_yticks(y); ax.set_yticklabels([f"{p} [{e}]" for p, e in zip(cpl.index, cpl.evidencia)], fontsize=7)
for yy, e in zip(y, cpl.evidencia):
    ax.plot(-0.02, yy, "s", color=EVC[e], ms=5, clip_on=False)
ax.set_xlim(0, 1); ax.set_xlabel("probabilidad de clase (priors iguales, D4)"); ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, fontsize=8)
fig.tight_layout(); fig.savefig(F("fig_v6_class_probabilities.png"), dpi=150); plt.close(fig)

# D2 validation: Pass (fully convective) vs Mors fluence across the 0.35 M_sun blend
sys.path.insert(0, os.path.join(ROOT, "src")); import xuv_history as X
g = X.mors_grid(); ms = np.linspace(0.1, 0.6, 51); fc = []; mo = {p: [] for p in (5, 50, 95)}
for m in ms:
    L = np.interp(5.0, g["ages"], g["Lbol"][np.argmin(abs(g["masses"] - m)), 4]) / 3.828e33
    fc.append(X.fluence_fc(np.array([L]), np.array([1.0]), np.array([5.0]), np.array([m]), -3.05, X.tsat_pass(np.array([m])), 1.25)[0][0] / L)
    for p in mo:
        mo[p].append(X.fluence_mors(np.array([L]), np.array([1.0]), np.array([5.0]), np.array([m]), np.array([float(p)]))[0][0] / L)
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.semilogy(ms, np.array(fc) / 2e18, "k-", label="Pass+2025 (FC, f_flare=1.25)")
for p, ls in zip(mo, [":", "-", "--"]):
    ax.semilogy(ms, np.array(mo[p]) / 2e18, ls, color="#d95f02", label=f"Johnstone+2021 / Mors, percentil {p}")
ax.axvline(0.35, color="grey", lw=0.8); ax.text(0.352, ax.get_ylim()[0] * 3, "mezcla σ((M-0.35)/0.03)", fontsize=8)
ax.set_xlabel("M* [Msun]"); ax.set_ylabel("fluencia XUV acumulada a 5 Gyr / I_E (por L_bol = 1 Lsun, a = 1 AU)")
ax.legend(fontsize=8); ax.set_title("D2: continuidad Pass+2025 ↔ Mors en 0.35 Msun"); fig.tight_layout(); fig.savefig(F("fig_v6_xuv_blend.png"), dpi=150); plt.close(fig)
print("ok", json.dumps(summ, default=float)[:2000])
