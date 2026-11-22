"""Resumen y figura del Sobol v6 (v6.1, fix 5).  python3 src/sobol_report_v6.py DIR_5S DIR_3S TAG
Cambios respecto a v6.0 (que promediaba salidas degeneradas y daba S_T de hasta 9e4):
  * se excluyen salidas constantes, no finitas, "casi siempre iguales" (< MIN_FRAC de evaluaciones distintas del valor modal;
    p. ej. el UHI Hycean de un planeta casi siempre 0) y las de estimador inestable (S_T puntual > ST_MAX);
  * S_T (Jansen) y S1 (Saltelli 2010) se recalculan sobre la matriz A, AB_i, B de SALib, y el IC95 de la MEDIA por rama
    sale de un bootstrap (B_BOOT remuestreos de las N filas base, comunes a todas las salidas);
  * se reportan ancho completo y semiancho del IC95 y el criterio de cierre (< 0.05) con ambas lecturas.
Escribe tables/sobol_<TAG>_<rama>_<t>s.csv, tables/sobol_<TAG>_closing.json, tables/sobol_<TAG>_outputs.csv y
figures/fig_<TAG>_sobol_ST.png (no sobreescribe las salidas de v6.0 si TAG es distinto)."""
import os, sys, json, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(ROOT, "src"))
MIN_FRAC, ST_MAX, B_BOOT, THR_CI = 0.10, 1.5, 500, 0.05


def jansen(A, Bm, AB):
    V = np.var(np.r_[A, Bm], axis=0)
    with np.errstate(all="ignore"):
        ST = 0.5 * np.mean((A[:, None, :] - AB) ** 2, axis=0) / V
        S1 = np.mean(Bm[:, None, :] * (AB - A[:, None, :]), axis=0) / V
    return ST, S1


def degenerate(Y, ST):
    out = []
    for k in range(Y.shape[1]):
        y = Y[:, k]
        if not np.isfinite(y).all() or np.std(y) == 0:
            out.append("constante"); continue
        vals, cnt = np.unique(y, return_counts=True); frac = 1 - cnt.max() / len(y)
        if frac < MIN_FRAC: out.append(f"casi_siempre_igual({frac:.3f})")
        elif not np.isfinite(ST[:, k]).all() or np.nanmax(ST[:, k]) > ST_MAX: out.append(f"S_T_inestable({np.nanmax(ST[:, k]):.3g})")
        else: out.append("")
    return out


def analyze(Y, names, D, branches=("rocky", "hycean", "combined"), seed=2026):
    N = Y.shape[0] // (D + 2); Yr = Y[: N * (D + 2)].reshape(N, D + 2, -1)
    A, Bm, AB = Yr[:, 0], Yr[:, -1], Yr[:, 1:-1]
    ST, S1 = jansen(A, Bm, AB); flags = degenerate(Y, ST)
    rng = np.random.default_rng(seed); idx = [rng.integers(0, N, N) for _ in range(B_BOOT)]
    boots = np.stack([jansen(A[i], Bm[i], AB[i])[0] for i in idx])
    res = {}
    for br in branches:
        ks = [k for k, on in enumerate(names) if on.startswith(br + ":") and not flags[k]]
        m = ST[:, ks].mean(1); lo, hi = np.percentile(boots[:, :, ks].mean(2), [2.5, 97.5], axis=0)
        res[br] = dict(ks=ks, m=m, lo=lo, hi=hi, s1=S1[:, ks].mean(1))
    return res, flags, N


if __name__ == "__main__":
    import sens_v6 as S
    d5, d3, tag = (sys.argv[1:4] if len(sys.argv) >= 4 else ("sens_out_v6_5s", "sens_out_v6_3s", "v6_corr"))
    names, P = S.out_names(), S.NAMES; D = len(P)
    fig, axs = plt.subplots(2, 3, figsize=(18, 11)); closing = {}; orows = []
    for i, (t, dd) in enumerate(((5, d5), (3, d3))):
        Y = np.load(os.path.join(ROOT, dd, "Y_outputs.npy"))
        res, flags, N = analyze(Y, names, D)
        orows += [dict(thr=t, output=on, excluida=bool(f), motivo=f) for on, f in zip(names, flags)]
        for j, br in enumerate(("rocky", "hycean", "combined")):
            r = res[br]; o = np.argsort(-r["m"])
            tab = pd.DataFrame(dict(rank=np.arange(1, D + 1), param=np.array(P)[o], ST_mean=r["m"][o], ci_lo=r["lo"][o], ci_hi=r["hi"][o], S1_mean=r["s1"][o]))
            tab["ci_full_width"] = tab.ci_hi - tab.ci_lo; tab["ci_half_width"] = tab.ci_full_width / 2
            tab["pass_full"] = tab.ci_full_width < THR_CI; tab["pass_half"] = tab.ci_half_width < THR_CI
            tab.to_csv(os.path.join(ROOT, "tables", f"sobol_{tag}_{br}_{t}s.csv"), index=False)
            top = tab.head(10)
            closing[f"{br}_{t}s"] = dict(N=int(N), n_outputs=len(r["ks"]), n_excluded=int(sum(1 for k, on in enumerate(names) if on.startswith(br + ":") and flags[k])),
                                         top10_pass_full_width=int(top.pass_full.sum()), top10_pass_half_width=int(top.pass_half.sum()),
                                         criterion_full_met=bool(top.pass_full.all()), criterion_half_met=bool(top.pass_half.all()),
                                         top10=top[["param", "ST_mean", "ci_lo", "ci_hi", "ci_full_width", "ci_half_width"]].round(4).values.tolist())
            s = tab.head(15)[::-1]; ax = axs[i, j]
            ax.barh(s.param, s.ST_mean, xerr=[s.ST_mean - s.ci_lo, s.ci_hi - s.ST_mean], color=np.where(s.pass_full, "#31a354", "#de2d26"), alpha=.85, capsize=2)
            ax.barh(s.param, s.S1_mean.clip(lower=0), color="#08306b", height=.35); ax.tick_params(axis="y", labelsize=8)
            ax.set_title(f"{br} {t}σ (N={N}; {len(r['ks'])} salidas): S_T medio, IC95 bootstrap; verde = ancho < 0.05", fontsize=9)
    fig.suptitle(f"UHI {tag}: Sobol (salidas degeneradas excluidas); S1 en azul oscuro"); fig.tight_layout(rect=(0, 0, 1, .97))
    fig.savefig(os.path.join(ROOT, "figures", f"fig_{tag}_sobol_ST.png"), dpi=130)
    pd.DataFrame(orows).to_csv(os.path.join(ROOT, "tables", f"sobol_{tag}_outputs.csv"), index=False)
    json.dump(closing, open(os.path.join(ROOT, "tables", f"sobol_{tag}_closing.json"), "w"), indent=1)
    for k, v in closing.items():
        print(k, {kk: vv for kk, vv in v.items() if kk != "top10"})
