r"""Summarize the released v6.1 Sobol samples and draw the S_T figure (v6.1 fix 5).

Usage: python3 src/sobol_report_v6.py DIR_5S DIR_3S TAG, with directories relative to the repository
root (for the release: sens_out_v6_5s_n512 sens_out_v6_3s_n512 v6.1_n512). Only each Y_outputs.npy is
read; the model is not rerun.

Differences from v6.0: outputs that are non-finite, constant, nearly always equal (fewer than
MIN_FRAC of evaluations differ from the modal value; for example a Hycean UHI that is almost always
zero), or unstable (any point S_T above ST_MAX) are excluded rather than averaged. Total effects use
Jansen (1999) and first-order effects Saltelli et al. (2010), computed from SALib's A, AB_i, B blocks.
Branch-mean S_T 95% intervals use B_BOOT bootstrap resamples of the N base rows, shared across outputs.
The closing criterion is met when every top-10 parameter has an interval width below THR_CI; both full
width and half width are reported because the stricter interpretation should remain visible.

Writes tables/sobol_<TAG>_<branch>_<thr>s.csv, tables/sobol_<TAG>_closing.json,
tables/sobol_<TAG>_outputs.csv and figures/fig_<TAG>_sobol_ST.png. A repeated TAG overwrites those files."""
import os, sys, json, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(ROOT, "src"))
# Analysis settings: minimum non-modal fraction, maximum stable point S_T, bootstrap resamples, and the 0.05
# interval-width closing threshold.
MIN_FRAC, ST_MAX, B_BOOT, THR_CI = 0.10, 1.5, 500, 0.05


def jansen(A, Bm, AB):
    r"""Return Jansen total-effect and Saltelli first-order indices.

    A and Bm have shape (N, outputs); AB has shape (N, D, outputs). Variances combine A and B. Outputs
    have shape (D, outputs) and are dimensionless. Division-by-zero warnings are suppressed because
    degenerate outputs are flagged separately."""
    V = np.var(np.r_[A, Bm], axis=0)
    with np.errstate(all="ignore"):
        ST = 0.5 * np.mean((A[:, None, :] - AB) ** 2, axis=0) / V
        S1 = np.mean(Bm[:, None, :] * (AB - A[:, None, :]), axis=0) / V
    return ST, S1


def degenerate(Y, ST):
    r"""Return one exclusion reason per output column ('' means retained).

    Y is the full (evaluations, outputs) matrix and ST the point estimates. Thresholds MIN_FRAC=0.10 and
    ST_MAX=1.5 are documented analysis choices (v6.1 fix 5), not properties of SALib."""
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
    r"""Compute retained branch-mean S_T, bootstrap intervals and mean S1.

    Y is reshaped as (N, D+2, outputs), names gives labels, D is number of factors and seed fixes the
    bootstrap. Returns (per-branch dictionary, exclusion reasons, N). Retained outputs are equally weighted."""
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
