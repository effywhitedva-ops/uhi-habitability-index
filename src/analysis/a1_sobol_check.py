r"""LEGACY AUDIT SCRIPT (v6.0-to-v6.1 review). Kept as a record; it is not rerunnable from the release
as shipped because its input paths or inputs are not included (see README, "Legacy scripts").

Independent Phase-1 audit of the Sobol output (written during the v6.0 to v6.1 review).

Checks: (1) degenerate outputs; (2) robust branch-mean Jansen S_T after excluding them; and (3) the
closing criterion that the top-10 S_T 95% intervals have width (and half-width) below 0.05.
Bootstrap: B=500 resamples of the N base rows with seed 2026, shared across outputs.

Run from the repository root. Paths are relative; it reads sens_out_v6_{5,3}s/Y_outputs.npy, i.e.
the earlier N directories. The release N=512 data are in sens_out_v6_{5,3}s_n512 and are summarized
by src/sobol_report_v6.py, which implements the same criterion as the production report. Writes
results/analysis_phase1/sobol_robust_*.csv, sobol_output_diagnostics.csv and sobol_closing_criterion.json."""
import numpy as np, pandas as pd, json, sys, os
sys.path.insert(0, "src"); import sens_v6 as S
OUT = "results/analysis_phase1"; D = len(S.NAMES); names = S.out_names()
rng = np.random.default_rng(2026); B = 500
res = {}; rows = []
for t in (5, 3):
    Y = np.load(f"sens_out_v6_{t}s/Y_outputs.npy"); N = Y.shape[0] // (D + 2)
    Yr = Y.reshape(N, D + 2, -1); A = Yr[:, 0]; Bm = Yr[:, -1]; AB = Yr[:, 1:-1]   # (N,out),(N,out),(N,D,out)
    def st(idx):
        r"""Return (Jansen S_T, Saltelli S1) for selected base rows idx, each with shape (D, outputs)."""
        a, b, ab = A[idx], Bm[idx], AB[idx]
        V = np.var(np.r_[a, b], axis=0)
        num = 0.5 * np.mean((a[:, None, :] - ab) ** 2, axis=0)            # (D,out)
        s1 = np.mean(b[:, None, :] * (ab - a[:, None, :]), axis=0)
        with np.errstate(all="ignore"):
            return num / V, s1 / V
    ST0, S10 = st(np.arange(N))
    boots = np.stack([st(rng.integers(0, N, N))[0] for _ in range(B)])   # (B,D,out)
    # degeneracy diagnostic: fraction of evaluations whose output differs from the modal base value (implemented as the fraction of nonzero outputs, i.e. UHI>0)
    fnz = (np.abs(Y) > 0).mean(0); sd = Y.std(0)
    deg = []
    for k, on in enumerate(names):
        bad = (sd[k] == 0) or (not np.isfinite(ST0[:, k]).all()) or (np.nanmax(ST0[:, k]) > 1.5) or (fnz[k] < 0.10)
        deg.append(bad)
        rows.append(dict(thr=t, output=on, frac_nonzero=round(float(fnz[k]), 4), sd=float(sd[k]), max_ST=float(np.nanmax(ST0[:, k])) if sd[k] > 0 else np.nan,
                         degenerate=bad))
    deg = np.array(deg)
    for br in ["rocky", "hycean", "combined"]:
        ks = [k for k, on in enumerate(names) if on.startswith(br + ":") and not deg[k]]
        m = ST0[:, ks].mean(1); mb = boots[:, :, ks].mean(2)
        lo, hi = np.percentile(mb, [2.5, 97.5], axis=0)
        order = np.argsort(-m)
        tab = pd.DataFrame(dict(param=np.array(S.NAMES)[order], ST_mean=m[order], ci_lo=lo[order], ci_hi=hi[order],
                                ci_width=(hi - lo)[order], S1_mean=S10[:, ks].mean(1)[order]))
        tab["half_width"] = tab.ci_width / 2
        tab["pass_width_lt_0.05"] = tab.ci_width < 0.05
        tab["pass_halfwidth_lt_0.05"] = tab.half_width < 0.05
        tab.insert(0, "rank", np.arange(1, D + 1))
        tab.to_csv(f"{OUT}/sobol_robust_{br}_{t}s.csv", index=False)
        # per-output criterion: SALib-like half-width (1.96*bootstrap sd) over the top 10 of each output
        perout = []
        for k in ks:
            sdk = boots[:, :, k].std(0) * 1.96; top = np.argsort(-ST0[:, k])[:10]
            perout.append((sdk[top] < 0.05).mean())
        res[f"{br}_{t}s"] = dict(n_outputs_used=len(ks), n_degenerate=int(sum(1 for k, on in enumerate(names) if on.startswith(br + ':') and deg[k])),
                                 top10=tab.head(10)[["param", "ST_mean", "ci_lo", "ci_hi", "ci_width"]].round(4).values.tolist(),
                                 top10_pass_full_width=int(tab.head(10)["pass_width_lt_0.05"].sum()),
                                 top10_pass_half_width=int(tab.head(10)["pass_halfwidth_lt_0.05"].sum()),
                                 per_output_frac_top10_halfwidth_lt005=round(float(np.mean(perout)), 3))
pd.DataFrame(rows).to_csv(f"{OUT}/sobol_output_diagnostics.csv", index=False)
json.dump(res, open(f"{OUT}/sobol_closing_criterion.json", "w"), indent=1)
for k, v in res.items():
    print(k, {kk: vv for kk, vv in v.items() if kk != "top10"})
    for r in v["top10"]: print("   ", r)
