"""Fase 1 analisis independiente del Sobol: (1) salidas degeneradas, (2) S_T medio robusto, (3) criterio de cierre
CI(S_T) < 0.05 en el top 10. Bootstrap propio (Jansen) sobre las N muestras base, mismo remuestreo para todas las salidas."""
import numpy as np, pandas as pd, json, sys, os
sys.path.insert(0, "src"); import sens_v6 as S
OUT = "results/analysis_phase1"; D = len(S.NAMES); names = S.out_names()
rng = np.random.default_rng(2026); B = 500
res = {}; rows = []
for t in (5, 3):
    Y = np.load(f"sens_out_v6_{t}s/Y_outputs.npy"); N = Y.shape[0] // (D + 2)
    Yr = Y.reshape(N, D + 2, -1); A = Yr[:, 0]; Bm = Yr[:, -1]; AB = Yr[:, 1:-1]   # (N,out),(N,out),(N,D,out)
    def st(idx):
        a, b, ab = A[idx], Bm[idx], AB[idx]
        V = np.var(np.r_[a, b], axis=0)
        num = 0.5 * np.mean((a[:, None, :] - ab) ** 2, axis=0)            # (D,out)
        s1 = np.mean(b[:, None, :] * (ab - a[:, None, :]), axis=0)
        with np.errstate(all="ignore"):
            return num / V, s1 / V
    ST0, S10 = st(np.arange(N))
    boots = np.stack([st(rng.integers(0, N, N))[0] for _ in range(B)])   # (B,D,out)
    # diagnostico de degeneracion: fraccion de evaluaciones con salida != valor base modal (UHI>0)
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
        # criterio por salida individual: half-width SALib-like (1.96*sd boot) en el top 10 de cada salida
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
