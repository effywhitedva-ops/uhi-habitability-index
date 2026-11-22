#!/usr/bin/env python3
"""Top-N membership robustness for UHI v6.1 (standalone analysis; does not modify UHI code or outputs).

Why a re-run: the released CSVs keep only per-planet summaries (P_top10 over the whole branch list, mixing
facilities); the per-draw score matrices are discarded inside uhi_code_v6.finalize. This script re-runs the
released v6.1 'main' configuration (seed 42, 10,000 draws) in memory and captures those matrices by wrapping
finalize at run time (the source file is not edited). Validation: the CSVs the wrapped run produces must match
results/v6.1/uhi_v6_{branch}_main_{thr}s.csv exactly (UHI_rank, P_alive, P_top10, rank_med), proving the draws
are the released ones.

Definition: per branch (rocky, hycean), per facility list (JWST = transiting == 1, ELT/PCS = transiting == 0),
per threshold: in each draw, rank that facility's planets by the draw's UHI score (dead draw = 0, never counted);
P(top-N) = fraction of draws in which the planet is in the top N with a positive score. N = 5 and 10.
Nominal top-N = the published order (UHI_rank, as in README); the median-per-draw-rank order is reported too.
Supplementary: P(top-N | alive) = P(top-N) / P_alive, because a planet cannot be in a list in draws where it is vetoed.

Usage: run from the repository root (about 3 min per threshold on one core):
  .venv/bin/python analysis/topN_membership_v61.py --thr 5 --scratch repro/topN_scratch
  .venv/bin/python analysis/topN_membership_v61.py --thr 3 --scratch repro/topN_scratch
  .venv/bin/python analysis/topN_membership_v61.py --summarize --scratch repro/topN_scratch --tables repro/tables
Defaults: --scratch analysis/_topN_scratch (per-draw matrices, not shipped) and --tables tables (the shipped
tables/v6.1_topN_membership_*); pass a separate --tables directory to compare instead of overwriting.
Relative paths are taken from the repository root.
"""
import argparse, os, sys, json, hashlib, warnings
warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(ROOT, "src"))
import numpy as np, pandas as pd

ap = argparse.ArgumentParser(); ap.add_argument("--thr", type=int); ap.add_argument("--summarize", action="store_true")
ap.add_argument("--scratch", default=os.path.join("analysis", "_topN_scratch")); ap.add_argument("--tables", default="tables")
a = ap.parse_args()
TAB = os.path.join(ROOT, a.tables); SCR = os.path.join(ROOT, a.scratch)   # os.path.join keeps absolute paths as given
os.makedirs(SCR, exist_ok=True); os.makedirs(TAB, exist_ok=True)
BR = ["rocky", "hycean"]; FAC = {"JWST": 1, "ELT": 0}; NS = [5, 10]

if not a.summarize:
    import uhi_code as u, uhi_code_v6 as v6
    SNAP = os.path.join(ROOT, "data", "pscomppars_2026-11-18.csv")
    assert hashlib.md5(open(SNAP, "rb").read()).hexdigest() == "f2f70b24e67a9054940267e8081c81e8"
    u.V3["snr_thr"] = float(a.thr)
    P = {**u.V3, **v6.V6}; cfg = dict(v6.CFG_V6)             # variant 'main' as in run_uhi_v6.py
    captured = []; _orig = v6.finalize
    def finalize_capture(rows, draws, joint=True, strict=True):
        captured.append(({i: np.asarray(x, dtype=np.float64) for i, x in draws.items()}, [r["planet"] for r in rows]))
        return _orig(rows, draws, joint, strict)
    v6.finalize = finalize_capture
    res = v6.run_v6(pd.read_csv(SNAP), 10000, seed=42, cfg=cfg, P=P)
    names = list(res.keys())                                   # same order as the finalize calls
    rows_out = []; check = {}
    for b, (draws, planets) in zip(names, captured):
        if b not in BR: continue
        rel = pd.read_csv(os.path.join(ROOT, "results", "v6.1", f"uhi_v6_{b}_main_{a.thr}s.csv"))
        new = res[b]
        m = rel.merge(new, on="planet", suffixes=("_rel", "_new"))
        assert len(m) == len(rel) == len(new)
        check[b] = {c: float(np.nanmax(np.abs(m[c + "_rel"] - m[c + "_new"]))) for c in ["UHI_rank", "P_alive", "P_top10", "rank_med"]}
        assert all(v <= 1e-12 for v in check[b].values()), (b, check[b])
        idx = np.array(sorted(draws)); M = np.vstack([draws[i] for i in idx]); pl = np.array(planets)[idx]
        tr = rel.set_index("planet").loc[pl, "transiting"].values
        np.savez_compressed(os.path.join(SCR, f"draws_{b}_{a.thr}s.npz"), M=M.astype(np.float32), planet=pl, transiting=tr)
        for fac, fv in FAC.items():
            sel = tr == fv; Mf = M[sel]; pf = pl[sel]
            order = np.argsort(-Mf, axis=0, kind="stable"); rk = np.empty_like(order)
            np.put_along_axis(rk, order, np.arange(1, Mf.shape[0] + 1)[:, None].repeat(Mf.shape[1], 1), axis=0)
            alive = Mf > 0; npos = alive.sum(0); rk_eff = np.where(alive, rk, npos[None, :] + 1)
            sub = rel.set_index("planet").loc[pf]
            df = pd.DataFrame({"branch": b, "thr": f"{a.thr}s", "facility": fac, "planet": pf,
                               "UHI_rank": sub.UHI_rank.values, "P_alive": sub.P_alive.values,
                               "rank_med_fac": np.median(rk_eff, 1), "P_top5": (alive & (rk <= 5)).mean(1),
                               "P_top10": (alive & (rk <= 10)).mean(1)})
            # facility planets with no stored draws (P_alive = 0) never enter a list
            df = df.sort_values(["UHI_rank"], ascending=False).reset_index(drop=True)
            df["nominal_pos"] = np.arange(1, len(df) + 1)
            df["medrank_pos"] = df.rank_med_fac.rank(method="first").astype(int)
            rows_out.append(df)
    out = pd.concat(rows_out)
    out.to_csv(os.path.join(SCR, f"topN_{a.thr}s.csv"), index=False)
    json.dump(check, open(os.path.join(SCR, f"validation_{a.thr}s.json"), "w"), indent=1)
    print("validation max|diff| vs released CSVs:", check)
else:
    out = pd.concat([pd.read_csv(os.path.join(SCR, f"topN_{t}s.csv")) for t in (5, 3)])
    val = {f"{t}s": json.load(open(os.path.join(SCR, f"validation_{t}s.json"))) for t in (5, 3)}
    for N in NS:   # supplementary: membership conditional on the planet surviving the draw (P_topN is capped by P_alive)
        out[f"P_top{N}_given_alive"] = np.where(out.P_alive > 0, out[f"P_top{N}"] / out.P_alive.where(out.P_alive > 0, 1), 0.0)
    keep = out[(out.nominal_pos <= 10) | (out.medrank_pos <= 10)].sort_values(["branch", "thr", "facility", "nominal_pos"])
    keep.to_csv(os.path.join(TAB, "v6.1_topN_membership_targets.csv"), index=False)
    summ = {"definition": __doc__.split("Definition:")[1].split("Usage:")[0].strip(), "validation_max_abs_diff": val, "lists": []}
    for (b, t, f), g in out.groupby(["branch", "thr", "facility"], sort=False):
        for N in NS:
            pc = f"P_top{N}"
            for basis, pos in [("nominal_UHI_rank", "nominal_pos"), ("median_draw_rank", "medrank_pos")]:
                top = g[g[pos] <= N]
                summ["lists"].append(dict(branch=b, thr=t, facility=f, N=N, basis=basis, n_in_list=int(len(top)),
                    n_ge_0p5=int((top[pc] >= 0.5).sum()), n_ge_0p8=int((top[pc] >= 0.8).sum()),
                    n_cond_ge_0p5=int((top[pc + "_given_alive"] >= 0.5).sum()), n_cond_ge_0p8=int((top[pc + "_given_alive"] >= 0.8).sum()),
                    min_P=round(float(top[pc].min()), 4), members=[f"{p} ({v:.3f})" for p, v in zip(top.planet, top[pc])]))
    json.dump(summ, open(os.path.join(TAB, "v6.1_topN_membership_summary.json"), "w"), indent=1)
    pd.DataFrame(summ["lists"]).drop(columns="members").to_csv(os.path.join(TAB, "v6.1_topN_membership_counts.csv"), index=False)
    print("written")
