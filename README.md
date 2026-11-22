[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24358421.svg)](https://doi.org/10.5281/zenodo.24358421)

# UHI v6.1: facility-specific, probabilistic priority lists of temperate exoplanets

**Author:** Dr. Effy Viridis, Independent Researcher, ORCID [0009-0004-5620-5352](https://orcid.org/0009-0004-5620-5352)
**Release:** `v6.1.1` (2026-11-21), version DOI [10.5281/zenodo.24431744](https://doi.org/10.5281/zenodo.24431744), a documentation release of v6.1: English code comments, README corrections and the top-N membership analysis. Model code and numerical outputs are identical to `v6.1` (version DOI [10.5281/zenodo.24430693](https://doi.org/10.5281/zenodo.24430693)). · **Code:** MIT · **Data, results, figures:** CC BY 4.0

## In one paragraph
This package ranks planets from the NASA Exoplanet Archive by how worthwhile they are as targets for atmospheric characterization. Every planet gets 10,000 Monte Carlo realizations that carry the catalogue uncertainties through four layers: climate (with a smooth tidal-locking transition), atmospheric retention under the star's cumulative XUV output, a per-draw composition class (rocky, Hycean or gas-rich), and a detection model for the one facility that suits it (JWST for transiting planets, ELT/PCS reflected light for the rest). For each planet the code reports the probability of passing every veto (`P_alive`), the score conditional on passing, the probability of landing in the top 10 (`P_top10`), and the composite `UHI_rank`. Lists are produced **per facility and per composition class**. They are not meant to be merged into one cross-facility ranking.

**Status of the branches.** The rocky branch is the main product. Its global sensitivity analysis meets the closing criterion on the half-width reading (10/10 leading parameters at both thresholds) but not on the stricter full-width reading (9/10 at both thresholds). The Hycean branch is released as **INDICATIVE (not converged)**; see "Global sensitivity" below.

The previous model, v5.2, stays frozen and citable on Zenodo (version DOI [10.5281/zenodo.24358996](https://doi.org/10.5281/zenodo.24358996), Git tags `v5.2-final`, `v5.2.1`, `v5.2.2`). v6.1 is a different model, so its numbers should not be compared one to one with v5.2 tables.

## Concrete use cases
1. **Choosing JWST transmission targets among temperate rocky planets.** Open `results/v6.1/uhi_v6_rocky_main_5s.csv`, keep `transiting == 1`, and sort by `UHI_rank`. At 5σ the first five are TRAPPIST-1 e, Kepler-442 b, TRAPPIST-1 f, TRAPPIST-1 d and TOI-1452 b. Each row also gives `P_alive`, `P_top10`, `P_R` (probability of being rocky) and the median sub-scores, so a proposer can see *why* a planet ranks where it does.
2. **Building an ELT reflected-light shortlist.** Same file, `transiting == 0`. At 5σ the first five are GJ 273 b, Wolf 1061 c, GJ 1061 d, Proxima Cen b and GJ 3378 b.
3. **Deciding which measurement or theory input to improve.** `tables/sobol_v6.1_n512_rocky_5s.csv` lists the assumptions that move the rocky lists the most. The position and the width of the cosmic shoreline lead by a clear margin, so better constraints on XUV-driven atmospheric escape would change target priorities more than any other single input we sampled.
4. **Rerunning with your own assumptions.** Change a threshold or a prior on the command line (`run_uhi_v6.py --thr 3 --variant main`) and regenerate the lists in a few minutes (see "Reproduce").

```python
import pandas as pd
cols = ["planet", "UHI_rank", "P_alive", "P_top10", "P_R"]
for thr in ("5s", "3s"):                                   # 5σ (main) and 3σ detection thresholds
    d = pd.read_csv(f"results/v6.1/uhi_v6_rocky_main_{thr}.csv")
    jwst = d[d.transiting == 1].sort_values("UHI_rank", ascending=False)   # JWST: transiting planets
    elt = d[d.transiting == 0].sort_values("UHI_rank", ascending=False)    # ELT/PCS: non-transiting planets
    print(thr, "JWST\n", jwst[cols].head(10), "\n", thr, "ELT/PCS\n", elt[cols].head(10))
```
Replace `rocky` by `hycean` for the (indicative) Hycean lists. `P_top10` in these files is computed over the whole branch list, mixing both facilities; for membership within one facility list use the top-N tables below.

## Main results (10,000 draws per planet, 6,112 planets scored)
| Branch | Alive (P_alive > 0), 5σ and 3σ | Alive with P_alive ≥ 0.01 (5σ) | Of the alive: JWST / ELT |
|---|---|---|---|
| Rocky | 255 | 130 | 167 / 88 |
| Hycean (indicative) | 71 | 22 | 26 / 45 |
| Combined | 287 | 148 | |

Top five per facility, rocky branch (`UHI_rank` in parentheses):

| Threshold | JWST (transiting) | ELT/PCS (non-transiting) |
|---|---|---|
| 5σ | TRAPPIST-1 e (0.082), Kepler-442 b (0.067), TRAPPIST-1 f (0.061), TRAPPIST-1 d (0.054), TOI-1452 b (0.045) | GJ 273 b (0.287), Wolf 1061 c (0.249), GJ 1061 d (0.222), Proxima Cen b (0.202), GJ 3378 b (0.164) |
| 3σ | Kepler-442 b (0.230), TRAPPIST-1 e (0.145), Kepler-452 b (0.139), Kepler-1544 b (0.135), Gliese 12 b (0.125) | GJ 273 b (0.287), Wolf 1061 c (0.249), GJ 1061 d (0.226), Proxima Cen b (0.202), GJ 3378 b (0.199) |

The Hycean lists are in `results/v6.1/uhi_v6_hycean_main_{5s,3s}.csv`. Treat them as indicative. The prior-robustness tests (22 class-prior and mass-radius variants) are in `tables/v6.1/prior_sens_*` and `results/analysis_phase1/v61_prior_robustness.csv`.

## Global sensitivity (Saltelli–Sobol, N = 512)
32 sampled assumptions, 17,408 model evaluations per detection threshold, 100 draws for each of 1,214 planets. Total-order indices S_T come with a 95% bootstrap interval (B = 500) on their mean over the non-degenerate outputs. **Closing criterion:** the interval for each of the ten leading parameters must be narrower than 0.05. Both readings are reported. Phase 1 adopts the half-width reading; the full-width counts are given for transparency.

| Branch | Outputs used / excluded | Half-width < 0.05 (5σ, 3σ) | Full width < 0.05 (5σ, 3σ) | Leading S_T (5σ; 3σ) | Verdict |
|---|---|---|---|---|---|
| Rocky | 15 / 2 | 10/10, 10/10 | 9/10, 9/10 | shoreline offset `shore_logc` 0.22; 0.25 · shoreline width `shore_sig` 0.12; 0.14 · mass-route `mv_c` 0.11; 0.13 | Half-width criterion met; full-width criterion not met (9/10) |
| Combined | 16 / 1 | 10/10, 10/10 | 7/10, 8/10 | `shore_logc` 0.20; 0.23 · `lph2_c` 0.13; 0.11 · `shore_sig` 0.11; 0.13 | Half-width criterion met; full-width criterion not met (7/10, 8/10) |
| Hycean | 6 / 11 | 6/10, 6/10 | 2/10, 3/10 | H2 pressure `lph2_c` 0.71; 0.71 · Innes tolerance `fcl_c` 0.38; 0.37 · Hycean albedo `hyc_alb` 0.29; 0.30 | **Neither criterion met: indicative** |

Why the Hycean branch is indicative: 11 of its 17 outputs are degenerate (the score is zero in almost every evaluation) and are excluded, so the branch estimate rests on six outputs. A larger N cannot recover those 11 outputs, because their degeneracy comes from the model (almost no Hycean survivors), not from sampling. Whether a larger N would narrow the intervals of the six retained outputs has not been tested. Files: `tables/sobol_v6.1_n512_closing.json`, `tables/sobol_v6.1_n512_{rocky,hycean,combined}_{5s,3s}.csv`, `tables/sobol_v6.1_n512_outputs.csv`, figure `figures/fig_v6.1_n512_sobol_ST.png`.

## Reproduce
**Prerequisites.** Linux with Python 3.13 available as `python3.13` (tested with 3.13.5; the pinned versions in `requirements.txt` target 3.13) and Git. Run every command from the repository root. Nothing below needs network access after `pip install`.

The commands write to a separate folder, `repro/` (ignored by Git), so the shipped results stay untouched and can be compared byte for byte.
```bash
python3.13 -m venv .venv
.venv/bin/pip install -r requirements.txt        # exact versions used for v6.1
.venv/bin/python tests/test_v61.py                # 21 unit and regression checks (seconds)
.venv/bin/python tests/test_v6.py                 # v5.2 compatibility on the full catalogue (~4 min)

# Main lists at 5σ and 3σ (~3 min each on 1 core); --outdir is relative to the repository root
.venv/bin/python run_uhi_v6.py --variant main --thr 5 --outdir repro/results
.venv/bin/python run_uhi_v6.py --variant main --thr 3 --outdir repro/results
for f in repro/results/*.csv; do cmp "$f" "results/v6.1/$(basename "$f")" && echo "identical $(basename "$f")"; done
# summary_*.json differ only in the "seconds" runtime field

# Sobol report from the shipped Sobol outputs; a new tag writes new files instead of overwriting v6.1_n512
.venv/bin/python src/sobol_report_v6.py sens_out_v6_5s_n512 sens_out_v6_3s_n512 repro_n512
for b in rocky hycean combined; do for t in 5s 3s; do cmp tables/sobol_repro_n512_${b}_${t}.csv tables/sobol_v6.1_n512_${b}_${t}.csv && echo "identical ${b}_${t}"; done; done
cmp tables/sobol_repro_n512_closing.json tables/sobol_v6.1_n512_closing.json && echo "identical closing"
# the figure differs only because its title contains the tag; the repro_n512 files are ignored by Git and can be deleted

# Prior-robustness tables from the shipped variant results (the three variables are required)
UHI_RES=results/v6.1 UHI_TAB=repro/tables/v6.1 UHI_FIG=repro/figures/v6.1 .venv/bin/python src/prior_sensitivity_v6.py

# Top-N membership within each facility list (see "Top-N membership" below; ~3 min per threshold)
.venv/bin/python analysis/topN_membership_v61.py --thr 5 --scratch repro/topN_scratch
.venv/bin/python analysis/topN_membership_v61.py --thr 3 --scratch repro/topN_scratch
.venv/bin/python analysis/topN_membership_v61.py --summarize --scratch repro/topN_scratch --tables repro/tables
for f in counts.csv targets.csv summary.json; do cmp repro/tables/v6.1_topN_membership_$f tables/v6.1_topN_membership_$f && echo "identical $f"; done
```
Without `UHI_RES`, `src/prior_sensitivity_v6.py` looks for results in `results/`, which holds no CSVs in the release, and stops with `FileNotFoundError`. `UHI_TAB` and `UHI_FIG` only choose where tables and figures are written (defaults `tables/` and `figures/`; the shipped versions are in `tables/v6.1/` and `figures/v6.1/`).

**Full batch scripts.** These write into the shipped locations; run them on a separate copy of the repository if you want to compare against the release.
- `bash run_v61_mc.sh`: all 24 variant runs plus the legacy regression (8 processes), into `results/v6.1` and `results/v6.1_legacy_check`.
- `bash run_v61_sobol.sh`: full Sobol N=512 from scratch into `sens_out_v6_{5s,3s}_n512` (on our machine, roughly a day per threshold with 4 workers), then the report with tag `v6.1_n512`.
- `bash pipeline_v61.sh`: chains the MC batch, prior-robustness tables (`tables/v6.1`, `figures/v6.1`), `tests/test_v6.py` and the Sobol run.

The scripts expect the virtual environment at `./.venv` and write their logs to `logs/v6.1/` (one file per run, e.g. `mc_main_5s.log`, `sobol_5s_n512.log`, `pipeline.log`). That folder ships with the release logs and the frozen environment listing (`env_resume_venv_freeze.txt`); reruns overwrite the log files of the same name. Sobol samples use SALib with seed 11; the Monte Carlo uses fixed seeds, so reruns in the same environment are bit-identical.

**Checking that the commented code equals v6.1.** v6.1.1 changes only comments and docstrings in the Python sources (plus two development-machine paths in the scripts). To verify:
```bash
mkdir -p repro/v61 && git archive v6.1 src run_uhi_v6.py | tar -x -C repro/v61
.venv/bin/python tools/ast_identity_check.py repro/v61 . $(cd repro/v61 && ls src/*.py src/analysis/*.py run_uhi_v6.py)
```
It compares the syntax trees of both versions with docstrings removed and prints `ALL_IDENTICAL` when no code differs.

## Top-N membership (added in v6.1.1)
How stable is a published facility list? `analysis/topN_membership_v61.py` reruns the released `main` configuration in memory (seed 42, 10,000 draws), checks that it reproduces `results/v6.1` (UHI_rank, P_alive, P_top10 and rank_med agree to within 1e-12), and then ranks each facility list (JWST = transiting, ELT/PCS = non-transiting) separately in every draw. P(top-N) is the fraction of draws in which a planet is among the first N of its list with a positive score, for N = 5 and 10; P(top-N | alive) divides by P_alive.

| Rocky branch, nominal top-N | P(top-N) ≥ 0.5 (≥ 0.8) | P(top-N \| alive) ≥ 0.5 (≥ 0.8) |
|---|---|---|
| 5σ JWST top 10 | 3 (0) | 10 (2) |
| 5σ ELT/PCS top 10 | 4 (0) | 10 (9) |
| 3σ JWST top 10 | 1 (0) | 9 (2) |
| 3σ ELT/PCS top 10 | 3 (0) | 10 (7) |
| Top 5, any list | 0 (0) | 1–5 (0–2) |

Membership is limited mainly by survival: in the draws where a planet survives, the published top 10 is largely stable, while positions within the top 5 swap often. For Hycean lists P(top-N) equals P_alive (≤ 0.15), so their order carries no information beyond survival. Files: `tables/v6.1_topN_membership_{targets,counts}.csv`, `tables/v6.1_topN_membership_summary.{json,md}`.

## Legacy scripts (kept as records, not rerunnable from the release)
- `src/make_outputs_v6.py` reads `results/` and needs the `teff`, `bt3d`, `jackson`, `only_clim`, `only_xuv`, `only_d3` and `single` variant runs, which are not shipped; its tables were produced during development.
- `src/analysis/a1_*` to `a7_*` are the Phase 1 audit scripts of the v6.0-to-v6.1 review. They use v6.0-era paths (`sens_out_v6_{5,3}s`, `results/uhi_v6_*.csv`); `a7` needs v6.0 results that are not shipped, and `a4` hard-codes the preprint Innes value (0.0852) instead of the published 0.0860. Their outputs are in `results/analysis_phase1/`.
- `src/build_mors_grid.py` needs the third-party Mors code and stellar tracks (`ext/`, not included); the grid it produces ships as `data/mors_grid.npz`.

## Contents
| Path | What |
|---|---|
| `run_uhi_v6.py` | Runs one variant at one threshold and writes the per-branch CSVs and a summary JSON |
| `analysis/topN_membership_v61.py` | Top-N membership within each facility list (v6.1.1) |
| `tools/ast_identity_check.py` | Verifies that two versions of the sources differ only in comments and docstrings |
| `src/uhi_code_v6.py`, `src/xuv_history.py`, `src/composition.py` | The v6 model: climate transition, cumulative XUV and retention, composition classes |
| `src/uhi_code.py`, `src/anchors_v5.py` | v5.2 core and synthetic anchors, reused through the compatibility switch |
| `src/sens_v6.py`, `src/sobol_report_v6.py` | Sobol sampling (resumable) and the report with degenerate-output filtering and bootstrap intervals |
| `src/prior_sensitivity_v6.py`, `src/fit_empirical_prior.py` | Prior-robustness tables and the empirical mass-route prior fit |
| `src/make_outputs_v6.py`, `src/build_mors_grid.py`, `src/analysis/` | Legacy and build scripts (see "Legacy scripts") |
| `tests/` | `test_v61.py` (fixes and regression) and `test_v6.py` (v5.2 compatibility, manual ARM values, anchors) |
| `backups/pre_v6.1_20261120_1625/` | The v6.0 source, used by `test_v61.py` as the regression reference |
| `data/` | Catalogue snapshot `pscomppars_2026-11-18.csv`, v5.2 reference results, Zeng et al. (2019) mass-radius curves, the precomputed Mors grid `mors_grid.npz`, E. Mamajek's mean dwarf colour and T_eff sequence |
| `results/v6.1/` | Per-branch results for every variant and both thresholds, plus summaries and anchors |
| `results/v6.1_legacy_check/` | Regression run with the v6.1 fixes switched off |
| `results/analysis_phase1/` | Independent checks: hand recomputation, anchors, Innes vs Madhusudhan, prior robustness, Sobol diagnostics |
| `tables/`, `figures/` | v6.1 tables and figures (Sobol N=512, top-N membership, prior sensitivity, Innes vs Madhusudhan) |
| `sens_out_v6_{5s,3s}_n512/` | Sobol outputs needed to regenerate the report (`Y_outputs.npy`, `base_raw.npy`, `meta.json`, `sobol_results.csv`) |
| `logs/v6.1/` | Run logs written by the batch scripts and the frozen environment listing; `logs/release_v6.1/` has the pre-release test logs |
| `docs/assumptions_v6_es.md` | Every assumption with its source or flagged as an assumption (Spanish) |
| `docs/AI_DISCLOSURE.md` | How AI agents were used and how their output was checked |
| `patches/` | Small Python 3.13 fix applied to the Mors code when building the grid |

**Not included:** the raw Sobol arrays `Y_raw.npy` (about 250 MB per threshold, above the GitHub file limit; only needed to resume an interrupted run), and the third-party Mors code and stellar-track grid (`ext/`), which are only needed to rebuild `data/mors_grid.npz`. Get Mors from Johnstone et al. (2021) and its stellar models from the Zenodo record cited in `src/build_mors_grid.py`.

## Known limitations
- The Hycean branch is indicative (see above). Its lists also depend strongly on the climate criterion: Innes et al. (2023) and Madhusudhan et al. (2021) give 71 vs 1,021 alive planets (`tables/v6.1/innes_vs_madhu.md`).
- The tidal-locking rule is a smoothed assumption; a Gladman-type locking time is planned for the next phase.
- `docs/assumptions_v6_es.md` is in Spanish, and so are some printed messages, figure labels and a few output column names (e.g. `evidencia`). Code comments and docstrings are in English.

## Licenses and data sources
- Code (`src/`, `tests/`, scripts): [MIT](LICENSE). Data products, results, tables and figures: [CC BY 4.0](LICENSE-DATA).
- NASA Exoplanet Archive, Planetary Systems Composite Parameters (`pscomppars`), snapshot of 2026-11-18. This research has made use of the NASA Exoplanet Archive, operated by Caltech under contract with NASA under the Exoplanet Exploration Program.

## Cite
Use `CITATION.cff`. Concept DOI (always resolves to the newest version): [10.5281/zenodo.24358421](https://doi.org/10.5281/zenodo.24358421). For exact reproducibility, cite the version DOI of the release you used: v6.1 is [10.5281/zenodo.24430693](https://doi.org/10.5281/zenodo.24430693); v6.1.1, with identical numerical outputs, is [10.5281/zenodo.24431744](https://doi.org/10.5281/zenodo.24431744).
