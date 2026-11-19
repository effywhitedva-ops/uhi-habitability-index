[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24358422.svg)](https://doi.org/10.5281/zenodo.24358422)

# UHI v5.2-final: Unified Habitability Index, reproducibility package

**Lead author:** Effy Viridis · **Release:** `v5.2-final` (2026-11-18) · **Code license:** MIT · **Data/results/figures:** CC-BY 4.0

UHI is a probabilistic habitability-and-observability index: a Dirichlet-weighted geometric mean of 7 sub-indices (S_fis, S_clima, S_atm, S_estrella, S_geo, S_mag, S_obs). It is evaluated by Monte Carlo (10 000 draws per planet) over the NASA Exoplanet Archive `pscomppars` table. Planets are ranked by UHI_rank = P(UHI > 0) × median(UHI | UHI > 0).

This is the **final** frozen version of the model. Repository: https://github.com/effywhitedva-ops/uhi-habitability-index. Zenodo DOI for `v5.2-final`: [10.5281/zenodo.24358422](https://doi.org/10.5281/zenodo.24358422) (all versions: [10.5281/zenodo.24358421](https://doi.org/10.5281/zenodo.24358421)).

## Licenses
- **Code** (`src/`, `notebooks/`, `run_uhi_v52.py` and the textual documentation): [MIT License](LICENSE), copyright Effy Viridis, 2026.
- **Data, results and figures** (`data/`, `results/`, `tables/`, `figures/`): [Creative Commons Attribution 4.0 International (CC BY 4.0)](LICENSE-DATA).
- **PHL data are not redistributed.** The Habitable Worlds Catalog is downloaded at run time from its public URL, under PHL's own terms. Check those terms before publishing any PHL-derived columns. The NASA Exoplanet Archive snapshot in `data/` is public NASA/IPAC data; follow the Archive's acknowledgment policy.

See also `CITATION.cff`.

## Contents
| Path | What |
|---|---|
| `run_uhi_v52.py` | Master script: runs everything with no manual steps |
| `src/uhi_code.py` | Model (v5.2-final) |
| `src/anchors_v5.py` | Synthetic Earth / Venus / Mars at 10 pc (controls) |
| `src/sweep_v52.py` | 1-D sweep of the runaway-greenhouse veto (F_rg, T_rg) |
| `src/make_tables.py`, `src/make_figures.py`, `src/make_mrt.py` | All tables/statistics, figures (PDF + PNG) and the AAS MRT table |
| `src/uhi_sens.py` | Sobol/LHS sensitivity model (post-v5; not re-run for v5.1 or v5.2) |
| `data/pscomppars_2026-11-18.csv` / `.parquet` | Dated catalog snapshot (6470 rows; md5 of CSV `f2f70b24e67a9054940267e8081c81e8`) |
| `results/uhi_results_v52.csv` / `.mrt` | Full results, N = 6112 scored planets (5σ main) |
| `results/uhi_results_v52_3sigma.csv` | 3σ variant |
| `results/anchors_v52*.csv`, `results/sweep_v52/` | Controls and sweep runs |
| `tables/` | `tables_v52.md`, `stats_v52.json` and CSV tables (generated) |
| `figures/` | (a) Sobol S_i/S_Ti post-v5, (b) P_alive vs conditional median, (c) UHI 5σ/3σ vs ESI/PHL; plus anchors and top-20 previews |
| `sens_out_v5/` | Sobol post-v5 summary outputs (inputs of figure a) |
| `notebooks/UHI_v5.2_colab.ipynb` | Colab notebook (runs the master script) |
| `docs/` | Assumptions (Spanish), the 5σ vs 3σ discussion (Spanish) and the Software/Acknowledgments note (English) |
| `LICENSE`, `LICENSE-DATA`, `CITATION.cff` | Licenses and citation metadata |

## Reproduce
```bash
pip install -r requirements.txt
python run_uhi_v52.py                # ~7 min on 2 CPUs without the sweep; the sweep adds 14 full runs
python run_uhi_v52.py --skip-sweep   # main runs, anchors, benchmark, tables and figures only
```
Outputs go to `output/`. The script checks the snapshot md5 before running. The certified run used Google Colab (CPU). A local re-run with newer library versions gave the same results to within 1.9e-9.

## Certified v5.2-final results (5σ main; 3σ variant)
- 6112 planets scored; 515 alive (P_alive > 0) in both variants; Spearman 5σ vs 3σ (live) = 0.985.
- Anchors at 10 pc: Earth at the 99.30 / 99.26 percentile (non-transiting / transiting) at 5σ and 99.85 / 99.79 at 3σ. Venus and Mars have P_alive = 0. **All pass.**
- Top 10 (5σ): TOI-4336 A b, TOI-7166 b, GJ 273 b, K2-332 b, TRAPPIST-1 d, Wolf 1061 c, Proxima Cen d, TOI-2094 b, GJ 1061 d, K2-125 b.

## Accepted limitation (not tuned)
The 1-D veto sweep passes for F_rg = 275–295 W m⁻² (including Goldblatt et al. 2013 at 282 W m⁻² and the 2014 corrigendum absorbed-solar value of 286 W m⁻²) and for T_rg = 3850–3900 K. It **fails the “top-10 composition invariant” criterion at T_rg = 3700 and 3750 K**: K2-125 b (T_eff = 3654 K) leaves #10 and Ross 128 b enters. At the default T_rg = 3800 K the margin is 0.333 vs 0.323 in UHI_rank; the other 9 members of the top 10 are invariant. Venus remains at P_alive = 0 at every primary sweep value, and Spearman stays above 0.99. This boundary sensitivity is documented as an accepted limitation. The default T_rg was **not** adjusted to make the test pass.

## Notes
- **ESI erratum.** Versions v1–v5.1 used a total exponent of w/8 instead of w/4, which compressed ESI toward 1. The formula is corrected in v5.2 (Mars 0.70 and Venus 0.44 with the standard reference properties, as in Schulze-Makuch et al. 2011). This does **not** affect the UHI. ESI values quoted in the earlier v4 and v5 documents are inflated.
- **Sobol coverage.** The Sobol analysis is the post-v5 run (51 assumptions). It does **not** include F_rg (added in v5.1) or the Yang albedo ramp w_eff (added in v5.2). It was not re-run, because the top 10 did not change from v5 to v5.2.

## External data
- NASA Exoplanet Archive, Planetary Systems Composite Parameters (`pscomppars`), TAP query of 2026-11-18 (exact query in `src/fetch_catalog.py`).
- PHL @ UPR Arecibo, Habitable Worlds Catalog (`hwc.csv`, last updated 2024-03-12, md5 `e971b8dcc305a2fb142d2840aaa8ca37`). It is **downloaded at run time and not redistributed**. `tables/uhi_results_v52_with_benchmark.csv` and the benchmark tables contain the PHL ESI and habitability flag for matched planets; check PHL's terms before publishing them.

## Pending (outside this package)
- Cite with `CITATION.cff` and the Zenodo DOI of the `v5.2-final` release.
