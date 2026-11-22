# Changelog

## v6.1 (2026-11-21)
First public release of the v6 model line. Archived on Zenodo under the same concept DOI as v5.2 ([10.5281/zenodo.24358421](https://doi.org/10.5281/zenodo.24358421)). The v5.2 package remains available from its Git tags and Zenodo versions.

**Decision record.** On 2026-11-21 the author approved this release after the Sobol N=512 run finished, with the Hycean branch labelled indicative (not converged). Phase 1 of the v6 plan is closed for the rocky branch, which leads the analysis. Nothing in the model or its outputs was changed for the release; this step only documents and archives.

### Model changes relative to v5.2
- Climate: the synchronous/non-synchronous regime is drawn per realization from a smooth tidal-locking probability and an orbital-period-dependent efficiency, replacing the hard T_eff cut.
- Atmospheric retention: the atmosphere sub-score now comes from cumulative XUV exposure against a cosmic shoreline (Pass et al. 2025 for fully convective stars, Johnstone et al. 2021 Mors tracks above 0.35 M_sun, Jackson et al. 2012 as an alternative).
- Composition: every draw is classed as rocky, Hycean or gas-rich using Zeng et al. (2019), Otegi et al. (2020), Venturini et al. (2024) and Madhusudhan et al. (2021). Separate rocky, Hycean and combined lists are produced. Unmeasured masses come from class-specific Otegi relations.
- `P_top10` counts only surviving draws.
- A compatibility switch reproduces v5.2 on the full catalogue (maximum difference 9×10⁻¹⁵, which the v5.2 code itself shows when rerun in the same numerical environment).

### Fixes adopted in v6.1 (after the Phase 1 audit, 2026-11-20)
- Physical floor on stellar radius, R* ≥ 0.5 R_MS(M*) (white dwarfs exempt). It removes unphysical tails that inflated the number of surviving planets.
- The class prior now acts in the radius and mass routes as well; with a flat prior the output is bit-identical to v6.0.
- Innes et al. (2023) critical instellation for M stars at 10 bar set to 0.0860.
- Sobol report rewritten: degenerate outputs are excluded, the mean S_T gets a bootstrap 95% interval (B = 500), and both full width and half-width are reported. Validated on the Ishigami function. The v6.0 Hycean and combined Sobol bars were invalid for this reason and are not shipped.
- New flags: `run_uhi_v6.py --outdir` and `--legacy` (switches the fixes off for regression).

### Results
- Alive (P_alive > 0) at 5σ and 3σ: rocky 255, Hycean 71, combined 287 (v6.0: 495, 242, 530).
- Sobol N=512 (32 parameters, 17,408 evaluations per threshold): rocky converged (half-width 10/10 at both thresholds, full width 9/10); combined converged on the half-width reading (10/10; full width 7/10 and 8/10); Hycean not converged (half-width 6/10, full width 2/10 and 3/10; 11 of 17 outputs degenerate). Rocky is driven by the shoreline offset and width (S_T 0.22 and 0.12 at 5σ, 0.25 and 0.14 at 3σ), then the mass-route parameter `mv_c`. Hycean is driven by the H2 pressure `lph2_c` (S_T 0.71), then `fcl_c` and `hyc_alb`.

### Verification before tagging
- `tests/test_v61.py`: 21/21 pass. `tests/test_v6.py` (v5.2 compatibility, manual ARM values, anchors): pass. Log: `logs/release_v6.1/tests_release_20261121.log`.

### Deferred
- Gladman-type tidal-locking time (next phase).
- English code comments in a separate copy (no logic changes), and a clean-clone reproducibility check.

## v5.2.x (2026-11-19)
Frozen v5.2 package, available from tags `v5.2-final`, `v5.2.1`, `v5.2.2`. Version DOI of v5.2.1: [10.5281/zenodo.24358996](https://doi.org/10.5281/zenodo.24358996). Its full changelog is in `CHANGELOG.md` at tag `v5.2.2`.
