# CHANGELOG: UHI (Unified Habitability Index)

## v5.2-final (2026-11-18), frozen release `v5.2-final`
This is the final version of the model. Licenses: code MIT (copyright Effy Viridis, 2026); data, results and figures CC-BY 4.0. PHL data are not redistributed (downloaded at run time under PHL terms). Archived on Zenodo from the GitHub release `v5.2-final` (`CITATION.cff`).

- **Smooth Yang et al. (2013) cloud-albedo blend.** The hard cutoff at S0 = 1000 W m⁻² in v5.1 is replaced by a linear ramp w_eff(S0) = clip((S0 − 800)/200, 0, 1). The model does not extrapolate below the tabulated range. The top 10 is identical and in the same order as v5.1; Spearman v5.1 vs v5.2 (live planets) = 0.999998.
- **Runaway-greenhouse veto (v5.1, Goldblatt et al. 2013, 282 W m⁻²).** The Teff boundary of the veto is now a separate parameter (`teff_rg`, default 3800 K = `teff_md`).
- **1-D sensitivity sweep of the veto** (`src/sweep_v52.py`). It passes for F_rg = 275–295 W m⁻² and for T_rg = 3850–3900 K. It fails the "top-10 composition invariant" criterion at T_rg = 3700 and 3750 K, where K2-125 b (Teff 3654 K) drops out of #10 and Ross 128 b enters. At the default T_rg = 3800 K the margin is UHI_rank 0.333 (K2-125 b) vs 0.323 (Ross 128 b); the other 9 members of the top 10 are invariant. Synthetic Venus has P_alive = 0 and Spearman > 0.99 at every primary sweep value. **Accepted as a limitation; T_rg was not retuned.**
- **ESI erratum.** v1–v5.1 used a total exponent of w/8 instead of w/4, which compressed ESI toward 1. This does not affect the UHI. ESI values in the earlier v4 and v5 documents are inflated.
- **Sobol coverage.** The post-v5 Sobol run (51 assumptions) does not include F_rg or the w_eff ramp; it was not re-run for v5.1 or v5.2.
- Automated tables (`src/make_tables.py`), vector figures (`src/make_figures.py`), an AAS MRT table (`src/make_mrt.py`), a master script (`run_uhi_v52.py`) and a Colab notebook.
- Certified run: Google Colab CPU, catalog md5 f2f70b24e67a9054940267e8081c81e8, 10 000 draws, seed 42. A local re-run is identical to within 1.9e-9.

## v5.1 (2026-11-18)
- Hybrid runaway-greenhouse veto: fHZ = 0 if S·S0/4·(1−A) > 282 W m⁻² (Goldblatt et al. 2013), except for tidally locked M dwarfs (Teff < 3800 K). Synthetic Venus moves from the 98.7th percentile to 0.

## v5 (2026-11-18): audit by R-Aetherios
- Orthogonal S_atm / S_obs, log-sigmoid μ(Rp) with a density override, an M-dwarf climate model (Yang+2013, Kopparapu+2016 k_in, Turbet+2018 CO2 floor), a hurdle ranking (P_alive, conditional median, UHI_rank, P(top10), rank CI), ε = 1e-12, Solar System anchors, and Sobol post-v5 (51 assumptions).

## v4: ELT/PCS reflected-light observability (Kasper et al. 2021). v3: detection probability Φ(SNR − u), cosmic shoreline. v2: Teff-dependent albedo (Shields et al. 2013). v1: first implementation.
