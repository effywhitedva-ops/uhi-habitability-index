# A Probabilistic Unified Habitability Index Coupling Planetary Viability with JWST and ELT/PCS Observability across 6,112 Exoplanets

**Effy Viridis** (ORCID: [0009-0004-5620-5352](https://orcid.org/0009-0004-5620-5352))
*Independent Researcher, El Salvador*
Email: effywhitedva@gmail.com

*Keywords:* Astrobiology (74) — Computational methods (1965) — Direct imaging (387) — Exoplanets (498) — Habitable planets (695) — Transmission spectroscopy (2133)

## Abstract

Traditional habitability metrics such as the Earth Similarity Index (ESI) rank exoplanets by bulk physical analogy without accounting for atmospheric characterization feasibility. We present version 5.2-final of the Unified Habitability Index (UHI), a probabilistic, hurdle-based framework evaluated via 10,000 Monte Carlo draws per planet across 6,112 scored exoplanets from the NASA Exoplanet Archive. UHI combines seven physical, climatic, stellar, and instrumental sub-indices with Dirichlet weight propagation, synchronous M-dwarf climate parameterizations, a hybrid runaway greenhouse veto, and physical observability models for JWST and ELT/PCS. At a 5σ detection threshold, 515 planets remain viable (P_alive > 0). Synthetic Solar System anchors at 10 pc pass cleanly (Earth at the 99.30th (non-transiting) and 99.26th (transiting) percentiles; Venus and Mars = 0). Compared against ESI, UHI filters out 71.2% of high-ESI (≥ 0.8) targets that are currently unobservable while elevating seven nearby non-transiting worlds into the Top 20 via reflected-light imaging.

## 1. Introduction

Prioritizing exoplanets for atmospheric characterization requires balancing intrinsic planetary habitability against realistic instrumental detectability. Widely used similarity indices such as the Earth Similarity Index (ESI; Schulze-Makuch et al. 2011) and classical Habitable Zone (HZ) boundaries (Kopparapu et al. 2013, 2014) successfully identify Earth-sized temperate candidates, yet they treat characterization feasibility as an afterthought. Consequently, many top-ranked analogs lie beyond the reach of current or near-future observatories. Here we summarize version 5.2-final of the Unified Habitability Index (UHI), a fully open, reproducible Monte Carlo framework designed to bridge intrinsic habitability and observational priority for both transiting (JWST) and non-transiting (ELT/PCS) exoplanets.

## 2. Methodology

We evaluate 6,112 scored exoplanets from a dated snapshot (2026-11-18) of the NASA Exoplanet Archive Planetary Systems Composite Parameters table (pscomppars; 6,470 rows, 358 excluded for missing core parameters; Akeson et al. 2013) using 10,000 Monte Carlo draws per planet (seed 42). Measurement uncertainties are sampled from split-normal distributions, with missing masses or radii propagated via the probabilistic mass–radius relation of Chen & Kipping (2017).

The composite index is defined as a weighted geometric mean, UHI = ∏ S_i^(w_i), over seven sub-indices—physical (S_fis), climate (S_clima), atmospheric retention (S_atm), stellar environment (S_estrella), geophysical (S_geo), magnetic (S_mag), and observability (S_obs)—with weights drawn per realization from Dirichlet(α = [2, 2, 2, 2, 1, 1, 2]). A numerical floor ε = 10⁻¹² precedes log-aggregation, imposing an exact veto (UHI = 0) whenever any S_i ≤ ε. Key architectural features of v5.2-final include:

- **Orthogonal Atmospheric–Observational Decoupling:** S_atm models intrinsic atmospheric retention via the bolometric Cosmic Shoreline (Zahnle & Catling 2017), while all instrumental detection terms reside strictly in S_obs = P_obs × Availability(tier). Transiting observability combines JWST transmission (TSM-scaled with a 20 ppm NIRISS systematic floor; Kempton et al. 2018; Greene et al. 2016) and secondary eclipse emission (ESM), whereas non-transiting targets use a residual-halo photon-noise limited ELT/PCS reflected-light model (Kasper et al. 2021) with inner working angle IWA ~ U(2, 4) λ/D at 0.75 μm. The primary threshold is 5σ, with a 3σ comparison variant.

- **Synchronous M-Dwarf Climate and Hybrid Runaway Veto:** S_clima couples a 1-D two-mode energy balance model (North et al. 1981) and T_eff-dependent Bond albedos (Shields et al. 2013) with 3-D GCM-calibrated parameterizations for tidally locked M-dwarf planets (T_eff < 3800 K): substellar cloud albedos (Yang et al. 2013) smoothly blended via a linear ramp w_eff(S_0) = clip((S_0 − 800)/200, 0, 1), inner-edge expansion k_in ~ U(1.5, 2.0) (Kopparapu et al. 2016), and an outer-edge CO₂ surface-pressure floor (Turbet et al. 2018). For non-synchronous or T_eff ≥ 3800 K planets, a hybrid runaway greenhouse veto enforces f_HZ = 0 whenever absorbed flux exceeds 282 W m⁻² (Goldblatt et al. 2013, unchanged by the 2014 corrigendum).

- **Smooth Mean Molecular Weight and Hurdle Ranking:** Across the radius valley (Fulton et al. 2017), μ(R_p) transitions log-sigmoidally between secondary (28 or 44 g mol⁻¹) and H₂/He (2.3 g mol⁻¹) envelopes centered at R_0 ~ N(1.65, 0.08) R⊕ (ΔR = 0.12 R⊕), overridden to secondary atmospheres when measured densities exceed 5.5 g cm⁻³ (M/σ_M > 3). Because multiplicative vetoes yield bimodal posteriors, planets are ranked by a continuous hurdle metric: UHI_rank = P_alive × median(UHI | UHI > 0), accompanied by P(top 10) and 90% rank credible intervals.

## 3. Results, Validation, and Benchmarking

Across the 6,112 planets, 515 remain viable (P_alive > 0) at both 5σ and 3σ (Spearman ρ = 0.985 between variants). Synthetic Solar System twins injected at 10 pc validate the pipeline: Earth ranks at the 99.30th percentile (#44, non-transiting) and 99.26th (#46, transiting) with P_alive = 1.00, whereas Venus and Mars yield UHI = 0 (P_alive = 0).

Table 1 summarizes the Top 10 targets at 5σ alongside their 3σ rankings and corrected ESI values (Schulze-Makuch et al. 2011). At 5σ, temperate sub-Neptunes around M dwarfs—TOI-4336 A b (UHI_rank = 0.562) and TOI-7166 b (0.547)—lead due to large H₂-rich scale heights and an upper radius trapezoid limit of f_R,d = 2.5 R⊕. Relaxing to 3σ (a proxy for deeper JWST campaigns) elevates rocky worlds: TRAPPIST-1 d climbs from #5 (0.436) to #1 (0.609), and TRAPPIST-1 e rises from #17 (0.272) to #9 (0.428). This shift directly reflects the observational trade-off identified by a 51-parameter Sobol sensitivity analysis (SALib; Herman & Usher 2017) performed after v5 and not rerun for v5.1 or v5.2. The quoted total-effect indices are means across that run's 13 per-planet outputs and 3 aggregate metrics: the detection threshold snr_thr (mean S_T = 0.45) and f_R,d (mean S_T = 0.15) dominate variance, while the newly added transition parameters (R_0, ΔR, k_in) have mean S_T < 0.01. The Sobol analysis does not cover the runaway-greenhouse flux threshold F_rg or the substellar-cloud mixing ramp w_eff.

**Table 1.** *Top 10 UHI v5.2-final target exoplanets at 5σ detection threshold with 3σ comparison metrics, conditional posterior distributions, and classical ESI.*

| Rank (5σ) | Rank (3σ) | Planet | UHI_rank (5σ) | P_alive | Cond. Median [IQR] | P(top 10) | Tier | ESI_SM11 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 2 | TOI-4336 A b | 0.562 | 0.99 | 0.567 [0.460–0.665] | 0.60 | T0 | 0.741 |
| 2 | 3 | TOI-7166 b | 0.547 | 0.99 | 0.554 [0.456–0.647] | 0.56 | T1 | 0.816 |
| 3 | 4 | GJ 273 b | 0.477 | 0.87 | 0.546 [0.424–0.642] | 0.47 | T2 | 0.847 |
| 4 | 6 | K2-332 b | 0.449 | 0.91 | 0.494 [0.393–0.584] | 0.36 | T2 | 0.786 |
| 5 | 1 | TRAPPIST-1 d | 0.436 | 1.00 | 0.436 [0.261–0.602] | 0.35 | T0 | 0.905 |
| 6 | 8 | Wolf 1061 c | 0.434 | 0.82 | 0.529 [0.426–0.625] | 0.41 | T2 | 0.823 |
| 7 | 11 | Proxima Cen d | 0.389 | 0.84 | 0.463 [0.356–0.566] | 0.28 | T2 | 0.876 |
| 8 | 5 | TOI-2094 b | 0.377 | 0.99 | 0.382 [0.229–0.521] | 0.23 | T2 | 0.813 |
| 9 | 17 | GJ 1061 d | 0.342 | 0.64 | 0.537 [0.429–0.611] | 0.33 | T2 | 0.884 |
| 10 | 15 | K2-125 b | 0.333 | 0.84 | 0.396 [0.271–0.502] | 0.17 | T2 | 0.767 |

The 90% rank credible intervals are wide (e.g., 1–30 for the leading target) and are provided in the archive.

Benchmarking against ESI_SM11 and the Habitable Worlds Catalog (PHL @ UPR Arecibo 2024, file last updated 2024-03-12) confirms that UHI provides distinct prioritization (Spearman ρ = 0.466 vs. ESI_SM11 and 0.444 vs. PHL ESI among viable worlds). Specifically, 71.2% (47/66) of planets with ESI_SM11 ≥ 0.8 have median detection probability P_det < 0.001 at 5σ (72.0%, 18/25 for PHL ESI); UHI demotes these unobservable analogs so that none enter the top 29 (best: Kepler-1653 b at #30). Simultaneously, UHI rescues 26 nearby non-transiting viable worlds with median P_det ≥ 0.5 via ELT/PCS, placing seven in the Top 20 (GJ 273 b, Wolf 1061 c, Proxima Cen d, GJ 1061 d, Ross 128 b, Proxima Cen b, and GJ 433 d).

A 1-D sweep of the runaway veto across F_rg ∈ [275, 295] W m⁻² preserves Venus P_alive = 0, an identical Top 10, and Spearman ρ ≥ 0.996 among viable worlds. Sweeping the thermal boundary to T_rg = 3700–3750 K swaps #10 K2-125 b (T_eff = 3654 K, UHI_rank = 0.333) with #11 Ross 128 b (0.323) while leaving the other nine unchanged; we retain T_rg = 3800 K without post-hoc tuning and report this boundary sensitivity transparently.

## 4. Reproducibility, Software, and Acknowledgments

All code (MIT license), dated catalog snapshots, machine-readable tables, and figures (CC BY 4.0) for UHI v5.2-final are archived at https://github.com/effywhitedva-ops/uhi-habitability-index and Zenodo (doi:10.5281/zenodo.24358422). The PHL catalog is fetched at runtime and not redistributed.

Software: Python 3.13, NumPy, SciPy, pandas, SALib (Herman & Usher 2017), Matplotlib, Astropy, PyArrow. This research has made use of the NASA Exoplanet Archive, operated by the California Institute of Technology under contract with NASA under the Exoplanet Exploration Program.

**AI disclosure:** AI agents based on Grok (via Grok Bot; agents "Shoggothsito" and "Aetherios"), with computing tools and under the author's supervision, were used to implement and refactor the Python code, run Monte Carlo and sensitivity simulations on local and Google Colab runtimes, generate tables and figures from the author's data with Matplotlib, audit methods, assist in drafting, editing, and formatting the manuscript, and supplement the literature search. Gemini was used as a reviewer to provide critical feedback on the methodology and manuscript drafts; all suggestions were evaluated by the author. All references were checked against the original sources and NASA ADS. The author reviewed and approved every change and takes full responsibility for the content. AI tools are not authors.

## References

- Akeson, R. L., Chen, X., Ciardi, D., et al. 2013, PASP, 125, 989
- Chen, J., & Kipping, D. 2017, ApJ, 834, 17
- Fulton, B. J., Petigura, E. A., Howard, A. W., et al. 2017, AJ, 154, 109
- Goldblatt, C., Robinson, T. D., Zahnle, K. J., & Crisp, D. 2013, NatGeo, 6, 661
- Goldblatt, C., Robinson, T. D., Zahnle, K. J., & Crisp, D. 2014, NatGeo, 7, 931 (Corrigendum)
- Greene, T. P., Line, M. R., Montero, C., et al. 2016, ApJ, 817, 17
- Herman, J., & Usher, W. 2017, JOSS, 2, 97
- Kasper, M., Cerpa Urra, N., Pathak, P., et al. 2021, Msngr, 182, 38
- Kempton, E. M.-R., Bean, J. L., Louie, D. R., et al. 2018, PASP, 130, 114401
- Kopparapu, R. K., Ramirez, R., Kasting, J. F., et al. 2013, ApJ, 765, 131
- Kopparapu, R. K., Ramirez, R. R., SchottelKotte, J., et al. 2014, ApJ, 787, L29
- Kopparapu, R. K., Wolf, E. T., Haqq-Misra, J., et al. 2016, ApJ, 819, 84
- PHL @ UPR Arecibo 2024, The Habitable Worlds Catalog (file last updated 2024-03-12), https://phl.upr.edu/hwc
- North, G. R., Cahalan, R. F., & Coakley, J. A. 1981, RGeo, 19, 91
- Schulze-Makuch, D., Méndez, A., Fairén, A. G., et al. 2011, Astrobiology, 11, 1041
- Shields, A. L., Meadows, V. S., Bitz, C. M., et al. 2013, Astrobiology, 13, 715
- Turbet, M., Bolmont, E., Leconte, J., et al. 2018, A&A, 612, A86
- Yang, J., Cowan, N. B., & Abbot, D. S. 2013, ApJL, 771, L45
- Zahnle, K. J., & Catling, D. C. 2017, ApJ, 843, 122
