# Análisis Fase 1 UHI v6 (2026-11-20, hora de El Salvador)
Scripts en `src/analysis/` (a1 Sobol, a2 prior, a3 Innes vs Madhusudhan, a4 recálculo a mano, a6 fallbacks). El resumen está en `CHANGELOG_v6.md` («Análisis de Fase 1»).
- `sobol_robust_{rama}_{5,3}s.csv`, `sobol_closing_criterion.json`, `sobol_output_diagnostics.csv`
- `prior_robustness_all_branches.csv`
- `innes_vs_madhu_formal.json`, `innes_tail_artifact.json`
- `hand_recompute_vs_pipeline.csv`, `hand_recompute_values.json`
- `anchors_recheck.json`, `fallback_checks.json`, `zeng_md5_check.txt`, `compat_origin_check.txt`
- Figuras: `figures/analysis_phase1/fig_a1_sobol_ST_robusto.png`, `figures/analysis_phase1/fig_a3_innes_vs_madhu.png`
