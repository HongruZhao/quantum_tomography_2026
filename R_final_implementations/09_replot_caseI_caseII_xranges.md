# 09_replot_caseI_caseII_xranges.R

## Purpose
Rebuild summary plots from an existing `full_results_caseI_caseII.rds` without rerunning simulation.

## Behavior
- Uses central config for x-window defaults:
  - 1q: `plot_1q_start..n_total`
  - 2q: `plot_2q_start..n_total`
- Keeps compatibility with legacy result schema using policy key `uniform`.
- User-facing legend text uses `Permutation`.

## CLI
- `--config_profile=<profile>`
- `--config_file=<path>`
- `--results_dir=<path>`
