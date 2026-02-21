# 07_run_caseI_caseII_relative_oracle_multinomial.R

## Purpose
Primary all-15-case runner for Case I/Case II alpha studies with multinomial batch sampling.

## Coverage
- 15 cases: 1q(L1/L2/L3) + 2q(A/B) crossed with (frobenius, bures, observable).
- Policies: `uniform` (permutation baseline), `exact`, `GI1`, `oracle_GI1`.
- Per-step MLE update and trajectory logging support.

## Baseline policy behavior
- Policy key remains `uniform`.
- Actual baseline action rule is permutation-cyclic by default:
  - after initialization, `a_t = ((t - n_init - 1) %% k) + 1`.
- Optional random fallback is available through `baseline_policy_mode`.

## CLI
- `--config_profile=full|smoke_alpha0p9|case11`
- `--config_file=<path>` (optional `.rds` or list expression file)
- Any parameter override from `00_simulation_parameters.R` can be passed as `--name=value`.

## Main outputs
- `oracle_limits_caseI_caseII.csv`
- `curve_summary_caseI_caseII.csv`
- `eq23_verification_frobenius_bures.csv`
- `policy_comparison_comments.csv`
- `plot_manifest.csv`
- `trajectory_manifest.csv`
- `full_results_caseI_caseII.rds`
- `plots/*.png`
- `trajectories/*.rds`
