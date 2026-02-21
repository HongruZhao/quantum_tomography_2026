# 10_run_case11_bures_rep1000_parallel.R

## Purpose
Dedicated high-replicate runner for Case 11 (`2q_A_bures`) with parallel replicate execution.

## Behavior
- Policies: `uniform`(permutation baseline), `exact`, `GI1`, `oracle_GI1`.
- Supports CVX fallback, PGD-only, and dual-compare solver modes.
- Outputs case11 proxy, relative-oracle, and Bures-vs-oracle plots.

## Baseline policy
- Default `baseline_policy_mode=permutation_cyclic`.
- Legacy alias `--uniform_policy_mode=` is still accepted.

## CLI
- `--config_profile=case11` (default)
- `--config_file=<path>`
- direct overrides like `--n_rep`, `--n_total`, `--alpha`, `--n_trials`, `--n_workers`.

## Outputs
- `case11_rep1000_proxy_curves.csv`
- `case11_rep1000_final_summary.csv`
- `case11_rep1000_full_results.rds`
- `case11_rep1000_config.rds`
- `plots/alpha_*_case11_*.png`
