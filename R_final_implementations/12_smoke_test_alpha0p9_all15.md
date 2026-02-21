# 12_smoke_test_alpha0p9_all15.R

## Purpose
Fast runnable smoke test for alpha=0.9 all-15 workflow.

## Pipeline executed
1. Run `07` in smoke profile.
2. Build `alpha_0p9_all15_plus_summary.pdf` via `11`.
3. Run structural verification via `08`.

## Default smoke profile
- `n_rep=2`
- `n_total=40`
- `n_trials=100`
- permutation baseline enabled

## Success criteria
Confirms existence of:
- `alpha_0p9_all15_relative_oracle_proxy.png`
- `alpha_0p9_frobenius_mse_vs_oracle.png`
- `alpha_0p9_bures_sq_vs_oracle.png`
- `alpha_0p9_all15_plus_summary.pdf`
