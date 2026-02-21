# R_final_implementations

## What this package contains
This folder is the finalized simulation package for:
- all-15-case alpha studies (Case I and Case II),
- case11 dedicated runs,
- summary PDF generation,
- output verification,
- runnable smoke testing.

Internal policy key `uniform` is preserved for compatibility, but baseline behavior is now permutation-cyclic by default.

## Numbered structure
- `00_*`: central parameter control
- `01_*` to `06_*`: shared math/model libraries
- `07_*`: main all-15 runner
- `08_*`: verification
- `09_*`: replot existing results
- `10_*`: case11 dedicated runner
- `11_*`: 3-page alpha summary PDF builder
- `12_*`: fast smoke test orchestrator

Each numbered `.R` has a matching numbered `.md`.

## Central parameter control
Edit or override values from `00_simulation_parameters.R`.

High-impact parameters:
- simulation size: `alphas`, `n_rep`, `n_total`, `n_trials`
- baseline: `baseline_policy_mode`
- numerical controls: `eta_mle`, `ridge_selection`, `fisher_eps`, `nll_eps_log`, `bures_tol`
- solver controls: `solver_selection_mode`, `cvx_solver_preference`, `pgd_max_iter`
- runtime controls: `n_workers`, `resume`
- plotting: `plot_1q_start`, `plot_2q_start`, `plot_t_end`

## How to run
Run from this folder or provide absolute script paths.

### 1) Full all-15 alpha run
```bash
Rscript 07_run_caseI_caseII_relative_oracle_multinomial.R \
  --config_profile=full \
  --output_dir=/Users/cyberslave/GitHub/quantum_tomography_2026/R_final_implementations/results_final_caseI_caseII_all15
```

### 2) Build 3-page summary PDF for alpha 0.9
```bash
Rscript 11_build_alpha_plus_summary_pdf.R \
  --results_dir=/Users/cyberslave/GitHub/quantum_tomography_2026/R_final_implementations/results_final_caseI_caseII_all15 \
  --alpha=0.9
```

### 3) Verify outputs
```bash
Rscript 08_verify_caseI_caseII_outputs.R \
  --results_dir=/Users/cyberslave/GitHub/quantum_tomography_2026/R_final_implementations/results_final_caseI_caseII_all15
```

### 4) Fast smoke test
```bash
Rscript 12_smoke_test_alpha0p9_all15.R \
  --output_dir=/Users/cyberslave/GitHub/quantum_tomography_2026/R_final_implementations/results_smoke_alpha0p9_all15
```

## How it works (short)
1. Build basis and measurement libraries (`01`-`04`).
2. Simulate multinomial sampling and adaptive selection with per-step MLE (`07`/`10`).
3. Compute scaled proxy and loss metrics and oracle references.
4. Write trajectories, summaries, and plots.
5. Build compact summary PDFs (`11`) and validate outputs (`08`).

## Troubleshooting
- If CVX is unavailable, use `--solver_selection_mode=single_pgd`.
- If run is too slow, reduce `n_rep` and `n_total`.
- If memory is high, set `--save_full_trajectory=0`.
- For reproducibility, keep `seed_base` fixed.
