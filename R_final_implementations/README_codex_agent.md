# README_codex_agent.md

## Purpose
Prompt guide for using coding agents (for example Codex) to modify and rerun this simulation package safely.

## Ground rules for prompts
1. Always specify target folder: `R_final_implementations`.
2. State whether you want code edits, rerun, or both.
3. Include exact parameter overrides and output directory.
4. Ask the agent to keep policy key `uniform` while using permutation baseline.

## Example prompts

### A) Regenerate alpha 0.9 summary PDF only
```text
In /Users/cyberslave/GitHub/quantum_tomography_2026/R_final_implementations,
run 11_build_alpha_plus_summary_pdf.R using results_dir=<path> and alpha=0.9.
Do not modify code.
```

### B) Run full all15 with custom parameters
```text
Use 07_run_caseI_caseII_relative_oracle_multinomial.R with:
config_profile=full,
alphas=0.9,
n_rep=500,
n_total=300,
n_trials=100,
solver_selection_mode=single_pgd,
baseline_policy_mode=permutation_cyclic,
output_dir=<path>.
Then run 08_verify_caseI_caseII_outputs.R on the same output_dir.
```

### C) Case11 long run
```text
Use 10_run_case11_bures_rep1000_parallel.R with:
config_profile=case11,
alpha=0.9,
n_rep=1000,
n_total=1000,
n_trials=100,
solver_selection_mode=single_pgd,
baseline_policy_mode=permutation_cyclic,
output_dir=<path>.
After completion, report the two plot paths.
```

### D) Change stabilization parameters
```text
Edit only /Users/cyberslave/GitHub/quantum_tomography_2026/R_final_implementations/00_simulation_parameters.R:
set eta_mle=1e-4, ridge_selection=0, fisher_eps=0, nll_eps_log=0, bures_tol=0.
Do not touch other files. Then run the smoke test.
```

## Safe structural-edit prompt
```text
Refactor only files inside R_final_implementations.
Keep numbered .R/.md pairs intact.
Do not modify R_implementations_Feb.
After edits, run parse checks and 12_smoke_test_alpha0p9_all15.R.
```

## Quick checklist for agent output
Ask the agent to always provide:
1. Exact files modified.
2. Exact command(s) run.
3. Output directory.
4. Plot/PDF artifact paths.
5. Any failed checks.
