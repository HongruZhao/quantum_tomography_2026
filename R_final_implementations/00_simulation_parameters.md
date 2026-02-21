# 00_simulation_parameters.R

## Purpose
Single source of truth for simulation parameters, profiles, CLI override handling, baseline policy mode resolution, and optional external config loading.

## Core entrypoints
- `load_simulation_config(profile, config_file, args)`
- `get_caseI_caseII_config(profile, config_file, args)`
- `get_case11_config(profile, config_file, args)`

## Profiles
- `full`: production-style all15 defaults (`alphas=0.5,0.9`, large `n_rep`, `n_total=200`, `n_trials=100`).
- `smoke_alpha0p9`: fast run (`alpha=0.9`, small `n_rep`, `n_total=40`).
- `case11`: dedicated case11 settings (`alpha=0.9`, larger `n_total` defaults).

## Baseline policy mode
- Internal policy key remains `uniform` for compatibility.
- User-facing baseline is permutation sampling.
- `baseline_policy_mode` accepted values:
  - `permutation_cyclic` (default): deterministic cyclic action order after initialization.
  - `uniform_random`: legacy random baseline.
- Legacy alias mapping is supported (`cyclic`, `random`, `uniform_policy_mode`).

## Parameter catalog

### Monte Carlo / design scale
- `alphas`: alpha values for true state mixture.
- `n_rep`: number of Monte Carlo replicates.
- `n_total`: adaptive steps (selection + sampling + MLE updates).
- `n_trials`: multinomial trials per adaptive step.
- `check_every`: MLE update interval (forced to 1 in final runs).

### Solver and selection strategy
- `solver_selection_mode`: `pilot_choose|cvx_fallback_pgd|always_dual|single_pgd`.
- `solver_selection_mode_resolved`: runtime-resolved mode.
- `production_solver`: selected production solver.
- `cvx_solver_preference`: preferred CVX backend (usually `SCS`).
- `solver`: legacy solver field for compatibility.
- `oracle_solver`: solver used for oracle design optimization.

### Pilot solver parameters
- `pilot_n_rep`, `pilot_n_total`: pilot workload.
- `pilot_cases`: case subset used in pilot.
- `pilot_alpha`: alpha used in pilot.

### Numerical stabilization parameters
- `eta_mle`: density eigenvalue floor in MLE.
- `ridge_selection`: ridge term for matrix inversions in selection/proxy.
- `fisher_eps`: epsilon inside FI probability denominator paths.
- `nll_eps_log`: additive epsilon in log-likelihood.
- `bures_tol`: tolerance checks in Bures square-root/eigenspectrum routines.
- `truncate_density`: whether to project estimated density back to feasible region.
- `allow_eta_zero`: allow zero floor (still only recommended for diagnostics).

### PGD controls
- `pgd_max_iter`: max PGD iterations per MLE call.
- `pgd_step`: optional explicit PGD step size (if enabled in solver).
- `pgd_tol`: optional PGD convergence tolerance.

### Runtime / reproducibility
- `n_workers`: parallel workers for replicate loops.
- `seed_base`: deterministic seed base.
- `resume`: enable checkpoint resume behavior.
- `l3_seed`: one-qubit random library seed.

### Artifacts and storage
- `save_full_trajectory`: store per-replicate trajectories.
- `trajectory_chunk_size`: trajectory chunk write size.
- `save_full_matrices`: keep full replicate matrices in results object.

### Plot controls
- `plot_1q_start`: x-axis start for one-qubit panels.
- `plot_2q_start`: x-axis start for two-qubit panels.
- `plot_t_end`: forced to `n_total` in final scripts.

### Output path controls
- `output_dir`: run output folder.
- `plots_dir`: plots output path (derived).
- `trajectories_dir`: trajectories output path (derived).

## Mathematical note for stabilization terms
Let `p(theta)` denote modeled probabilities and `J` total Fisher information.
- `nll_eps_log > 0`: uses `log(p + eps)` instead of `log(p)`; adds bias but avoids `-Inf` when probabilities become tiny.
- `fisher_eps > 0`: replaces `1/p` with `1/(p+eps)` in FI-related paths; reduces instability but distorts asymptotic scale.
- `ridge_selection > 0`: uses `(J + ridge I)^(-1)` to stabilize inversion; improves conditioning but changes proxy objective.
- `eta_mle > 0`: enforces density floor `rho >= eta I`; this is retained as the primary stabilization for estimator feasibility.

Final default philosophy in this package: keep non-essential stabilizers at zero unless explicitly requested.
