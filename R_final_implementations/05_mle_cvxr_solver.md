# 05_mle_cvxr_solver.R

## Purpose
MLE solver stack for density estimation under linear measurement model.

## Main functionality
- CVX-based solver path and PGD path.
- Dual compare helper (CVX vs PGD on same counts) for solver diagnostics.
- Density floor enforcement utilities (`rho >= eta I`) retained as primary stabilization.
- Solver status and train-loss diagnostics used in trajectory logging.

## Key controlled parameters
- `eta_mle`, `nll_eps_log`, `pgd_max_iter`, `truncate_density`, `allow_eta_zero`.
