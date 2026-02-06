# Stabilized MLE solver (CVXR) for Bloch-parametrized tomography

This file specifies the maximum-likelihood estimator used in the adaptive loop.

The PDF recommends enforcing a uniform eigenvalue lower bound to stabilize:
- log-likelihood terms,
- Fisher information denominators,
- Bures metric evaluation.

## 1) Likelihood model

Given a measurement library \(\{Q_{a,b}\}\) and counts \(N_{a,b}\), define probabilities
\[
p_{a,b}(\rho)=\operatorname{tr}(\rho Q_{a,b}).
\]
The likelihood for the dataset \(D_n\) is:
\[
L(\rho) = \prod_{a=1}^k\prod_{b=1}^{r_a} p_{a,b}(\rho)^{N_{a,b}},
\]
and the negative log-likelihood (NLL) is:
\[
\mathcal L_{\mathrm{MLE}}(\rho) = -\sum_{a,b}N_{a,b}\log\operatorname{tr}(\rho Q_{a,b}).
\]

## 2) Stabilized feasible set (eigenvalue floor)

Fix \(\eta\) such that \(0<\eta<1/N\) and enforce:
\[
\rho\succeq \eta I,\quad \operatorname{tr}(\rho)=1,\quad \rho=\rho^\dagger.
\]
Denote this convex set by \(\mathcal D_\eta\).

The stabilized MLE is:
\[
\hat\rho_{\mathrm{MLE},\eta}\in\arg\min_{\rho\in\mathcal D_\eta}\mathcal L_{\mathrm{MLE}}(\rho).
\]

## 3) Implementation strategy used in CVXR1121.R

Rather than optimizing directly over complex Hermitian \(\rho\), CVXR1121.R optimizes over Bloch
coordinates \(\theta\) and enforces PSD by embedding \(\rho(\theta)\) as a real PSD matrix.

### 3.1 Theta-domain probability model

With `build_Sab_cab(sigmas, Q_list, ...)` we get:

- `c_ab` where \(c_{a,b}=\mathrm{tr}(Q_{a,b})/N\)
- `S_ab` where \([S_{ab}]_{(a,b),j} = s_j(a,b)=\mathrm{tr}(\sigma_j Q_{a,b})\)

Then:
\[
p_{a,b}(\theta)=c_{a,b}+\frac12 s(a,b)^\top\theta
\]
and the NLL becomes:
\[
\mathcal L(\theta) = -\sum_{a,b}N_{a,b}\log\Big(c_{a,b}+\frac12 s(a,b)^\top\theta\Big).
\]

### 3.2 PSD constraint via real embedding

CVXR1121.R builds:
- `Slist[[j]] = real_embed(sigmas[[j]])`
- `SI = diag(2*N)` (which is `real_embed(I_N)`)

and the embedded affine map:
\[
\rho_{\mathbb R}(\theta)=\frac{SI}{N}+\frac12\sum_j\theta_j Slist_j.
\]

It then introduces a CVXR PSD variable `rho` and constrains `rho == rho_R(theta)`.
This enforces \(\rho(\theta)\succeq 0\) in the complex sense.

## 4) Required modification: eigenvalue floor

We must additionally enforce:
\[
\rho(\theta)\succeq \eta I.
\]
In the real embedding, this becomes:
\[
\rho_{\mathbb R}(\theta) - \eta\,SI \succeq 0.
\]

### Implementation detail
Inside CVXR, add the constraint:

```r
constraints <- list(
  rho == A_affine,
  rho - eta*SI >> 0  # PSD constraint
)
```

If you keep `rho` as PSD already, the second line is the eigenvalue-floor margin.

### Proposed API change
Modify or wrap `fit_theta_cvxr`:

```r
fit_theta_cvxr <- function(..., eta = 0, eps_log = 1e-12, ...) {
  # if eta > 0: add rho - eta*SI >> 0
}
```

## 5) Numerical safeguards

### Log barrier epsilon (`eps_log`)
The code currently uses `log(p + eps_log)`. This is a pragmatic modeling choice to avoid `log(0)` if the
solver briefly evaluates out-of-domain points.

**Recommendation:** keep `eps_log` very small (e.g., 1e-12) and rely on the eigenvalue floor plus
POVM positivity to prevent true probabilities near zero.

### Probability positivity constraint (optional)
Optionally enforce `p >= p_min` directly in CVXR:
```r
constraints <- c(constraints, list(p_expr >= p_min))
```
but this is usually redundant if `eta > 0` and all POVM effects are nonzero.

### Solver choice
CVXR1121.R supports `MOSEK`, `ECOS`, `SCS`. For semidefinite constraints, `SCS` works but may require:
- tighter tolerances,
- more iterations,
- post-solve projection / hermitianization.

In simulation, prefer:
- MOSEK if available,
- otherwise SCS.

## 6) Outputs required from the estimator

Each MLE call should return:

- `theta_hat` (length d)
- `rho_hat` (N×N Hermitian)
- solver status, objective value, runtime (optional)

We will use:
- `theta_hat` for Fisher information and metric computations,
- `rho_hat` for Bures metric computations and probability sampling sanity checks.

