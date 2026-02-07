# Validation and sanity checks (recommended)

This file lists checks that should be implemented early to avoid silent math/implementation errors.

## 1) State parametrization checks

### 1.1 `rho_of_theta` produces Hermitian trace-1 matrices
For random theta vectors:
- `rho <- rho_of_theta(theta, sigmas, N)`
- Check `max(Mod(rho - Conj(t(rho)))) < tol`
- Check `abs(Re(sum(diag(rho))) - 1) < tol`

### 1.2 PSD and eigenvalue floor feasibility
For stabilized MLE with standard floor `eta = 1e-3`:
- ensure `eta < 1/N`
- after solving, check `min(eigen(rho_hat)$values) >= 1e-3 - tol`

## 2) Measurement library checks

For each setting a:

### 2.1 POVM completeness
Check:
\[
\sum_{b=1}^{r_a} Q_{a,b} = I_N
\]
numerically (Frobenius norm of difference small).

### 2.2 PSD of each effect
Check `min(eigen(Q)$values) >= -tol`.

### 2.3 For PVMs: idempotence and orthogonality (optional)
For rank-1 or projective settings:
- \(Q_{a,b}^2 \approx Q_{a,b}\)
- \(Q_{a,b}Q_{a,b'}\approx 0\) for \(b\neq b'\)

Library A should yield rank-2 projectors.
Library B should yield rank-1 projectors.

## 3) Born probabilities sanity

For random rho and each setting a:
- compute `p <- born_probs_list(rho, Q_list)[[a]]`
- check:
  - `all(p >= -tol)`
  - `abs(sum(p) - 1) < tol`

If using stabilized states (\(\rho\succeq\eta I\) with \(\eta=10^{-3}\)), also check:
- `min(p) >= 1e-3 * min(tr(Q_{a,b}))` up to tolerances

## 4) Fisher information sanity

### 4.1 FI is PSD
For random theta, compute \(I_a(\theta)\) and check eigenvalues are nonnegative up to tolerance.

### 4.2 Total FI monotonicity
When you add a shot at setting a, \(J_n\) should update by adding \(I_a(\hat\theta_n)\) (PSD), so:
- \(J_{n+1}-J_n \succeq 0\) approximately.
This is a sanity check on the update logic.

## 5) Bures metric checks

### 5.1 Symmetry and PSD
`G_bures` should be symmetric.
It should be positive definite for full-rank rho (in exact math). Numerically:
- eigenvalues > 0 up to tolerance.

### 5.2 Consistency with Lyapunov equation
Verify that `lyapunov_inv` returns X such that:
\[
\frac12(\rho X + X\rho)\approx Y.
\]

## 6) MUB checks (Library B)

Check mutual unbiasedness numerically:

For each pair of bases \(B_a, B_{a'}\), \(a\neq a'\):
- for each |ψ> in Ba and |φ> in Ba':
  - compute `abs(Conj(t(psi)) %*% phi)^2`
  - should be close to `1/4`.

Also check each basis is orthonormal:
- `Conj(t(Psi)) %*% Psi ≈ I_4`, where Psi columns are basis vectors.

## 7) Selector sanity

For a fixed theta and counts:

- exact score should be finite once J invertible.
- GI1 score should be finite once J invertible.
- If G is identity, the Frobenius exact selector should match the existing unweighted A-opt selector.

## 8) Regression test: uniform policy invariance

Under uniform policy, the mean risk curve should be invariant to:
- relabeling settings,
- reordering outcomes b within each setting,
provided the implementation uses the mapping `ab_df` consistently.

This is a powerful way to catch indexing bugs.
