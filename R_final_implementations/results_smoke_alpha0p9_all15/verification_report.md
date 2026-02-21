# Verification Report: Case I/II Multinomial Study

Results directory: `/Users/cyberslave/GitHub/quantum_tomography_2026/R_final_implementations/results_smoke_alpha0p9_all15`

## 1) Asymptotic formula alignment (Bures + Frobenius MSE)
Formula checked: E[D(rho_hat, rho)^2] ≈ (1/n) tr( G_D(theta) I(theta)^(-1) ).
Using `eq23_verification_frobenius_bures.csv`.

| metric | policy | mean abs gap (final) | q90 abs gap (final) | max abs gap (final) | mean abs gap (tail) | q90 abs gap (tail) | max abs gap (tail) |
|---|---|---:|---:|---:|---:|---:|---:|
| scaled_bures_squared | exact | 0.7094 | 1.6506 | 2.3690 | 0.6174 | 1.3345 | 1.8538 |
| scaled_frobenius_mse | exact | 0.2005 | 0.2819 | 0.3132 | 0.2965 | 0.5088 | 0.6090 |
| scaled_bures_squared | GI1 | 0.6077 | 0.8473 | 0.9750 | 0.6083 | 0.9726 | 0.9738 |
| scaled_frobenius_mse | GI1 | 0.1613 | 0.2886 | 0.4003 | 0.2096 | 0.4282 | 0.4660 |
| scaled_bures_squared | oracle_GI1 | 0.8538 | 1.8683 | 2.5995 | 0.6536 | 1.5276 | 2.1166 |
| scaled_frobenius_mse | oracle_GI1 | 0.4824 | 0.6974 | 0.8005 | 0.5180 | 0.6530 | 0.6701 |
| scaled_bures_squared | uniform | 0.6201 | 1.1686 | 1.3531 | 0.8089 | 1.5661 | 1.9428 |
| scaled_frobenius_mse | uniform | 0.6817 | 1.1549 | 1.1640 | 0.8911 | 1.9850 | 2.8008 |

## 2) Plot generation checks

- Total manifest plots: 18
- Missing manifest plots: 0

Required summary plots:
- `alpha_0p9_all15_relative_oracle_proxy.png`: OK
- `alpha_0p9_frobenius_mse_vs_oracle.png`: OK
- `alpha_0p9_bures_sq_vs_oracle.png`: OK

## 3) Policy ranking checks (Exact/GI1/OracleGI1 vs Permutation)

- Total alpha/case entries: 15
- Rank failures: 1

Failures:
- alpha=0.9 case=2q_A_bures comment=Exact<Permutation: NO | GI1<Permutation: NO | OracleGI1<Permutation: YES

## 4) Trajectory artifact checks

- Trajectory rows: 60
- Missing trajectory files: 0
- Missing policy trajectory rows: 0
- Permutation cyclic check rows: 15
- Permutation cyclic failures: 0
