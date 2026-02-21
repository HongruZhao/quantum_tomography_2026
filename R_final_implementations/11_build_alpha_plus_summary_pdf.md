# 11_build_alpha_plus_summary_pdf.R

## Purpose
Create the compact 3-page `alpha_*_all15_plus_summary.pdf` from existing plot PNG files.

## Required inputs
For selected alpha tag (`0p5` or `0p9`), these files in `<results_dir>/plots`:
1. `alpha_*_all15_relative_oracle_proxy.png`
2. `alpha_*_frobenius_mse_vs_oracle.png`
3. `alpha_*_bures_sq_vs_oracle.png`

## Output
- `<results_dir>/alpha_*_all15_plus_summary.pdf`

## Notes
Uses macOS built-ins (`sips` + Automator `join`) to avoid extra package dependencies.
