## ============================================================================
## 12_smoke_test_alpha0p9_all15.R
## Fast end-to-end smoke test for alpha=0.9 all-15 workflow.
## Runs:
##   - 07_run_caseI_caseII_relative_oracle_multinomial.R
##   - 11_build_alpha_plus_summary_pdf.R
##   - 08_verify_caseI_caseII_outputs.R
## ============================================================================

parse_cli_value <- function(prefix, args = commandArgs(trailingOnly = TRUE), default = NULL) {
  hit <- grep(paste0("^", prefix), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", prefix), "", hit[1])
}

script_dir <- tryCatch({
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  NULL
}, error = function(e) NULL)
if (is.null(script_dir) || !nzchar(script_dir)) script_dir <- getwd()
if (!file.exists(file.path(script_dir, "00_simulation_parameters.R"))) {
  candidate <- file.path(script_dir, "R_final_implementations")
  if (file.exists(file.path(candidate, "00_simulation_parameters.R"))) script_dir <- candidate
}
setwd(script_dir)

source(file.path(script_dir, "00_simulation_parameters.R"))

args <- commandArgs(trailingOnly = TRUE)
config_file <- parse_cli_value("--config_file=", args = args, default = NULL)
cfg <- load_simulation_config(profile = "smoke_alpha0p9", config_file = config_file, args = args)

output_dir <- parse_cli_value("--output_dir=", args = args, default = cfg$output_dir %||% "results_smoke_alpha0p9_all15")
run_r <- file.path(script_dir, "07_run_caseI_caseII_relative_oracle_multinomial.R")
pdf_r <- file.path(script_dir, "11_build_alpha_plus_summary_pdf.R")
verify_r <- file.path(script_dir, "08_verify_caseI_caseII_outputs.R")

if (!file.exists(run_r)) stop(sprintf("Missing runner: %s", run_r))
if (!file.exists(pdf_r)) stop(sprintf("Missing pdf builder: %s", pdf_r))
if (!file.exists(verify_r)) stop(sprintf("Missing verifier: %s", verify_r))

run_args <- c(
  run_r,
  "--config_profile=smoke_alpha0p9",
  sprintf("--output_dir=%s", output_dir),
  "--alphas=0.9",
  "--baseline_policy_mode=permutation_cyclic",
  sprintf("--n_rep=%d", as.integer(cfg$n_rep)),
  sprintf("--n_total=%d", as.integer(cfg$n_total)),
  sprintf("--n_trials=%d", as.integer(cfg$n_trials))
)

cat("[1/3] Running smoke simulation...\n")
st1 <- system2("Rscript", run_args)
if (!identical(st1, 0L)) stop("Smoke simulation failed.")

cat("[2/3] Building alpha_0p9 summary PDF...\n")
st2 <- system2("Rscript", c(
  pdf_r,
  "--config_profile=smoke_alpha0p9",
  sprintf("--results_dir=%s", output_dir),
  "--alpha=0.9"
))
if (!identical(st2, 0L)) stop("Smoke PDF build failed.")

cat("[3/3] Running structural verification...\n")
st3 <- system2("Rscript", c(
  verify_r,
  "--config_profile=smoke_alpha0p9",
  sprintf("--results_dir=%s", output_dir)
))
if (!identical(st3, 0L)) stop("Smoke verification failed.")

required <- c(
  file.path(output_dir, "plots", "alpha_0p9_all15_relative_oracle_proxy.png"),
  file.path(output_dir, "plots", "alpha_0p9_frobenius_mse_vs_oracle.png"),
  file.path(output_dir, "plots", "alpha_0p9_bures_sq_vs_oracle.png"),
  file.path(output_dir, "alpha_0p9_all15_plus_summary.pdf")
)
missing <- required[!file.exists(required)]
if (length(missing) > 0) {
  stop(sprintf("Smoke test completed but missing artifacts:\n%s", paste(missing, collapse = "\n")))
}

cat("Smoke test completed. Key artifacts:\n")
for (p in required) cat(sprintf("  - %s\n", p))
