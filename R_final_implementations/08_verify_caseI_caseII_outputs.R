## ============================================================================
## 08_verify_caseI_caseII_outputs.R
## Verification checks for Case I/II multinomial tomography outputs.
## ============================================================================

parse_cli_value <- function(prefix, args = commandArgs(trailingOnly = TRUE), default = NULL) {
  hit <- grep(paste0("^", prefix), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", prefix), "", hit[1])
}

read_csv_or_empty <- function(path) {
  tryCatch({
    lines <- readLines(path, warn = FALSE)
    lines <- lines[nzchar(trimws(lines))]
    if (length(lines) == 0L) return(data.frame())
    if (length(lines) == 1L && trimws(lines[1]) %in% c("\"\"", "")) return(data.frame())
    read.csv(path, stringsAsFactors = FALSE)
  }, error = function(e) data.frame())
}

results_dir <- parse_cli_value(
  "--results_dir=",
  default = NA_character_
)

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

cfg_verify <- load_simulation_config(
  profile = parse_cli_value("--config_profile=", default = "full"),
  config_file = parse_cli_value("--config_file=", default = NULL),
  args = commandArgs(trailingOnly = TRUE)
)

if (is.na(results_dir) || !nzchar(results_dir)) {
  results_dir <- cfg_verify$output_dir %||% "/Users/cyberslave/GitHub/quantum_tomography_2026/R_final_implementations/results_final_caseI_caseII_all15"
}

if (!dir.exists(results_dir)) {
  stop(sprintf("Missing results_dir: %s", results_dir))
}

oracle_file <- file.path(results_dir, "oracle_limits_caseI_caseII.csv")
eq23_file <- file.path(results_dir, "eq23_verification_frobenius_bures.csv")
comment_file <- file.path(results_dir, "policy_comparison_comments.csv")
plot_manifest_file <- file.path(results_dir, "plot_manifest.csv")
traj_manifest_file <- file.path(results_dir, "trajectory_manifest.csv")

for (f in c(oracle_file, eq23_file, comment_file, plot_manifest_file, traj_manifest_file)) {
  if (!file.exists(f)) stop(sprintf("Missing required file: %s", f))
}

oracle_df <- read_csv_or_empty(oracle_file)
eq23_df <- read_csv_or_empty(eq23_file)
comment_df <- read_csv_or_empty(comment_file)
plot_df <- read_csv_or_empty(plot_manifest_file)
traj_df <- read_csv_or_empty(traj_manifest_file)

## 1) Asymptotic formula matching checks (Bures + Frobenius MSE)
eq23_rows <- data.frame(
  metric = character(0),
  policy = character(0),
  mean_abs_rel_gap_final = numeric(0),
  q90_abs_rel_gap_final = numeric(0),
  max_abs_rel_gap_final = numeric(0),
  mean_abs_rel_gap_tail = numeric(0),
  q90_abs_rel_gap_tail = numeric(0),
  max_abs_rel_gap_tail = numeric(0),
  stringsAsFactors = FALSE
)

if (nrow(eq23_df) > 0 && all(c("rel_gap_final", "rel_gap_tail", "metric", "policy") %in% names(eq23_df))) {
  eq23_df$abs_rel_gap_final <- abs(eq23_df$rel_gap_final)
  eq23_df$abs_rel_gap_tail <- abs(eq23_df$rel_gap_tail)

  eq23_summary <- aggregate(
    cbind(abs_rel_gap_final, abs_rel_gap_tail) ~ metric + policy,
    data = eq23_df,
    FUN = function(x) c(
      mean = mean(x, na.rm = TRUE),
      q90 = as.numeric(quantile(x, probs = 0.9, na.rm = TRUE)),
      max = max(x, na.rm = TRUE)
    )
  )

  eq23_rows <- do.call(rbind, lapply(seq_len(nrow(eq23_summary)), function(i) {
    data.frame(
      metric = eq23_summary$metric[i],
      policy = eq23_summary$policy[i],
      mean_abs_rel_gap_final = eq23_summary$abs_rel_gap_final[i, "mean"],
      q90_abs_rel_gap_final = eq23_summary$abs_rel_gap_final[i, "q90"],
      max_abs_rel_gap_final = eq23_summary$abs_rel_gap_final[i, "max"],
      mean_abs_rel_gap_tail = eq23_summary$abs_rel_gap_tail[i, "mean"],
      q90_abs_rel_gap_tail = eq23_summary$abs_rel_gap_tail[i, "q90"],
      max_abs_rel_gap_tail = eq23_summary$abs_rel_gap_tail[i, "max"],
      stringsAsFactors = FALSE
    )
  }))
}

## 2) Plot existence checks
if (nrow(plot_df) > 0 && "path" %in% names(plot_df)) {
  plot_df$exists <- file.exists(plot_df$path)
  missing_plots <- plot_df[!plot_df$exists, , drop = FALSE]
} else {
  missing_plots <- data.frame()
}

alphas_present <- numeric()
if (nrow(oracle_df) > 0 && "alpha" %in% names(oracle_df)) {
  alphas_present <- sort(unique(as.numeric(oracle_df$alpha)))
}
if (length(alphas_present) == 0L) {
  alphas_present <- c(0.5, 0.9)
}
alpha_tags <- gsub("\\.", "p", format(alphas_present, nsmall = 1, trim = TRUE))

required_named_plots <- as.vector(rbind(
  sprintf("alpha_%s_all15_relative_oracle_proxy.png", alpha_tags),
  sprintf("alpha_%s_frobenius_mse_vs_oracle.png", alpha_tags),
  sprintf("alpha_%s_bures_sq_vs_oracle.png", alpha_tags)
))

named_plot_exists <- data.frame(
  filename = required_named_plots,
  path = file.path(results_dir, "plots", required_named_plots),
  exists = file.exists(file.path(results_dir, "plots", required_named_plots)),
  stringsAsFactors = FALSE
)

## 3) Policy ranking checks (Exact / GI1 / OracleGI1 vs Permutation)
rank_fail <- data.frame()
if (nrow(comment_df) > 0) {
  if (!("exact_better" %in% names(comment_df))) comment_df$exact_better <- FALSE
  if (!("gi1_better" %in% names(comment_df))) comment_df$gi1_better <- FALSE
  if (!("oracle_gi1_better" %in% names(comment_df))) comment_df$oracle_gi1_better <- FALSE
  comment_df$exact_better <- as.logical(comment_df$exact_better)
  comment_df$gi1_better <- as.logical(comment_df$gi1_better)
  comment_df$oracle_gi1_better <- as.logical(comment_df$oracle_gi1_better)
  rank_fail <- comment_df[!(comment_df$exact_better & comment_df$gi1_better & comment_df$oracle_gi1_better), , drop = FALSE]
}

## 4) Trajectory artifact checks
traj_missing <- data.frame()
missing_policy_df <- data.frame()
baseline_seq_check <- data.frame(
  alpha = numeric(0),
  case_id = character(0),
  policy = character(0),
  pass_cyclic = logical(0),
  checked_from_t = integer(0),
  checked_to_t = integer(0),
  detail = character(0),
  stringsAsFactors = FALSE
)

if (nrow(traj_df) > 0 && all(c("path", "alpha", "case_id", "policy") %in% names(traj_df))) {
  traj_df$exists <- file.exists(traj_df$path)
  traj_missing <- traj_df[!traj_df$exists, , drop = FALSE]

  required_policies <- c("uniform", "exact", "GI1", "oracle_GI1")
  missing_policy_rows <- list()
  mr <- 0L

  keys <- unique(traj_df[, c("alpha", "case_id")])
  for (i in seq_len(nrow(keys))) {
    a <- keys$alpha[i]
    cid <- keys$case_id[i]
    sub <- traj_df[traj_df$alpha == a & traj_df$case_id == cid, , drop = FALSE]
    have <- unique(sub$policy)
    miss <- setdiff(required_policies, have)
    if (length(miss) > 0) {
      for (m in miss) {
        mr <- mr + 1L
        missing_policy_rows[[mr]] <- data.frame(alpha = a, case_id = cid, missing_policy = m, stringsAsFactors = FALSE)
      }
    }
  }
  missing_policy_df <- if (length(missing_policy_rows) > 0) do.call(rbind, missing_policy_rows) else data.frame()

  base_rows <- traj_df[traj_df$policy == "uniform" & traj_df$exists, , drop = FALSE]
  if (nrow(base_rows) > 0) {
    b_rows <- list()
    bi <- 0L
    for (i in seq_len(nrow(base_rows))) {
      tr <- base_rows[i, , drop = FALSE]
      pass <- NA
      detail <- "ok"
      from_t <- NA_integer_
      to_t <- NA_integer_
      tryCatch({
        meta <- readRDS(tr$path)
        if (is.null(meta$part_files) || length(meta$part_files) == 0L) stop("no trajectory chunk files")
        chunk <- readRDS(meta$part_files[1])
        if (is.null(chunk$action_idx) || nrow(chunk$action_idx) < 1L) stop("missing action_idx in first chunk")

        action_vec <- as.integer(chunk$action_idx[1, ])
        n_total <- length(action_vec)
        n_init <- suppressWarnings(as.integer(tr$n_init[1]))
        k <- suppressWarnings(as.integer(tr$n_settings[1]))
        if (!is.finite(n_init) || !is.finite(k) || k < 1L) stop("missing n_init/n_settings in manifest")

        from_t <- min(n_total, n_init + 1L)
        to_t <- n_total
        if (from_t <= to_t) {
          idx <- seq.int(from_t, to_t)
          expected <- ((idx - n_init - 1L) %% k) + 1L
          observed <- action_vec[idx]
          pass <- isTRUE(all(observed == expected))
          if (!pass) detail <- "observed sequence differs from cyclic permutation rule"
        } else {
          pass <- TRUE
          detail <- "n_total <= n_init (nothing to check after initialization)"
        }
      }, error = function(e) {
        pass <<- FALSE
        detail <<- as.character(e$message)
      })

      bi <- bi + 1L
      b_rows[[bi]] <- data.frame(
        alpha = as.numeric(tr$alpha[1]),
        case_id = as.character(tr$case_id[1]),
        policy = as.character(tr$policy[1]),
        pass_cyclic = as.logical(pass),
        checked_from_t = as.integer(from_t),
        checked_to_t = as.integer(to_t),
        detail = as.character(detail),
        stringsAsFactors = FALSE
      )
    }
    baseline_seq_check <- do.call(rbind, b_rows)
  }
}

## Write machine-friendly outputs
write.csv(eq23_rows, file.path(results_dir, "verification_eq23_summary.csv"), row.names = FALSE)
write.csv(named_plot_exists, file.path(results_dir, "verification_required_plots.csv"), row.names = FALSE)
write.csv(rank_fail, file.path(results_dir, "verification_policy_rank_failures.csv"), row.names = FALSE)
write.csv(traj_missing, file.path(results_dir, "verification_missing_trajectories.csv"), row.names = FALSE)
write.csv(missing_policy_df, file.path(results_dir, "verification_missing_policy_trajectories.csv"), row.names = FALSE)
write.csv(baseline_seq_check, file.path(results_dir, "verification_permutation_baseline_checks.csv"), row.names = FALSE)

## Write markdown report
md <- c(
  "# Verification Report: Case I/II Multinomial Study",
  "",
  sprintf("Results directory: `%s`", results_dir),
  "",
  "## 1) Asymptotic formula alignment (Bures + Frobenius MSE)",
  "Formula checked: E[D(rho_hat, rho)^2] ≈ (1/n) tr( G_D(theta) I(theta)^(-1) ).",
  "Using `eq23_verification_frobenius_bures.csv`.",
  "",
  "| metric | policy | mean abs gap (final) | q90 abs gap (final) | max abs gap (final) | mean abs gap (tail) | q90 abs gap (tail) | max abs gap (tail) |",
  "|---|---|---:|---:|---:|---:|---:|---:|"
)

for (i in seq_len(nrow(eq23_rows))) {
  md <- c(md, sprintf(
    "| %s | %s | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f |",
    eq23_rows$metric[i],
    eq23_rows$policy[i],
    eq23_rows$mean_abs_rel_gap_final[i],
    eq23_rows$q90_abs_rel_gap_final[i],
    eq23_rows$max_abs_rel_gap_final[i],
    eq23_rows$mean_abs_rel_gap_tail[i],
    eq23_rows$q90_abs_rel_gap_tail[i],
    eq23_rows$max_abs_rel_gap_tail[i]
  ))
}

md <- c(md, "", "## 2) Plot generation checks", "")
md <- c(md, sprintf("- Total manifest plots: %d", ifelse(nrow(plot_df) > 0, nrow(plot_df), 0L)))
md <- c(md, sprintf("- Missing manifest plots: %d", nrow(missing_plots)))
md <- c(md, "")
md <- c(md, "Required summary plots:")
for (i in seq_len(nrow(named_plot_exists))) {
  md <- c(md, sprintf("- `%s`: %s", named_plot_exists$filename[i], ifelse(named_plot_exists$exists[i], "OK", "MISSING")))
}

md <- c(md, "", "## 3) Policy ranking checks (Exact/GI1/OracleGI1 vs Permutation)", "")
md <- c(md, sprintf("- Total alpha/case entries: %d", ifelse(nrow(comment_df) > 0, nrow(comment_df), 0L)))
md <- c(md, sprintf("- Rank failures: %d", nrow(rank_fail)))
if (nrow(rank_fail) > 0) {
  md <- c(md, "", "Failures:")
  for (i in seq_len(nrow(rank_fail))) {
    md <- c(md, sprintf("- alpha=%.1f case=%s comment=%s",
                         rank_fail$alpha[i], rank_fail$case_id[i], rank_fail$comment[i]))
  }
} else {
  md <- c(md, "- All cases satisfy Exact<Permutation, GI1<Permutation, OracleGI1<Permutation at final step.")
}

md <- c(md, "", "## 4) Trajectory artifact checks", "")
md <- c(md, sprintf("- Trajectory rows: %d", ifelse(nrow(traj_df) > 0, nrow(traj_df), 0L)))
md <- c(md, sprintf("- Missing trajectory files: %d", nrow(traj_missing)))
md <- c(md, sprintf("- Missing policy trajectory rows: %d", nrow(missing_policy_df)))
if (nrow(baseline_seq_check) > 0) {
  md <- c(md, sprintf("- Permutation cyclic check rows: %d", nrow(baseline_seq_check)))
  md <- c(md, sprintf("- Permutation cyclic failures: %d", sum(!baseline_seq_check$pass_cyclic)))
} else {
  md <- c(md, "- Permutation cyclic check rows: 0 (no baseline trajectories found).")
}

report_file <- file.path(results_dir, "verification_report.md")
writeLines(md, report_file)

cat("Wrote verification outputs:\n")
cat(sprintf("  - %s\n", file.path(results_dir, "verification_eq23_summary.csv")))
cat(sprintf("  - %s\n", file.path(results_dir, "verification_required_plots.csv")))
cat(sprintf("  - %s\n", file.path(results_dir, "verification_policy_rank_failures.csv")))
cat(sprintf("  - %s\n", file.path(results_dir, "verification_missing_trajectories.csv")))
cat(sprintf("  - %s\n", file.path(results_dir, "verification_missing_policy_trajectories.csv")))
cat(sprintf("  - %s\n", file.path(results_dir, "verification_permutation_baseline_checks.csv")))
cat(sprintf("  - %s\n", report_file))
