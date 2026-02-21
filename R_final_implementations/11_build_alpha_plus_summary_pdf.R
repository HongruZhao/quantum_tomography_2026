## ============================================================================
## 11_build_alpha_plus_summary_pdf.R
## Build 3-page alpha summary PDF from existing PNG plots.
## Pages:
##   1) alpha_*_all15_relative_oracle_proxy.png
##   2) alpha_*_frobenius_mse_vs_oracle.png
##   3) alpha_*_bures_sq_vs_oracle.png
## ============================================================================

parse_cli_value <- function(prefix, args = commandArgs(trailingOnly = TRUE), default = NULL) {
  hit <- grep(paste0("^", prefix), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", prefix), "", hit[1])
}

alpha_tag <- function(alpha) gsub("\\.", "p", format(alpha, nsmall = 1, trim = TRUE))

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

cfg <- load_simulation_config(
  profile = parse_cli_value("--config_profile=", default = "full"),
  config_file = parse_cli_value("--config_file=", default = NULL),
  args = commandArgs(trailingOnly = TRUE)
)

results_dir <- parse_cli_value("--results_dir=", default = cfg$output_dir %||% "results_final_caseI_caseII_all15")
alpha <- suppressWarnings(as.numeric(parse_cli_value("--alpha=", default = "0.9")))
if (!is.finite(alpha)) alpha <- 0.9
a_tag <- alpha_tag(alpha)

if (!dir.exists(results_dir)) stop(sprintf("Missing results_dir: %s", results_dir))
plots_dir <- file.path(results_dir, "plots")
if (!dir.exists(plots_dir)) stop(sprintf("Missing plots_dir: %s", plots_dir))

png_files <- c(
  file.path(plots_dir, sprintf("alpha_%s_all15_relative_oracle_proxy.png", a_tag)),
  file.path(plots_dir, sprintf("alpha_%s_frobenius_mse_vs_oracle.png", a_tag)),
  file.path(plots_dir, sprintf("alpha_%s_bures_sq_vs_oracle.png", a_tag))
)

for (p in png_files) {
  if (!file.exists(p)) stop(sprintf("Missing required PNG: %s", p))
}

out_pdf <- parse_cli_value(
  "--out_pdf=",
  default = file.path(results_dir, sprintf("alpha_%s_all15_plus_summary.pdf", a_tag))
)

join_bin <- "/System/Library/Automator/Combine PDF Pages.action/Contents/MacOS/join"
if (!file.exists(join_bin)) stop(sprintf("Missing join binary: %s", join_bin))

tmp_dir <- file.path(results_dir, sprintf(".pdf_pages_tmp_alpha_%s", a_tag))
if (dir.exists(tmp_dir)) unlink(tmp_dir, recursive = TRUE, force = TRUE)
dir.create(tmp_dir, recursive = TRUE)

pdf_pages <- character(length(png_files))
for (i in seq_along(png_files)) {
  in_png <- png_files[i]
  out_page <- file.path(tmp_dir, sprintf("%02d_%s.pdf", i, sub("\\.png$", "", basename(in_png))))
  status <- system2("sips", c("-s", "format", "pdf", in_png, "--out", out_page), stdout = FALSE, stderr = FALSE)
  if (!identical(status, 0L) || !file.exists(out_page)) {
    stop(sprintf("Failed PNG->PDF conversion: %s", in_png))
  }
  pdf_pages[i] <- out_page
}

join_status <- system2(join_bin, c("--output", out_pdf, pdf_pages), stdout = FALSE, stderr = FALSE)
if (!identical(join_status, 0L) || !file.exists(out_pdf)) {
  stop(sprintf("Failed to build summary PDF: %s", out_pdf))
}

unlink(tmp_dir, recursive = TRUE, force = TRUE)

cat("Built alpha summary PDF:\n")
cat(sprintf("  - %s\n", out_pdf))
cat("Using pages:\n")
for (p in png_files) cat(sprintf("  - %s\n", p))
