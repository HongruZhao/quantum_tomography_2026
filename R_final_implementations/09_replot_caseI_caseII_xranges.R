## ============================================================================
## 09_replot_caseI_caseII_xranges.R
## Replot existing Case I/II outputs with system-specific x-ranges:
##   - 1q: 10..n_total
##   - 2q: 20..n_total
## Includes: Permutation / Exact / GI1 / Oracle GI1.
## ============================================================================

parse_cli_value <- function(prefix, args = commandArgs(trailingOnly = TRUE), default = NULL) {
  hit <- grep(paste0("^", prefix), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", prefix), "", hit[1])
}

alpha_tag <- function(alpha) gsub("\\.", "p", format(alpha, nsmall = 1, trim = TRUE))

policy_style <- function(policy) {
  switch(policy,
    "uniform" = list(label = "Permutation", col = "gray40", lty = 1),
    "exact" = list(label = "Exact", col = "blue", lty = 1),
    "GI1" = list(label = "GI1", col = "red", lty = 2),
    "oracle_GI1" = list(label = "Oracle GI1", col = "darkgreen", lty = 4),
    list(label = policy, col = "black", lty = 1)
  )
}

policy_order <- function(case_result) {
  desired <- c("uniform", "exact", "GI1", "oracle_GI1")
  desired[desired %in% names(case_result$results)]
}

x_idx_for_system <- function(system, n_total,
                             one_start = plot_1q_start,
                             two_start = plot_2q_start,
                             x_end_cap = 500L) {
  x_start <- if (identical(system, "1q")) one_start else two_start
  x_end <- as.integer(n_total)
  if (x_start > x_end) x_start <- 1L
  seq.int(x_start, x_end)
}

plot_case_scaled_proxy <- function(case_result, out_file) {
  n_idx <- x_idx_for_system(case_result$system, case_result$config$n_total)
  n_seq <- n_idx
  p_ord <- policy_order(case_result)

  vals <- unlist(lapply(p_ord, function(p) case_result$results[[p]]$proxy$mean[n_idx]), use.names = FALSE)
  vals <- c(vals, case_result$oracle_limit)
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) vals <- c(0, 1)
  y_lim <- range(vals)
  span <- y_lim[2] - y_lim[1]
  pad <- if (span > 0) 0.05 * span else max(1e-6, abs(y_lim[1]) * 0.05 + 1e-6)
  y_lim <- c(max(0, y_lim[1] - pad), y_lim[2] + pad)

  png(out_file, width = 1100, height = 760, res = 120)
  par(mar = c(5, 5, 4, 2))

  first <- p_ord[1]
  plot(n_seq, case_result$results[[first]]$proxy$mean[n_idx], type = "l",
       col = policy_style(first)$col, lty = policy_style(first)$lty, lwd = 2.5,
       xlab = "Adaptive step (t)", ylab = expression(n[eff] %.% "Proxy risk"),
       main = sprintf("Case %02d | alpha=%.1f | %s | %s",
                      as.integer(case_result$case_num),
                      as.numeric(case_result$alpha),
                      case_result$case_id,
                      case_result$loss),
       ylim = y_lim)

  if (length(p_ord) > 1L) {
    for (p in p_ord[-1]) {
      st <- policy_style(p)
      lines(n_seq, case_result$results[[p]]$proxy$mean[n_idx], col = st$col, lty = st$lty, lwd = 2.5)
    }
  }

  abline(h = case_result$oracle_limit, col = "black", lwd = 2, lty = 3)
  grid(lty = "dotted", col = "gray85")

  legend("topright",
         legend = c(vapply(p_ord, function(p) policy_style(p)$label, character(1)), "Oracle limit"),
         col = c(vapply(p_ord, function(p) policy_style(p)$col, character(1)), "black"),
         lty = c(vapply(p_ord, function(p) policy_style(p)$lty, numeric(1)), 3),
         lwd = c(rep(2.5, length(p_ord)), 2),
         bg = "white")

  mtext(sprintf("Window: t=%d..%d | %s", min(n_idx), max(n_idx), case_result$comment),
        side = 1, line = 3.5, cex = 0.9)

  dev.off()
}

plot_alpha_relative_panel <- function(alpha_results, out_file) {
  ordered <- alpha_results[order(vapply(alpha_results, function(x) as.numeric(x$case_num), numeric(1)))]

  png(out_file, width = 2500, height = 2000, res = 140)
  par(mfrow = c(5, 3), mar = c(4, 4, 3, 1), oma = c(2, 2, 3, 1))

  for (cr in ordered) {
    idx <- x_idx_for_system(cr$system, cr$config$n_total)
    n_seq <- idx
    oracle <- cr$oracle_limit
    p_ord <- policy_order(cr)

    rel_map <- lapply(p_ord, function(p) cr$results[[p]]$proxy$mean[idx] / oracle)
    names(rel_map) <- p_ord

    y_all <- unlist(rel_map, use.names = FALSE)
    y_all <- y_all[is.finite(y_all)]
    if (length(y_all) == 0) y_all <- c(1, 1.2)
    y_max <- max(y_all, na.rm = TRUE)
    pad <- if (is.finite(y_max) && y_max > 1) 0.08 * (y_max - 1) else 0.1
    y_lim <- c(1, y_max + max(pad, 0.05))

    first <- p_ord[1]
    plot(n_seq, rel_map[[first]], type = "l", col = policy_style(first)$col, lwd = 1.8, lty = policy_style(first)$lty,
         xlab = "t", ylab = expression((n[eff] %.% "Proxy risk") / "Oracle limit"),
         ylim = y_lim,
         main = sprintf("%02d %s | %s\n%s",
                        as.integer(cr$case_num), cr$library, cr$loss,
                        ifelse(isTRUE(cr$exact_better) && isTRUE(cr$gi1_better) && isTRUE(cr$oracle_gi1_better),
                               "Adaptive better", "Check ranking")),
         cex.main = 0.95)

    if (length(p_ord) > 1L) {
      for (p in p_ord[-1]) {
        st <- policy_style(p)
        lines(n_seq, rel_map[[p]], col = st$col, lwd = 1.8, lty = st$lty)
      }
    }

    abline(h = 1, col = "black", lty = 3, lwd = 1.4)
    grid(lty = "dotted", col = "gray85")
  }

  n_total_lbl <- as.integer(ordered[[1]]$config$n_total)
  mtext(sprintf("Relative Oracle Panel (alpha=%.1f) | 1q:10..%d, 2q:20..%d", ordered[[1]]$alpha, n_total_lbl, n_total_lbl),
        outer = TRUE, cex = 1.45, line = 0.7)
  par(fig = c(0, 1, 0, 0.04), new = TRUE, mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", horiz = TRUE,
         legend = c("Permutation", "Exact", "GI1", "Oracle GI1", "Oracle ratio = 1"),
         col = c("gray40", "blue", "red", "darkgreen", "black"),
         lty = c(1, 1, 2, 4, 3), lwd = c(2, 2, 2, 2, 1.5), bty = "n")

  dev.off()
}

plot_alpha_metric_panel <- function(alpha_results, loss_target,
                                   metric_slot,
                                   oracle_getter,
                                   ylab,
                                   title_text,
                                   out_file) {
  keep <- alpha_results[vapply(alpha_results, function(x) identical(x$loss, loss_target), logical(1))]
  keep <- keep[order(vapply(keep, function(x) as.numeric(x$case_num), numeric(1)))]
  if (length(keep) == 0) return(invisible(NULL))

  png(out_file, width = 2100, height = 1700, res = 140)
  par(mfrow = c(3, 2), mar = c(4, 4, 3, 1), oma = c(2, 2, 3, 1))

  for (cr in keep) {
    idx <- x_idx_for_system(cr$system, cr$config$n_total)
    n_seq <- idx
    p_ord <- policy_order(cr)

    curve_map <- lapply(p_ord, function(p) cr$results[[p]][[metric_slot]]$mean[idx])
    names(curve_map) <- p_ord
    oracle_line <- oracle_getter(cr)

    y_all <- c(unlist(curve_map, use.names = FALSE), oracle_line)
    y_all <- y_all[is.finite(y_all)]
    if (length(y_all) == 0) y_all <- c(0, 1)
    y_lim <- range(y_all)
    span <- y_lim[2] - y_lim[1]
    pad <- if (span > 0) 0.07 * span else 0.1
    y_lim <- c(max(0, y_lim[1] - pad), y_lim[2] + pad)

    first <- p_ord[1]
    plot(n_seq, curve_map[[first]], type = "l", col = policy_style(first)$col, lty = policy_style(first)$lty, lwd = 2,
         xlab = "Adaptive step (t)", ylab = ylab,
         ylim = y_lim,
         main = sprintf("%02d %s | %s\n%s",
                        as.integer(cr$case_num), cr$library, cr$loss, cr$comment),
         cex.main = 0.9)

    if (length(p_ord) > 1L) {
      for (p in p_ord[-1]) {
        st <- policy_style(p)
        lines(n_seq, curve_map[[p]], col = st$col, lty = st$lty, lwd = 2)
      }
    }

    abline(h = oracle_line, col = "black", lwd = 1.6, lty = 3)
    grid(lty = "dotted", col = "gray85")
  }

  n_total_lbl <- as.integer(keep[[1]]$config$n_total)
  mtext(sprintf("%s | 1q:10..%d, 2q:20..%d", title_text, n_total_lbl, n_total_lbl), outer = TRUE, cex = 1.35, line = 0.7)
  par(fig = c(0, 1, 0, 0.04), new = TRUE, mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", horiz = TRUE,
         legend = c("Permutation", "Exact", "GI1", "Oracle GI1", "Oracle limit"),
         col = c("gray40", "blue", "red", "darkgreen", "black"),
         lty = c(1, 1, 2, 4, 3), lwd = c(2, 2, 2, 2, 1.6), bty = "n")

  dev.off()
}

script_dir <- tryCatch({
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  NULL
}, error = function(e) NULL)
if (is.null(script_dir) || !nzchar(script_dir)) script_dir <- getwd()
if (!file.exists(file.path(script_dir, "00_simulation_parameters.R"))) {
  candidate <- file.path(script_dir, "R_final_implementations")
  if (file.exists(file.path(candidate, "00_simulation_parameters.R"))) script_dir <- candidate
}
setwd(script_dir)
source(file.path(script_dir, "00_simulation_parameters.R"))

cfg_replot <- load_simulation_config(
  profile = parse_cli_value("--config_profile=", default = "full"),
  config_file = parse_cli_value("--config_file=", default = NULL),
  args = commandArgs(trailingOnly = TRUE)
)
plot_1q_start <- as.integer(cfg_replot$plot_1q_start)
plot_2q_start <- as.integer(cfg_replot$plot_2q_start)

results_dir <- parse_cli_value(
  "--results_dir=",
  default = cfg_replot$output_dir %||% "/Users/cyberslave/GitHub/quantum_tomography_2026/R_final_implementations/results_final_caseI_caseII_all15"
)

full_rds <- file.path(results_dir, "full_results_caseI_caseII.rds")
if (!file.exists(full_rds)) stop(sprintf("Missing: %s", full_rds))

all_results <- readRDS(full_rds)
if (!is.list(all_results) || length(all_results) == 0) stop("Invalid full results object.")

plots_dir <- file.path(results_dir, "plots")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

for (nm in names(all_results)) {
  cr <- all_results[[nm]]
  a_tag <- alpha_tag(cr$alpha)
  out_file <- file.path(plots_dir, sprintf("alpha_%s_case%02d_%s_scaled_proxy.png",
                                           a_tag, as.integer(cr$case_num), cr$case_id))
  plot_case_scaled_proxy(cr, out_file)
}

alpha_vals <- sort(unique(vapply(all_results, function(x) as.numeric(x$alpha), numeric(1))))
for (alpha in alpha_vals) {
  a_tag <- alpha_tag(alpha)
  group <- all_results[vapply(all_results, function(x) isTRUE(all.equal(as.numeric(x$alpha), alpha)), logical(1))]

  rel_file <- file.path(plots_dir, sprintf("alpha_%s_all15_relative_oracle_proxy.png", a_tag))
  plot_alpha_relative_panel(group, rel_file)

  fro_file <- file.path(plots_dir, sprintf("alpha_%s_frobenius_mse_vs_oracle.png", a_tag))
  plot_alpha_metric_panel(
    group,
    loss_target = "frobenius",
    metric_slot = "mse",
    oracle_getter = function(cr) cr$oracle_mse_limit,
    ylab = expression(n[eff] %.% "Frobenius MSE"),
    title_text = sprintf("Frobenius MSE vs Oracle (alpha=%.1f)", alpha),
    out_file = fro_file
  )

  bures_file <- file.path(plots_dir, sprintf("alpha_%s_bures_sq_vs_oracle.png", a_tag))
  plot_alpha_metric_panel(
    group,
    loss_target = "bures",
    metric_slot = "bures_sq",
    oracle_getter = function(cr) cr$oracle_limit,
    ylab = expression(n[eff] %.% "Bures"^2),
    title_text = sprintf("Bures Squared Loss vs Oracle (alpha=%.1f)", alpha),
    out_file = bures_file
  )
}

if (length(all_results) > 0) {
  n_total_lbl <- as.integer(all_results[[1]]$config$n_total)
  cat(sprintf("Re-generated plots with windows: 1q=10..%d, 2q=20..%d\n", n_total_lbl, n_total_lbl))
} else {
  cat("Re-generated plots.\n")
}
cat(sprintf("Plots dir: %s\n", plots_dir))
