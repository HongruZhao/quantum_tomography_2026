## ============================================================================
## 00_simulation_parameters.R
## Centralized parameter defaults and profile loading for final simulations.
## ============================================================================

parse_cli_value <- function(prefix, args = commandArgs(trailingOnly = TRUE), default = NULL) {
  hit <- grep(paste0("^", prefix), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", prefix), "", hit[1])
}

as_int_or <- function(x, default) {
  v <- suppressWarnings(as.integer(x))
  if (is.na(v)) default else v
}

as_num_or <- function(x, default) {
  v <- suppressWarnings(as.numeric(x))
  if (is.na(v)) default else v
}

as_logical_or <- function(x, default) {
  if (is.logical(x) && length(x) == 1L && is.finite(as.integer(x))) return(x)
  if (is.null(x)) return(default)
  xs <- tolower(trimws(as.character(x)))
  if (xs %in% c("1", "true", "t", "yes", "y")) return(TRUE)
  if (xs %in% c("0", "false", "f", "no", "n")) return(FALSE)
  default
}

parse_num_vec <- function(x, default = numeric()) {
  vals <- suppressWarnings(as.numeric(strsplit(as.character(x), ",", fixed = TRUE)[[1]]))
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0L) default else unique(vals)
}

parse_int_vec <- function(x, default = integer()) {
  vals <- suppressWarnings(as.integer(strsplit(as.character(x), ",", fixed = TRUE)[[1]]))
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0L) default else unique(vals)
}

resolve_baseline_policy_mode <- function(mode) {
  x <- tolower(trimws(as.character(mode)))
  if (x %in% c("permutation_cyclic", "permutation", "cyclic", "uniform_cyclic")) return("permutation_cyclic")
  if (x %in% c("uniform_random", "random")) return("uniform_random")
  warning(sprintf("Unknown baseline_policy_mode='%s'; fallback to permutation_cyclic.", mode))
  "permutation_cyclic"
}

load_external_config <- function(config_file) {
  if (is.null(config_file) || !nzchar(config_file)) return(list())
  if (!file.exists(config_file)) stop(sprintf("Config file not found: %s", config_file))

  ext <- tolower(tools::file_ext(config_file))
  if (ext == "rds") {
    cfg <- readRDS(config_file)
  } else {
    cfg <- tryCatch(dget(config_file), error = function(e) NULL)
    if (is.null(cfg)) {
      env <- new.env(parent = emptyenv())
      sys.source(config_file, envir = env)
      if (exists("config", envir = env, inherits = FALSE)) {
        cfg <- get("config", envir = env, inherits = FALSE)
      } else {
        stop("Config file must be an .rds list or an R expression/list (dget) or define object 'config'.")
      }
    }
  }

  if (!is.list(cfg)) stop("Loaded config must be a list.")
  cfg
}

merge_config <- function(base, overlay) {
  out <- base
  if (is.null(overlay) || length(overlay) == 0L) return(out)
  nms <- intersect(names(overlay), names(base))
  for (nm in nms) out[[nm]] <- overlay[[nm]]
  out
}

base_final_config <- function() {
  list(
    alphas = c(0.5, 0.9),
    n_rep = 1000L,
    n_total = 200L,
    n_trials = 100L,
    check_every = 1L,

    eta_mle = 1e-4,
    ridge_selection = 0,
    fisher_eps = 0,
    nll_eps_log = 0,
    bures_tol = 0,
    truncate_density = TRUE,
    allow_eta_zero = TRUE,

    solver_selection_mode = "pilot_choose",
    solver_selection_mode_resolved = "cvx_fallback_pgd",
    production_solver = "CVX",
    cvx_solver_preference = "SCS",
    solver = "PGD",
    oracle_solver = "SCS",

    pilot_n_rep = 2L,
    pilot_n_total = 40L,
    pilot_cases = c(1L, 11L),
    pilot_alpha = 0.9,

    n_workers = 12L,
    seed_base = 20260212L,
    pgd_max_iter = 20L,
    pgd_step = NA_real_,
    pgd_tol = NA_real_,
    resume = TRUE,

    save_full_trajectory = TRUE,
    trajectory_chunk_size = 25L,
    save_full_matrices = FALSE,

    plot_1q_start = 10L,
    plot_2q_start = 20L,
    plot_t_end = 200L,

    baseline_policy_mode = "permutation_cyclic",
    uniform_policy_mode = "cyclic",

    output_dir = "",
    plots_dir = "",
    trajectories_dir = "",

    l3_seed = 42L,
    alpha = 0.9
  )
}

profile_overrides <- function(profile = "full") {
  p <- tolower(trimws(as.character(profile)))
  if (p == "full") {
    return(list(
      alphas = c(0.5, 0.9),
      n_rep = 1000L,
      n_total = 200L,
      n_trials = 100L,
      plot_1q_start = 10L,
      plot_2q_start = 20L,
      output_dir = "results_final_caseI_caseII_all15"
    ))
  }
  if (p == "smoke_alpha0p9") {
    return(list(
      alphas = c(0.9),
      n_rep = 2L,
      n_total = 40L,
      n_trials = 100L,
      n_workers = 2L,
      pilot_n_rep = 1L,
      pilot_n_total = 15L,
      save_full_trajectory = TRUE,
      save_full_matrices = FALSE,
      plot_1q_start = 10L,
      plot_2q_start = 20L,
      output_dir = "results_smoke_alpha0p9_all15"
    ))
  }
  if (p == "case11") {
    return(list(
      alpha = 0.9,
      n_rep = 1000L,
      n_total = 1000L,
      n_trials = 100L,
      plot_1q_start = 10L,
      plot_2q_start = 20L,
      output_dir = "results_case11A_bures_final",
      baseline_policy_mode = "permutation_cyclic"
    ))
  }
  warning(sprintf("Unknown config_profile='%s'; using full.", profile))
  profile_overrides("full")
}

coerce_config_types <- function(cfg) {
  int_keys <- c(
    "n_rep", "n_total", "n_trials", "check_every", "pilot_n_rep", "pilot_n_total",
    "n_workers", "seed_base", "pgd_max_iter", "trajectory_chunk_size", "plot_1q_start",
    "plot_2q_start", "plot_t_end", "l3_seed"
  )
  num_keys <- c("eta_mle", "ridge_selection", "fisher_eps", "nll_eps_log", "bures_tol", "pilot_alpha", "pgd_step", "pgd_tol", "alpha")
  logi_keys <- c("truncate_density", "allow_eta_zero", "resume", "save_full_trajectory", "save_full_matrices")

  for (k in int_keys) if (!is.null(cfg[[k]])) cfg[[k]] <- suppressWarnings(as.integer(cfg[[k]]))
  for (k in num_keys) if (!is.null(cfg[[k]])) cfg[[k]] <- suppressWarnings(as.numeric(cfg[[k]]))
  for (k in logi_keys) if (!is.null(cfg[[k]])) cfg[[k]] <- as_logical_or(cfg[[k]], FALSE)

  if (!is.null(cfg$alphas)) cfg$alphas <- parse_num_vec(paste(cfg$alphas, collapse = ","), default = c(0.5, 0.9))
  if (!is.null(cfg$pilot_cases)) cfg$pilot_cases <- parse_int_vec(paste(cfg$pilot_cases, collapse = ","), default = c(1L, 11L))

  cfg$baseline_policy_mode <- resolve_baseline_policy_mode(cfg$baseline_policy_mode %||% cfg$uniform_policy_mode %||% "permutation_cyclic")
  cfg$uniform_policy_mode <- if (identical(cfg$baseline_policy_mode, "permutation_cyclic")) "cyclic" else "random"

  cfg$plot_t_end <- as.integer(cfg$n_total)
  if (cfg$check_every != 1L) cfg$check_every <- 1L

  cfg
}

`%||%` <- function(x, y) if (is.null(x)) y else x

apply_common_cli_overrides <- function(cfg, args = commandArgs(trailingOnly = TRUE)) {
  scalar_int <- c("n_rep", "n_total", "n_trials", "check_every", "n_workers", "seed_base", "pgd_max_iter", "pilot_n_rep", "pilot_n_total", "trajectory_chunk_size", "plot_1q_start", "plot_2q_start", "plot_t_end")
  scalar_num <- c("eta_mle", "ridge_selection", "fisher_eps", "nll_eps_log", "bures_tol", "pilot_alpha", "pgd_step", "pgd_tol", "alpha")
  scalar_chr <- c("solver_selection_mode", "cvx_solver_preference", "solver", "oracle_solver", "output_dir", "baseline_policy_mode", "uniform_policy_mode")
  scalar_logi <- c("resume", "save_full_trajectory", "save_full_matrices", "truncate_density", "allow_eta_zero")

  for (k in scalar_int) {
    v <- parse_cli_value(paste0("--", k, "="), args = args, default = NA_character_)
    if (!is.na(v)) cfg[[k]] <- as_int_or(v, cfg[[k]])
  }
  for (k in scalar_num) {
    v <- parse_cli_value(paste0("--", k, "="), args = args, default = NA_character_)
    if (!is.na(v)) cfg[[k]] <- as_num_or(v, cfg[[k]])
  }
  for (k in scalar_chr) {
    v <- parse_cli_value(paste0("--", k, "="), args = args, default = NA_character_)
    if (!is.na(v) && nzchar(v)) cfg[[k]] <- v
  }
  for (k in scalar_logi) {
    v <- parse_cli_value(paste0("--", k, "="), args = args, default = NA_character_)
    if (!is.na(v)) cfg[[k]] <- as_logical_or(v, cfg[[k]])
  }

  alphas_cli <- parse_cli_value("--alphas=", args = args, default = NA_character_)
  if (!is.na(alphas_cli)) cfg$alphas <- parse_num_vec(alphas_cli, default = cfg$alphas)

  pilot_cases_cli <- parse_cli_value("--pilot_cases=", args = args, default = NA_character_)
  if (!is.na(pilot_cases_cli)) cfg$pilot_cases <- parse_int_vec(pilot_cases_cli, default = cfg$pilot_cases)

  coerce_config_types(cfg)
}

load_simulation_config <- function(profile = "full", config_file = NULL, args = commandArgs(trailingOnly = TRUE)) {
  cfg <- base_final_config()
  cfg <- merge_config(cfg, profile_overrides(profile))
  cfg <- merge_config(cfg, load_external_config(config_file))
  cfg <- coerce_config_types(cfg)
  cfg <- apply_common_cli_overrides(cfg, args = args)
  cfg
}

get_caseI_caseII_config <- function(profile = "full", config_file = NULL, args = commandArgs(trailingOnly = TRUE)) {
  cfg <- load_simulation_config(profile = profile, config_file = config_file, args = args)
  if (!is.null(cfg$plot_t_end) && as.integer(cfg$plot_t_end) != as.integer(cfg$n_total)) {
    warning("plot_t_end is forced to n_total by task rule; overriding provided value.")
  }
  cfg$plot_t_end <- as.integer(cfg$n_total)
  cfg
}

get_case11_config <- function(profile = "case11", config_file = NULL, args = commandArgs(trailingOnly = TRUE)) {
  cfg <- load_simulation_config(profile = profile, config_file = config_file, args = args)
  if (!is.null(cfg$plot_t_end) && as.integer(cfg$plot_t_end) != as.integer(cfg$n_total)) {
    warning("plot_t_end is forced to n_total by task rule; overriding provided value.")
  }
  cfg$plot_t_end <- as.integer(cfg$n_total)
  cfg
}
