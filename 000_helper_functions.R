# Helper functions for the species accumulation / asymptotic models.
# Adapted from Zizka et al. (2026) "How many Bromeliads are there?"
# https://doi.org/10.5281/zenodo.21789730
# Required packages: dplyr, tidyr, purrr, broom, minpack.lm, MASS

# --- Helpers ---------------------------------------------------------------

# map scaled time t (1,2,...) to calendar year, given the first year in the series
t_to_year <- function(t, year_min) { year_min + (t - 1) }

# curve predictors
pred_mm   <- function(x, fit) { cf <- coef(fit); (cf["a"] * x) / (cf["b"] + x) }
pred_logi <- function(x, fit) { cf <- coef(fit); cf["a"] / (1 + exp(-cf["k"] * (x - cf["t0"]))) }

# Percentile bootstrap CI for 'a' using MVN draws + physical filters
get_a_ci_boot <- function(fit, level = 0.95,
                          k_bounds = c(1e-4, 0.5),
                          t_bounds = NULL,    # c(tmin, tmax) to bound t0
                          ymin_cap = NULL,    # cap multiplier base (e.g. fitted max)
                          B = 1000) {

  if (is.null(fit)) return(tibble::tibble(a_estimate = NA_real_, a_lwr = NA_real_, a_upr = NA_real_))
  cf <- coef(fit)
  vc <- tryCatch(vcov(fit), error = function(e) NULL)
  if (is.null(vc) || !"a" %in% names(cf)) return(tibble::tibble(a_estimate = NA_real_, a_lwr = NA_real_, a_upr = NA_real_))

  # Draw parameters
  draws <- tryCatch(MASS::mvrnorm(B, mu = cf, Sigma = vc), error = function(e) NULL)
  if (is.null(draws)) return(tibble::tibble(a_estimate = unname(cf["a"]), a_lwr = NA_real_, a_upr = NA_real_))
  colnames(draws) <- names(cf)

  # Filters
  ok <- rep(TRUE, nrow(draws))
  if (all(c("k","t0") %in% colnames(draws))) {
    ok <- ok & draws[, "a"] > 0 & draws[, "k"] > k_bounds[1] & draws[, "k"] < k_bounds[2]
    if (!is.null(t_bounds) && length(t_bounds) == 2) {
      tmin <- t_bounds[1]; tmax <- t_bounds[2]; tpad <- 0.25 * (tmax - tmin)
      ok <- ok & draws[, "t0"] > (tmin - tpad) & draws[, "t0"] < (tmax + tpad)
    }
  } else {
    ok <- ok & draws[, "a"] > 0
  }
  draws <- draws[ok, , drop = FALSE]
  if (nrow(draws) < 100) return(tibble::tibble(a_estimate = unname(cf["a"]), a_lwr = NA_real_, a_upr = NA_real_))

  a_draws <- draws[, "a"]
  # Enforce 'a' not below observed max (since fit used bounds). Small tolerance:
  a_draws <- a_draws[is.finite(a_draws) & a_draws > 0]

  alpha <- 1 - level
  tibble::tibble(
    a_estimate = unname(cf["a"]),
    a_lwr      = unname(quantile(a_draws, probs = alpha/2, na.rm = TRUE)),
    a_upr      = unname(quantile(a_draws, probs = 1 - alpha/2, na.rm = TRUE))
  )
}


# Parametric bootstrap CI bands with strict filters; returns year and t
predict_curve_band <- function(fit,
                               t_min, t_max, year_min,
                               extra_years = 150,        # extrapolate
                               model = c("MichaelisMenten","Logistic"),
                               level = 0.95, B = 800,
                               max_mult = 2,             # cap relative to fitted max
                               k_bounds = c(1e-4, 0.5),  # logistic k bounds
                               t0_pad_frac = 0.25) {     # allow t0 slightly outside window
  model <- match.arg(model)
  if (is.null(fit)) return(NULL)

  cf <- coef(fit)
  vc <- tryCatch(vcov(fit), error = function(e) NULL)
  if (is.null(vc)) return(NULL)

  # Time grid: observed window … plus +extra_years
  t_seq <- seq(t_min, t_max + extra_years, length.out = 300)
  year_seq <- t_to_year(t_seq, year_min)

  f_mm   <- function(x, par) (par["a"] * x) / (par["b"] + x)
  f_logi <- function(x, par) par["a"] / (1 + exp(-par["k"] * (x - par["t0"])))

  # Cap based on fitted curve within observed window
  y_obs_max <- max(fitted(fit), na.rm = TRUE)
  y_cap     <- y_obs_max * max_mult

  # Draw parameter vectors
  draws <- tryCatch(MASS::mvrnorm(B, mu = cf, Sigma = vc), error = function(e) NULL)
  if (is.null(draws)) return(NULL)
  colnames(draws) <- names(cf)

  # Param filters
  t_pad <- (t_max - t_min) * t0_pad_frac
  ok <- rep(TRUE, nrow(draws))
  if (model == "MichaelisMenten") {
    ok <- ok & draws[,"a"] > 0 & draws[,"b"] > 0 & draws[,"b"] < 10 * (t_max - t_min + 1)
  } else {
    ok <- ok & draws[,"a"] > 0 &
      draws[,"k"] > k_bounds[1] & draws[,"k"] < k_bounds[2] &
      draws[,"t0"] > (t_min - t_pad) & draws[,"t0"] < (t_max + t_pad)
  }
  draws <- draws[ok, , drop = FALSE]
  if (nrow(draws) < 80) return(NULL)

  # Predict safely over the extended grid
  pred_mat <- matrix(NA_real_, nrow = nrow(draws), ncol = length(t_seq))
  for (i in seq_len(nrow(draws))) {
    pr <- draws[i, ]
    yi <- if (model == "MichaelisMenten") f_mm(t_seq, pr) else f_logi(t_seq, pr)
    if (any(!is.finite(yi))) next
    if (any(yi < 0) || any(yi > y_cap)) next
    if (any(diff(yi) < -1e-8)) next
    pred_mat[i, ] <- yi
  }
  keep <- rowSums(is.finite(pred_mat)) == ncol(pred_mat)
  pred_mat <- pred_mat[keep, , drop = FALSE]
  if (nrow(pred_mat) < 80) return(NULL)

  # Central curve (extended grid)
  y_hat <- if (model == "MichaelisMenten") f_mm(t_seq, cf) else f_logi(t_seq, cf)
  y_hat <- pmin(pmax(y_hat, 0), y_cap)

  alpha <- 1 - level
  lwr <- apply(pred_mat, 2, quantile, probs = alpha/2, na.rm = TRUE)
  upr <- apply(pred_mat, 2, quantile, probs = 1 - alpha/2, na.rm = TRUE)

  tibble::tibble(t = t_seq, year = year_seq, y_hat = y_hat, lwr = lwr, upr = upr, model = model)
}


# --- Quasi-Poisson GLM cross-validation (Bebber et al. 2007) ---------------
# Models annual counts of new descriptions (not cumulative) via a
# quasi-Poisson GLM with cubic polynomial in t. Does not assume a sigmoid
# shape and handles overdispersion from monograph bursts. K is estimated
# by integrating projected annual counts 300 years into the future.
# If QPois K ≈ Logistic K → logistic estimate is credible.

fit_quasipoisson <- function(df_annual, t_col = "t", y_col = "y",
                             project_years = 300) {
  df <- df_annual %>%
    mutate(t_num = suppressWarnings(as.numeric(.data[[t_col]])),
           y_ann = suppressWarnings(as.numeric(.data[[y_col]]))) %>%
    filter(!is.na(t_num), !is.na(y_ann))
  if (nrow(df) < 6) return(tibble::tibble(K_qp = NA_real_, completeness_qp = NA_real_))
  s_obs  <- max(cumsum(df$y_ann), na.rm = TRUE)
  qp_fit <- tryCatch(
    glm(y_ann ~ poly(t_num, 3), data = df, family = quasipoisson(link = "log")),
    error = function(e) NULL
  )
  if (is.null(qp_fit)) return(tibble::tibble(K_qp = NA_real_, completeness_qp = NA_real_))
  t_future <- seq(max(df$t_num) + 1, max(df$t_num) + project_years)
  pred_ann <- pmax(predict(qp_fit,
                           newdata = data.frame(t_num = t_future),
                           type = "response"), 0)
  K_qp <- s_obs + sum(pred_ann)
  tibble::tibble(K_qp = round(K_qp, 0), completeness_qp = round(s_obs / K_qp, 3))
}

# --- Model fitting ----------------------------------------------------------

fit_models <- function(df, t_col = "t", cum_col = "cum_y") {
  df <- df %>%
    mutate(t_num = suppressWarnings(as.numeric(.data[[t_col]])),
           y_num = suppressWarnings(as.numeric(.data[[cum_col]]))) %>%
    filter(!is.na(t_num), !is.na(y_num)) %>%
    arrange(t_num)

  if (nrow(df) < 6) {
    return(tibble::tibble(model = character(), fit = list(), params = list(), glance = list()))
  }

  ymax   <- max(df$y_num, na.rm = TRUE)          # observed max (must be ≤ asymptote)
  trange <- diff(range(df$t_num, na.rm = TRUE))
  tmin   <- min(df$t_num, na.rm = TRUE)
  tmax   <- max(df$t_num, na.rm = TRUE)

  # Starts
  start_a_mm  <- max(ymax * 1.05, ymax + 1)
  start_b_mm  <- max(1, round(median(df$t_num)))
  start_a_lo  <- max(ymax * 1.05, ymax + 1)
  start_k_lo  <- 0.05
  start_t0_lo <- median(df$t_num)

  # Bounds (named vectors) — keep them fairly generous but sane
  lower_mm <- c(a = ymax * 1.001, b = 1e-6)
  upper_mm <- c(a = ymax * 100,   b = max(10, 50 * (trange + 1)))

  lower_lo <- c(a = ymax * 1.001, k = 1e-5, t0 = tmin - 0.5 * trange)
  upper_lo <- c(a = ymax * 100,   k = 1.0,  t0 = tmax + 2.0 * trange)

  mm_fit <- tryCatch(
    minpack.lm::nlsLM(y_num ~ (a * t_num) / (b + t_num),
                      data = df,
                      start = list(a = start_a_mm, b = start_b_mm),
                      lower = lower_mm, upper = upper_mm,
                      control = minpack.lm::nls.lm.control(maxiter = 500)
    ),
    error = function(e) NULL
  )

  logi_fit <- tryCatch(
    minpack.lm::nlsLM(y_num ~ a / (1 + exp(-k * (t_num - t0))),
                      data = df,
                      start = list(a = start_a_lo, k = start_k_lo, t0 = start_t0_lo),
                      lower = lower_lo, upper = upper_lo,
                      control = minpack.lm::nls.lm.control(maxiter = 500)
    ),
    error = function(e) NULL
  )

  tibble::tibble(
    model = c("MichaelisMenten", "Logistic"),
    fit   = list(mm_fit, logi_fit)
  ) %>%
    mutate(params = purrr::map(fit, ~ if (!is.null(.x)) broom::tidy(.x) else NULL),
           glance = purrr::map(fit, ~ if (!is.null(.x)) broom::glance(.x) else NULL)) %>%
    filter(!purrr::map_lgl(fit, is.null))
}

# Bootstrap-filter parameters per group.
# Defaults for all groups; add stricter filters for specific groups if the
# CI bands explode (in Bromeliaceae this was needed for Tillandsioideae).
get_band_params <- function(group, strict_groups = character(0)) {
  if (group %in% strict_groups) {
    list(k_bounds = c(1e-3, 0.2), max_mult = 1.5)
  } else {
    list(k_bounds = c(1e-4, 0.5), max_mult = 2.0)
  }
}

# --- Cumulative series -------------------------------------------------------
# Counts first descriptions per year and builds the cumulative series.
# `group` (optional, character) = column to split by, e.g. "region" or "tax_genus".
# t = years since the first description in that series (t = 1 in the first year).
# Only years with at least one description are kept (as in the Bromeliaceae pipeline).
build_cumulative <- function(dat, group = NULL, year_col = "year_first_description") {
  d <- dat %>%
    dplyr::mutate(year = suppressWarnings(as.numeric(as.character(.data[[year_col]])))) %>%
    dplyr::filter(!is.na(year))
  if (!is.null(group)) d <- d %>% dplyr::filter(!is.na(.data[[group]]))
  d %>%
    dplyr::count(dplyr::across(dplyr::all_of(c(group, "year"))), name = "y") %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group))) %>%
    dplyr::arrange(year, .by_group = TRUE) %>%
    dplyr::mutate(cum_y = cumsum(y),
                  t = year - min(year) + 1) %>%
    dplyr::ungroup()
}
