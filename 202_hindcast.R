# 202 - Hindcast: which model predicts future descriptions best?
#
# Fits the models with the descriptions up to a cutoff year only, predicts how
# many species should have been described by the last year in the data, and
# compares that with what was actually described. Repeated for several cutoffs
# so that the result does not depend on one particular year.
#
# Models compared
#   - Logistic and Michaelis-Menten on the cumulative series (as in 201)
#   - Quasi-Poisson GLM on annual counts, two versions:
#       qp_nonzero : only years with >= 1 description (as in 201 / 000_helper)
#       qp_zerofill: all calendar years, years without descriptions = 0
#                    (as in Bebber et al. 2007)
#
# For each model and cutoff we report
#   - predicted cumulative number of species in the last year vs observed
#   - error (predicted - observed new species after the cutoff)
#   - estimated total (asymptote / K) with the data up to the cutoff, to see
#     how stable each model's estimate is as data accumulate
#
# Input : output/velloziaceae_data_for_analyses.rda   (from 101)
# Output: output/hindcast_results.rda, output/hindcast_results.xlsx,
#         output/diagnostic_hindcast_*.png

library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(minpack.lm)
library(writexl)
library(ggplot2)

source("000_helper_functions.R")
set.seed(4657)

load("output/velloziaceae_data_for_analyses.rda")   # object: outp

cutoffs     <- c(1980, 1990, 2000, 2005)  # 2000 = main hindcast, others = robustness
min_species <- 20                         # same threshold as 201
min_years_train <- 8                      # need some years with descriptions before the cutoff

# --- 1. Series --------------------------------------------------------------
peryear        <- build_cumulative(outp) %>% mutate(group = "Velloziaceae", level = "overall")
peryear_region <- build_cumulative(outp, "region") %>% rename(group = region) %>% mutate(level = "region")
peryear_genus  <- build_cumulative(outp, "tax_genus") %>% rename(group = tax_genus) %>% mutate(level = "genus")

series_all <- bind_rows(peryear, peryear_region, peryear_genus) %>%
  group_by(level, group) %>%
  filter(max(cum_y) >= min_species) %>%
  ungroup()

year_end <- max(series_all$year)

# --- 2. Quasi-Poisson with prediction over an arbitrary set of years --------
# Fits annual counts ~ poly(year, 3). With zero_fill = TRUE, years without any
# description between the first description and the cutoff are added as 0.
qp_fit_predict <- function(train, years_pred, zero_fill = TRUE, project_years = 300) {
  d <- train %>% select(year, y)
  if (zero_fill) {
    d <- tibble(year = seq(min(d$year), max(d$year))) %>%
      left_join(d, by = "year") %>%
      mutate(y = coalesce(y, 0L))
  }
  fit <- tryCatch(glm(y ~ poly(year, 3), data = d, family = quasipoisson(link = "log")),
                  error = function(e) NULL)
  if (is.null(fit)) return(list(new = NA_real_, K = NA_real_))
  pr  <- function(yrs) pmax(predict(fit, newdata = data.frame(year = yrs), type = "response"), 0)
  s   <- sum(train$y)
  new <- sum(pr(years_pred))
  K   <- s + sum(pr(seq(max(train$year) + 1, max(train$year) + project_years)))
  list(new = new, K = K)
}

# --- 3. One hindcast ---------------------------------------------------------
hindcast_one <- function(series, cutoff) {
  train <- series %>% filter(year <= cutoff)
  if (nrow(train) < min_years_train) return(NULL)

  year_min <- min(series$year)
  t_cut    <- cutoff - year_min + 1
  t_end    <- year_end - year_min + 1
  obs_cut  <- max(train$cum_y)
  obs_end  <- max(series$cum_y)

  # logistic + MM
  fits <- fit_models(train, t_col = "t", cum_col = "cum_y")
  curve_res <- fits %>%
    mutate(pred_end = map2_dbl(fit, model, function(f, m) {
             if (m == "Logistic") unname(pred_logi(t_end, f)) else unname(pred_mm(t_end, f))
           }),
           # prediction of new species = curve increase between cutoff and end,
           # added to the observed count at the cutoff (the curve need not pass
           # exactly through the last observed point)
           fit_cut  = map2_dbl(fit, model, function(f, m) {
             if (m == "Logistic") unname(pred_logi(t_cut, f)) else unname(pred_mm(t_cut, f))
           }),
           pred_new = pred_end - fit_cut,
           K        = map_dbl(fit, ~ unname(coef(.x)["a"]))) %>%
    select(model, pred_new, K)

  # quasi-Poisson
  yrs_pred <- seq(cutoff + 1, year_end)
  qp0 <- qp_fit_predict(train, yrs_pred, zero_fill = FALSE)
  qp1 <- qp_fit_predict(train, yrs_pred, zero_fill = TRUE)
  qp_res <- tibble(model    = c("QuasiPoisson_nonzero", "QuasiPoisson_zerofill"),
                   pred_new = c(qp0$new, qp1$new),
                   K        = c(qp0$K, qp1$K))

  bind_rows(curve_res, qp_res) %>%
    mutate(cutoff        = cutoff,
           observed_cut  = obs_cut,
           observed_end  = obs_end,
           observed_new  = obs_end - obs_cut,
           predicted_end = observed_cut + pred_new,
           error_new     = pred_new - observed_new,
           rel_error_new = error_new / observed_new,
           # logistic / MM asymptote stuck at the fitting bounds in fit_models()
           # (a = 1.001 x observed or a = 100 x observed): not identified
           K_at_bound    = model %in% c("Logistic", "MichaelisMenten") &
             (K <= observed_cut * 1.002 | K >= observed_cut * 99))
}

hindcast <- series_all %>%
  group_by(level, group) %>%
  group_modify(~ map_dfr(cutoffs, function(cu) hindcast_one(.x, cu))) %>%
  ungroup() %>%
  mutate(model = factor(model, levels = c("Logistic", "MichaelisMenten",
                                          "QuasiPoisson_nonzero", "QuasiPoisson_zerofill")))

# --- 4. Current estimates (all data) for the stability plot -----------------
current <- series_all %>%
  group_by(level, group) %>%
  group_modify(~ {
    f  <- fit_models(.x, t_col = "t", cum_col = "cum_y") %>%
      transmute(model, K = map_dbl(fit, ~ unname(coef(.x)["a"])))
    q0 <- qp_fit_predict(.x, integer(0), zero_fill = FALSE)$K
    q1 <- qp_fit_predict(.x, integer(0), zero_fill = TRUE)$K
    bind_rows(f, tibble(model = c("QuasiPoisson_nonzero", "QuasiPoisson_zerofill"),
                        K = c(q0, q1)))
  }) %>%
  ungroup() %>%
  mutate(cutoff = year_end,
         model  = factor(model, levels = levels(hindcast$model)))

# --- 5. Summary ---------------------------------------------------------------
# Mean absolute error of new species across cutoffs, per group and model
summary_tab <- hindcast %>%
  group_by(level, group, model) %>%
  summarise(n_cutoffs        = sum(is.finite(error_new)),
            MAE_new          = mean(abs(error_new), na.rm = TRUE),
            mean_rel_error   = mean(rel_error_new, na.rm = TRUE),
            K_range_min      = min(K, na.rm = TRUE),
            K_range_max      = max(K, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(level, group, MAE_new)

print(as.data.frame(hindcast %>% filter(cutoff == 2000) %>%
                      select(level, group, model, observed_new, pred_new, error_new, K) %>%
                      mutate(across(where(is.numeric), ~ round(.x, 1)))))
print(as.data.frame(summary_tab %>% mutate(across(where(is.numeric), ~ round(.x, 2)))))

save(hindcast, summary_tab, current, file = "output/hindcast_results.rda")

write_xlsx(list(
  hindcast = hindcast %>% mutate(across(where(is.numeric), ~ round(.x, 2))),
  summary  = summary_tab %>% mutate(across(where(is.numeric), ~ round(.x, 2)))),
  path = "output/hindcast_results.xlsx")

# --- 6. Diagnostic plots ----------------------------------------------------
cols <- c(Logistic = "#1b6ca8", MichaelisMenten = "#d1495b",
          QuasiPoisson_nonzero = "#8f8f8f", QuasiPoisson_zerofill = "#2a9d8f")

p_err <- ggplot(hindcast, aes(cutoff, pred_new, colour = model)) +
  geom_line() + geom_point() +
  geom_line(aes(y = observed_new), colour = "black", linetype = "dashed") +
  geom_point(aes(y = observed_new), colour = "black", shape = 4) +
  facet_wrap(~ group, scales = "free_y") +
  scale_colour_manual(values = cols) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = c(0, 10, 30, 100, 300, 1000, 3000)) +
  labs(x = "Cutoff year (models fitted with data up to this year)",
       y = paste0("New species described after cutoff until ", year_end),
       subtitle = "Dashed black = observed") +
  theme_light() + theme(legend.position = "bottom")
ggsave("output/diagnostic_hindcast_new_species.png", p_err, width = 10, height = 7, dpi = 200)

# quasi-Poisson K can diverge to absurd values (> 1e10); not plotted
p_K <- ggplot(bind_rows(hindcast %>% select(group, model, cutoff, K), current) %>%
                filter(K < 1e5),
              aes(cutoff, K, colour = model)) +
  geom_line() + geom_point() +
  geom_line(data = series_all %>% group_by(group) %>%
              reframe(cutoff = c(cutoffs, year_end),
                        K = sapply(c(cutoffs, year_end), function(cu) max(cum_y[year <= cu]))),
            colour = "black", linetype = "dashed") +
  facet_wrap(~ group, scales = "free_y") +
  scale_colour_manual(values = cols) +
  scale_y_log10() +
  labs(x = "Data up to year", y = "Estimated total species (log scale)",
       subtitle = "Stability of the estimate as data accumulate; dashed black = described so far") +
  theme_light() + theme(legend.position = "bottom")
ggsave("output/diagnostic_hindcast_K_stability.png", p_K, width = 10, height = 7, dpi = 200)
