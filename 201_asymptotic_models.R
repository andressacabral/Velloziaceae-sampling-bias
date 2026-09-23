# 201 - Asymptotic models of species description (how many Velloziaceae are there?)
#
# Fits logistic and Michaelis-Menten curves to the cumulative number of first
# descriptions (overall, per region, per genus), with parametric-bootstrap CIs
# for the asymptote (a = estimated total number of species), and a
# quasi-Poisson GLM on annual counts as an independent check (Bebber et al. 2007).
# Method as in Zizka et al. (2026) How many Bromeliads are there?
#
# Input : output/velloziaceae_data_for_analyses.rda   (from 101)
# Output: output/peryear*.rda, output/overall_models.rda, output/overall_table.rda,
#         output/band_*.rda, output/species_estimates_*.xlsx,
#         output/diagnostic_*.png

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

# Groups with fewer species than this are not modelled (too few points for a
# sigmoid; e.g. Asia = 2 spp., Barbaceniopsis = 4 spp.)
min_species <- 20

# Years to extrapolate the curves into the future
extra_years_overall <- 150
extra_years_groups  <- 80

# --- 1. Cumulative series --------------------------------------------------
peryear        <- build_cumulative(outp)
peryear_region <- build_cumulative(outp, "region")
peryear_genus  <- build_cumulative(outp, "tax_genus")

keep_groups <- function(series, group_col) {
  n <- series %>%
    group_by(across(all_of(group_col))) %>%
    summarise(n_species = max(cum_y), n_years = n(), .groups = "drop")
  skipped <- n %>% filter(n_species < min_species)
  if (nrow(skipped) > 0) {
    message("Not modelled (< ", min_species, " spp.): ",
            paste0(skipped[[group_col]], " (", skipped$n_species, ")", collapse = ", "))
  }
  series %>% filter(.data[[group_col]] %in% n[[group_col]][n$n_species >= min_species])
}

peryear_region <- keep_groups(peryear_region, "region")
peryear_genus  <- keep_groups(peryear_genus, "tax_genus")

save(peryear, peryear_region, peryear_genus, file = "output/peryear.rda")

# --- 2. Overall -------------------------------------------------------------
overall_models <- fit_models(peryear)

overall_table <- overall_models %>%
  mutate(a_ci = map(fit, ~ get_a_ci_boot(.x, level = 0.95,
                                         k_bounds = c(1e-4, 0.5),
                                         t_bounds = range(peryear$t))),
         AIC  = map_dbl(glance, ~ .x$AIC)) %>%
  unnest(a_ci) %>%
  mutate(observed     = max(peryear$cum_y),
         missing      = pmax(0, a_estimate - observed),
         complete_pct = 100 * observed / a_estimate,
         year_min     = min(peryear$year),
         year_latest  = max(peryear$year)) %>%
  select(model, a_estimate, a_lwr, a_upr, observed, missing, complete_pct,
         AIC, year_min, year_latest)

# Quasi-Poisson check (after overall_table exists)
qp_overall <- fit_quasipoisson(peryear, t_col = "t", y_col = "y")
message(sprintf(
  "Quasi-Poisson K (overall): %.0f | completeness: %.3f | ratio vs logistic: %.3f",
  qp_overall$K_qp, qp_overall$completeness_qp,
  qp_overall$K_qp / overall_table$a_estimate[overall_table$model == "Logistic"]))

overall_table <- overall_table %>%
  mutate(K_quasipoisson = qp_overall$K_qp)

print(overall_table)

save(overall_models, file = "output/overall_models.rda")
save(overall_table,  file = "output/overall_table.rda")

# --- 3. Per group (region, genus) -------------------------------------------
fit_groups <- function(series, group_col) {
  series %>%
    group_by(across(all_of(group_col))) %>%
    nest() %>%
    ungroup() %>%
    mutate(models = map(data, ~ fit_models(.x, t_col = "t", cum_col = "cum_y"))) %>%
    unnest(models)
}

group_table <- function(fits, series, group_col) {
  obs <- series %>%
    group_by(across(all_of(group_col))) %>%
    summarise(observed    = max(cum_y),
              year_min    = min(year),
              year_latest = max(year),
              .groups = "drop")

  qp <- series %>%
    group_by(across(all_of(group_col))) %>%
    group_modify(~ fit_quasipoisson(.x, t_col = "t", y_col = "y")) %>%
    ungroup() %>%
    select(all_of(group_col), K_quasipoisson = K_qp)

  fits %>%
    mutate(a_ci = map2(fit, data, ~ get_a_ci_boot(.x, level = 0.95,
                                                  k_bounds = c(1e-4, 0.5),
                                                  t_bounds = range(.y$t))),
           AIC  = map_dbl(glance, ~ .x$AIC)) %>%
    unnest(a_ci) %>%
    left_join(obs, by = group_col) %>%
    left_join(qp,  by = group_col) %>%
    mutate(missing      = pmax(0, a_estimate - observed),
           complete_pct = 100 * observed / a_estimate) %>%
    select(all_of(group_col), model, a_estimate, a_lwr, a_upr, observed, missing,
           complete_pct, AIC, K_quasipoisson, year_min, year_latest) %>%
    arrange(.data[[group_col]], model)
}

fits_region <- fit_groups(peryear_region, "region")
fits_genus  <- fit_groups(peryear_genus,  "tax_genus")

region_table <- group_table(fits_region, peryear_region, "region")
genus_table  <- group_table(fits_genus,  peryear_genus,  "tax_genus")

print(region_table)
print(genus_table)

# --- 4. Export tables (rounded, logistic + MM) ------------------------------
tidy_out <- function(tab) {
  tab %>%
    mutate(across(c(a_estimate, a_lwr, a_upr), ~ round(.x, 1)),
           missing = round(missing, 0),
           complete_pct = round(complete_pct, 1),
           AIC = round(AIC, 1)) %>%
    rename(estimated_species   = a_estimate,
           lower_CI            = a_lwr,
           upper_CI            = a_upr,
           described_species   = observed,
           missing_species     = missing,
           percent_described   = complete_pct)
}

write_xlsx(list(overall = tidy_out(overall_table),
                region  = tidy_out(region_table),
                genus   = tidy_out(genus_table)),
           path = "output/species_estimates.xlsx")

# --- 5. CI bands over calendar years ---------------------------------------
band_overall_df <- overall_models %>%
  mutate(band = map2(fit, model, ~ predict_curve_band(
    .x, t_min = min(peryear$t), t_max = max(peryear$t),
    year_min = min(peryear$year), extra_years = extra_years_overall,
    model = .y, level = 0.95, B = 800))) %>%
  pull(band) %>%
  bind_rows()

save(band_overall_df, file = "output/band_overall_df.rda")

band_groups <- function(fits, group_col, strict_groups = character(0)) {
  fits %>%
    filter(model == "Logistic") %>%
    mutate(band = pmap(list(.data[[group_col]], data, fit), function(g, df, fitobj) {
      pars <- get_band_params(g, strict_groups)
      res <- predict_curve_band(fitobj,
                                t_min = min(df$t), t_max = max(df$t),
                                year_min = min(df$year),
                                extra_years = extra_years_groups,
                                model = "Logistic", level = 0.95, B = 600,
                                k_bounds = pars$k_bounds, max_mult = pars$max_mult)
      if (is.null(res)) {
        message("No CI band for ", g, " (too few valid bootstrap draws)")
        return(NULL)
      }
      mutate(res, !!group_col := g)
    })) %>%
    pull(band) %>%
    compact() %>%
    bind_rows()
}

band_region <- band_groups(fits_region, "region")
band_genus  <- band_groups(fits_genus,  "tax_genus")

save(band_region, file = "output/band_region_df.rda")
save(band_genus,  file = "output/band_genus_df.rda")

# --- 6. Diagnostic plots (not for the manuscript) ---------------------------
p_overall <- ggplot() +
  geom_ribbon(data = band_overall_df,
              aes(year, ymin = lwr, ymax = upr, fill = model), alpha = 0.25) +
  geom_line(data = band_overall_df, aes(year, y_hat, color = model), linewidth = 1) +
  geom_point(data = peryear, aes(year, cum_y), size = 0.6) +
  labs(x = "Year", y = "Cumulative species",
       title = "Velloziaceae: accumulation of first descriptions",
       subtitle = paste0("95% CI bands, +", extra_years_overall, " years")) +
  theme_light()

ggsave("output/diagnostic_overall.png", p_overall, width = 7, height = 5, dpi = 200)

plot_groups <- function(series, bands, group_col) {
  ggplot() +
    geom_ribbon(data = bands, aes(year, ymin = lwr, ymax = upr), alpha = 0.25) +
    geom_line(data = bands, aes(year, y_hat), linewidth = 0.8) +
    geom_point(data = series, aes(year, cum_y), size = 0.6) +
    facet_wrap(as.formula(paste("~", group_col)), scales = "free_y") +
    labs(x = "Year", y = "Cumulative species",
         subtitle = paste0("Logistic model, 95% CI bands, +", extra_years_groups, " years")) +
    theme_light()
}

ggsave("output/diagnostic_region.png",
       plot_groups(peryear_region, band_region, "region"),
       width = 9, height = 4, dpi = 200)
ggsave("output/diagnostic_genus.png",
       plot_groups(peryear_genus, band_genus, "tax_genus"),
       width = 11, height = 4, dpi = 200)
