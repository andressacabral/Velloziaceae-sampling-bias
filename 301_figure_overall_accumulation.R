# 301 - Main figure: accumulation of species descriptions with model curves
#
# (a) Velloziaceae overall, (b) per region and genus.
# Points = observed cumulative number of species (by year of first description).
# Lines + bands = logistic and quasi-Poisson projections with 95% CI.
# Groups with too little information for an asymptote (Afrotropics, Xerophyta;
# see 202_hindcast.R) are shown without model curves.
#
# Input : output/peryear.rda, output/band_*_df.rda, output/band_qp_df.rda,
#         output/overall_table.rda   (all from 201)
# Output: figures/figure_1_accumulation.pdf / .png

library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)
library(readxl)

load("output/peryear.rda")
load("output/band_overall_df.rda")
load("output/band_region_df.rda")
load("output/band_genus_df.rda")
load("output/band_qp_df.rda")

# --- Settings ----------------------------------------------------------------
# Which model is the main estimate? "both" = equal weight (until decided),
# "Logistic" or "QuasiPoisson" = that one solid, the other drawn lighter.
main_model <- "both"

# Groups shown without an estimate (logistic asymptote = observed number,
# quasi-Poisson does not converge; see 202_hindcast.R)
no_estimate <- c("Afrotropics", "Xerophyta")

# Panels and order for (b). Xerophyta is not shown separately: the
# Afrotropics are Xerophyta (56 spp.) + Talbotia (1 sp.), so the two panels
# would be practically identical.
group_panels <- c("Neotropics", "Afrotropics", "Vellozia", "Barbacenia")

model_labels <- c(Logistic = "Logistic", QuasiPoisson = "Quasi-Poisson")
model_cols   <- c(Logistic = "#2a78d6", QuasiPoisson = "#eb6834")
model_lty    <- c(Logistic = "solid",   QuasiPoisson = "longdash")

alpha_line <- function(m) if (main_model %in% c("both", m)) 1 else 0.45
alpha_band <- function(m) if (main_model %in% c("both", m)) 0.18 else 0.07

# --- Data in one long format ------------------------------------------------
pts <- bind_rows(
  peryear %>% mutate(group = "Velloziaceae"),
  peryear_region %>% rename(group = region),
  peryear_genus  %>% rename(group = tax_genus)) %>%
  select(group, year, cum_y)

# Afrotropics comes from the region series, Xerophyta from the genus series;
# they were dropped from the model bands only below, not in 201.
bands <- bind_rows(
  band_overall_df %>% filter(model == "Logistic") %>% mutate(group = "Velloziaceae"),
  band_qp_overall %>% mutate(group = "Velloziaceae"),
  band_region %>% rename(group = region),
  band_qp_region %>% rename(group = region),
  band_genus %>% rename(group = tax_genus),
  band_qp_genus %>% rename(group = tax_genus)) %>%
  filter(!group %in% no_estimate) %>%
  mutate(model = recode(model, Logistic = "Logistic"),
         alpha_l = vapply(model, alpha_line, numeric(1)),
         alpha_b = vapply(model, alpha_band, numeric(1)))

# Direct labels: estimated total (asymptote / K) with 95% CI, from 201,
# placed at the end of each curve (logistic above, quasi-Poisson below)
est_sheets <- c(overall = "overall", region = "region", genus = "genus")
estimates <- bind_rows(lapply(est_sheets, function(sh) {
  x <- readxl::read_xlsx("output/species_estimates.xlsx", sheet = sh)
  names(x)[1] <- if (sh == "overall") "model" else "group"
  if (sh == "overall") x$group <- "Velloziaceae"
  x %>% filter(model == "Logistic")
})) %>%
  transmute(group,
            Logistic     = sprintf("%.0f (%.0f\u2013%.0f)", estimated_species, lower_CI, upper_CI),
            QuasiPoisson = sprintf("%.0f (%.0f\u2013%.0f)", K_quasipoisson, K_quasipoisson_lwr,
                                   K_quasipoisson_upr)) %>%
  tidyr::pivot_longer(-group, names_to = "model", values_to = "estimate")

end_labels <- bands %>%
  group_by(group, model) %>%
  filter(year == max(year)) %>%
  ungroup() %>%
  left_join(estimates, by = c("group", "model")) %>%
  mutate(label = estimate,
         vjust = ifelse(model == "Logistic", -0.7, 1.7))

# --- Plot function ----------------------------------------------------------
plot_accum <- function(g, base_size = 9, label_size = 2.6, point_size = 0.9) {
  p <- pts %>% filter(group == g)
  b <- bands %>% filter(group == g)
  l <- end_labels %>% filter(group == g)
  last_year <- max(p$year)

  gg <- ggplot() +
    geom_vline(xintercept = last_year, colour = "grey75", linewidth = 0.3, linetype = "dotted")

  if (nrow(b) > 0) {
    gg <- gg +
      geom_ribbon(data = b, aes(year, ymin = lwr, ymax = upr, fill = model, alpha = I(alpha_b)),
                  colour = NA) +
      geom_line(data = b, aes(year, y_hat, colour = model, linetype = model, alpha = I(alpha_l)),
                linewidth = 0.6) +
      geom_text(data = l, aes(year, y_hat, label = label, colour = model, vjust = vjust),
                hjust = 1, size = label_size, show.legend = FALSE)
  } else {
    gg <- gg +
      annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.6, size = label_size,
               colour = "grey35", label = "No reliable estimate")
  }

  gg +
    geom_point(data = p, aes(year, cum_y), size = point_size, colour = "grey15") +
    scale_colour_manual(values = model_cols, labels = model_labels, name = NULL) +
    scale_fill_manual(values = model_cols, labels = model_labels, name = NULL, guide = "none") +
    scale_linetype_manual(values = model_lty, labels = model_labels, name = NULL) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
    labs(title = g, x = NULL, y = NULL) +
    theme_classic(base_size = base_size) +
    theme(plot.title = element_text(face = if (g %in% c("Vellozia", "Barbacenia", "Xerophyta"))
                                      "italic" else "plain", size = base_size),
          axis.line = element_line(colour = "grey40", linewidth = 0.3),
          axis.ticks = element_line(colour = "grey40", linewidth = 0.3),
          panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
          legend.position = "bottom")
}

# --- Assemble ---------------------------------------------------------------
p_a <- plot_accum("Velloziaceae", base_size = 10, label_size = 3, point_size = 1.1) +
  labs(y = "Cumulative number of species")

p_b <- wrap_plots(lapply(group_panels, function(g) plot_accum(g) + theme(legend.position = "none")),
                  ncol = 2)

fig <- (p_a / p_b) +
  plot_layout(heights = c(1.1, 1), guides = "collect") +
  plot_annotation(tag_levels = list(c("a", "b", "", "", "")),
                  theme = theme(legend.position = "bottom")) &
  theme(legend.text = element_text(size = 9),
        legend.key.width = unit(10, "mm")) &
  guides(colour = guide_legend(override.aes = list(linewidth = 0.8, alpha = 1)))

dir.create("figures", showWarnings = FALSE)
ggsave("figures/figure_1_accumulation.pdf", fig, width = 180, height = 210, units = "mm",
       device = cairo_pdf)
ggsave("figures/figure_1_accumulation.png", fig, width = 180, height = 210, units = "mm", dpi = 300)
