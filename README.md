# Velloziaceae sampling bias

How many Velloziaceae species are there? Species description history and
asymptotic estimates of total species richness, following the workflow of
Zizka et al. (2026) *How many Bromeliads are there?*
([code](https://doi.org/10.5281/zenodo.21789730)).

## Structure

```
data/                      raw input (WCVP accepted species + year), do not edit
000_helper_functions.R     model fitting, bootstrap CIs, quasi-Poisson check
10*_                       data preparation
20*_                       analyses
30*_                       figures
40*_                       supplementary material
output/                    intermediate results (not tracked)
figures/                   figures (not tracked)
supplementary_material/    supplementary files (not tracked)
```

Run the scripts in numerical order from the project root
(open `Velloziaceae-sampling-bias.Rproj`).

## Workflow

| Script | Status | Description |
|---|---|---|
| `101_prepare_species_list.R` | done | Cleans species list, assigns region, year of first description (basionym year for recombinations) |
| `201_asymptotic_models.R` | to do | Logistic / Michaelis-Menten fits, total and per region / genus |
| `301_figure_overall_accumulation.R` | to do | Main accumulation figure |

## Notes on the data

- `data/vell_acc_spp_year_afrotropics` and `..._neotropics` have swapped names
  (Afrotropics file contains the Neotropical genera). Corrected in `101`.
- `year` in the raw data is the year of the accepted name. For recombinations
  the basionym year is taken from WCVP via the `rWCVPdata` package if installed:
  `install.packages("rWCVPdata", repos = c("https://matildabrown.github.io/drat", "https://cloud.r-project.org"))`

## Packages

dplyr, tidyr, stringr, readr, purrr, writexl, broom, minpack.lm, MASS, ggplot2
(optional: rWCVPdata)
