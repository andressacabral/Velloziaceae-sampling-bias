# 101 - Prepare the Velloziaceae species list for the analyses
#
# Input : data/vell_acc_spp_year              (all accepted species, WCVP)
#         data/vell_acc_spp_year_<region>      (same species split by region)
# Output: output/velloziaceae_data_for_analyses.rda  (object `outp`)
#         output/velloziaceae_data_for_analyses.xlsx
#
# Column names follow the Bromeliaceae pipeline (how_many_bromeliads) so the
# downstream scripts (20*, 30*) can be reused with minimal changes.

library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(writexl)

dir.create("output", showWarnings = FALSE)

# --- 1. Read data ------------------------------------------------------------
inp <- read_csv("data/vell_acc_spp_year", show_col_types = FALSE, name_repair = "unique_quiet") %>%
  select(-any_of("...1")) %>%
  mutate(plant_name_id = as.character(plant_name_id),
         taxon_name    = str_squish(taxon_name),
         authors       = str_squish(authors))

# Regions
# NOTE: in the raw data the file names for Neotropics and Afrotropics are swapped
# (checked 2026-09-23: "_afrotropics" contains Vellozia/Barbacenia, "_neotropics"
# contains Xerophyta/Talbotia). They are relabelled here. If the raw files are
# renamed, fix this vector - the check below will stop the script otherwise.
region_files <- c(Neotropics  = "data/vell_acc_spp_year_afrotropics",
                  Afrotropics = "data/vell_acc_spp_year_neotropics",
                  Asia        = "data/vell_acc_spp_year_asia")

regions <- bind_rows(lapply(names(region_files), function(r) {
  read_csv(region_files[[r]], show_col_types = FALSE, name_repair = "unique_quiet",
           col_select = "plant_name_id") %>%
    mutate(plant_name_id = as.character(plant_name_id), region = r)
}))

# each species must be in exactly one region
stopifnot(!any(duplicated(regions$plant_name_id)))

dat <- inp %>%
  left_join(regions, by = "plant_name_id") %>%
  mutate(tax_genus = word(taxon_name, 1))

# sanity check of the region assignment (Barbacenia/Vellozia are Neotropical)
chk <- dat %>% filter(tax_genus %in% c("Barbacenia", "Vellozia"))
if (any(chk$region != "Neotropics")) {
  stop("Region labels look wrong: Barbacenia/Vellozia not in Neotropics. Check `region_files`.")
}

# --- 2. Filters (same logic as Bromeliaceae) ----------------------------------
dat <- dat %>%
  mutate(is_hybrid     = str_detect(taxon_name, "\u00d7| x "),
         is_infraspecific = str_detect(taxon_name, " (subsp\\.|var\\.|f\\.) "))

message("Hybrids removed: ", sum(dat$is_hybrid),
        " | infraspecific taxa removed: ", sum(dat$is_infraspecific))

dat <- dat %>% filter(!is_hybrid, !is_infraspecific)

# --- 3. Year of first description -----------------------------------------
# For recombinations the year of the accepted name is the year of the new
# combination, not of the first description. We need the basionym year.
# Recombinations are recognised by a parenthetical author: "(Hochst.) Baker".
dat <- dat %>%
  mutate(is_recombination = str_detect(authors, "^\\("),
         parenthetical_author = ifelse(is_recombination,
                                       str_match(authors, "^\\(([^)]*)\\)")[, 2],
                                       NA_character_),
         author_recombination = ifelse(is_recombination,
                                       str_squish(str_remove(authors, "^\\([^)]*\\)")),
                                       NA_character_))

# Basionym information from WCVP (package rWCVPdata), if installed:
#   install.packages("rWCVPdata", repos = c("https://matildabrown.github.io/drat", "https://cloud.r-project.org"))
if (requireNamespace("rWCVPdata", quietly = TRUE)) {
  wcvp <- rWCVPdata::wcvp_names %>%
    mutate(plant_name_id = as.character(plant_name_id),
           basionym_plant_name_id = as.character(basionym_plant_name_id))

  bas <- wcvp %>%
    filter(plant_name_id %in% dat$plant_name_id) %>%
    select(plant_name_id, basionym_plant_name_id) %>%
    filter(!is.na(basionym_plant_name_id), basionym_plant_name_id != "") %>%
    left_join(wcvp %>%
                transmute(basionym_plant_name_id = plant_name_id,
                          basionym = str_squish(paste(taxon_name, taxon_authors)),
                          basionym_author = taxon_authors,
                          basionym_year = as.numeric(str_extract(first_published, "\\d{4}"))),
              by = "basionym_plant_name_id")

  dat <- dat %>% left_join(bas, by = "plant_name_id")
  message("Basionyms found in WCVP: ", sum(!is.na(dat$basionym)))
} else {
  message("rWCVPdata not installed: basionym years not added. ",
          "Recombinations keep the year of the combination (see is_recombination).")
  dat <- dat %>% mutate(basionym_plant_name_id = NA_character_,
                        basionym = NA_character_,
                        basionym_author = NA_character_,
                        basionym_year = NA_real_)
}

# recombinations without a basionym year -> flag, keep combination year
no_bas <- dat %>% filter(is_recombination, is.na(basionym_year))
if (nrow(no_bas) > 0) {
  warning(nrow(no_bas), " recombination(s) without basionym year; using the year of ",
          "the combination: ", paste(no_bas$taxon_name, collapse = ", "))
}

dat <- dat %>%
  mutate(year_first_description = ifelse(!is.na(basionym_year), basionym_year, year),
         author_first_description = case_when(
           !is.na(basionym_author)      ~ basionym_author,
           !is.na(parenthetical_author) ~ parenthetical_author,
           TRUE                         ~ authors),
         species_name_first_description = ifelse(!is.na(basionym), basionym,
                                                 paste(taxon_name, authors)))

# --- 4. Time bins -----------------------------------------------------------
dat <- dat %>%
  mutate(decade_first_description = cut(year_first_description,
                                        breaks = seq(1750, 2030, by = 10),
                                        labels = seq(1755, 2025, by = 10)),
         halfcentury_first_description = cut(year_first_description,
                                             breaks = c(1750, 1800, 1850, 1900, 1950, 2000, 2050),
                                             labels = c("1750-1800", "1801-1850", "1851-1900",
                                                        "1901-1950", "1951-2000", "2001-2025")))

# --- 5. Final table ---------------------------------------------------------
outp <- dat %>%
  transmute(plant_name_id,
            tax_accepted_name = taxon_name,
            tax_authors = authors,
            tax_genus,
            region,
            year,                               # year of the accepted name
            publication,
            volume_page,
            is_recombination,
            author_recombination,
            basionym,
            basionym_year,
            year_first_description,             # basionym year if recombination
            author_first_description,
            species_name_first_description,
            decade_first_description,
            halfcentury_first_description) %>%
  arrange(tax_accepted_name)

# --- 6. Checks --------------------------------------------------------------
stopifnot(!any(duplicated(outp$plant_name_id)),
          !any(is.na(outp$year_first_description)),
          !any(is.na(outp$region)))

message("Species retained: ", nrow(outp),
        " | first description: ", min(outp$year_first_description),
        "-", max(outp$year_first_description))
print(count(outp, region, tax_genus))

# --- 7. Write to disk -------------------------------------------------------
write_xlsx(outp, "output/velloziaceae_data_for_analyses.xlsx")
save(outp, file = "output/velloziaceae_data_for_analyses.rda")
