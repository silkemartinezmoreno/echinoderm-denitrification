# =============================================================================
# 00. Shared settings for scripts 04-09
# -----------------------------------------------------------------------------
# Sourced at the start of each analysis script. Paths are relative to the
# project folder (open echinoderm-denitrification.Rproj in RStudio).
# =============================================================================

library(readxl)
library(dplyr)

tables_dir <- file.path("results", "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# Samples excluded from all ddPCR analyses. In each of them, one of the three
# pseudoreplicates deviated strongly from the other two (relative deviation
# > 0.25) but was not flagged by the robust z-score, because |z| < 3.5:
#   12_PL2_S  nirS  248 vs 3247 and 3900 copies/ng-DNA (z = -3.09)
#   1_HT2_H   nosZ   27 vs  296 and  357 copies/ng-DNA (z = -3.02)
#   8_MG1_S   nosZ   14 vs  4.3 and  6.5 copies/ng-DNA (z =  2.40)
# These samples were considered unreliable and removed from the analyses
# (see results/tables/ddpcr_replicates.csv).
excluded_samples <- c("12_PL2_S", "8_MG1_S", "1_HT2_H")

# The 11 species with nirS and nosZ ddPCR data in both compartments
# (biofilm and subcuticle), shown in Figures 3-5
species_11 <- c(
  "Antedon_mediterranea",                                   # Crinoidea
  "Ophioderma_longicauda", "Ophiothrix_fragilis",           # Ophiuroidea
  "Marthasterias_glacialis",                                # Asteroidea
  "Paracentrotus_lividus", "Arbacia_lixula",                # Echinoidea
  "Arbaciella_elegans", "Sphaerechinus_granularis",
  "Holothuria_tubulosa", "Holothuria_arguinensis",          # Holothuroidea
  "Holothuria_sanctori"
)

# Class names used in the data (Spanish) and in the manuscript
class_names <- c(
  "crinoideo" = "Crinoidea",
  "ofiura"    = "Ophiuroidea",
  "estrella"  = "Asteroidea",
  "erizo"     = "Echinoidea",
  "holoturia" = "Holothuroidea",
  "Agua"      = "Seawater",
  "Sedimento" = "Sediment"
)

# Curated ddPCR table (script 04)
read_ddpcr <- function() {
  d <- read.csv(file.path(tables_dir, "ddpcr_curated.csv"))
  d[!d$sample %in% excluded_samples, ]
}

# PICRUSt2 abundance of each denitrification gene per sample
read_picrust2 <- function() {
  read_excel(file.path("data", "picrust2", "denitrification_KOs_contrib.xlsx")) %>%
    group_by(code = sample, gene = gen_funcional) %>%
    summarise(abundance_total = sum(taxon_function_abun, na.rm = TRUE),
              abundance_rel   = sum(taxon_rel_function_abun, na.rm = TRUE),
              .groups = "drop")
}

# Descriptive statistics reported in the text
describe <- function(x) {
  x <- x[!is.na(x)]
  data.frame(n = length(x), median = median(x),
             q25 = quantile(x, 0.25, names = FALSE), q75 = quantile(x, 0.75, names = FALSE),
             min = min(x), max = max(x))
}
