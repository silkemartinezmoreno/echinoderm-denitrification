# =============================================================================
# 09. Number of samples per species and analysis (Table S1)
# -----------------------------------------------------------------------------
# Counts the samples of each species and compartment that enter each analysis.
# Run scripts 05-08 first.
#
# Output: results/tables/TableS1_samples_per_analysis.csv
# =============================================================================

source(file.path("scripts", "00_settings.R"))

count_samples <- function(d, label) {
  d %>% count(type, sp, name = label)
}

fig2 <- read.csv(file.path(tables_dir, "Fig2_picrust2_genes_by_species.csv"))
fig35 <- read.csv(file.path(tables_dir, "Fig3-5_ddpcr_per_sample.csv"))
fig67 <- read.csv(file.path(tables_dir, "Fig6-7_picrust2_vs_ddpcr_per_sample.csv"))

counts <- list(
  fig2 %>% group_by(type, sp) %>% summarise(Fig2_picrust2 = max(n_samples), .groups = "drop"),
  count_samples(filter(fig35, !is.na(copiasnirs)), "Fig3_nirS_ddpcr"),
  count_samples(filter(fig35, !is.na(copiasnosz)), "Fig4_nosZ_ddpcr"),
  count_samples(filter(fig35, !is.na(rationirsnosz)), "Fig5_ratio"),
  count_samples(filter(fig67, gene == "nirS"), "Fig6_nirS_correlation"),
  count_samples(filter(fig67, gene == "nosZ"), "Fig7_nosZ_correlation")
)

species_class <- distinct(read_ddpcr(), sp, class)

table_s1 <- Reduce(function(a, b) full_join(a, b, by = c("type", "sp")), counts) %>%
  left_join(species_class, by = "sp") %>%
  mutate(across(starts_with("Fig"), ~ coalesce(.x, 0L)),
         class = factor(class, levels = class_names[1:5])) %>%
  select(type, class, sp, everything()) %>%
  arrange(type, class, sp)

totals <- table_s1 %>%
  group_by(type) %>%
  summarise(across(starts_with("Fig"), sum)) %>%
  mutate(sp = "Total")

write.csv(bind_rows(table_s1, totals), file.path(tables_dir, "TableS1_samples_per_analysis.csv"),
          row.names = FALSE, na = "")
