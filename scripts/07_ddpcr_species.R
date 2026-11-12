# =============================================================================
# 07. nirS, nosZ and nirS/nosZ ratio per species (Figures 3, 4 and 5)
# -----------------------------------------------------------------------------
# Manuscript sections: Results, "ddPCR quantification of nirS and nosZ genes"
# and "The nirS/nosZ ratio".
#
# Figure 5: the nirS/nosZ ratio is computed on the log10-transformed
# abundances, log10(nirS + 1) / log10(nosZ + 1). Values > 1 indicate relative
# predominance of nirS and values < 1 of nosZ. The ratio of the untransformed
# abundances (nirS / nosZ) is also given for reference.
#
# Input:  results/tables/ddpcr_curated.csv
# Output: results/tables/Fig3-5_ddpcr_per_sample.csv
#         results/tables/Fig3_nirS_by_species.csv
#         results/tables/Fig4_nosZ_by_species.csv
#         results/tables/Fig5_ratio_by_species.csv
# =============================================================================

source(file.path("scripts", "00_settings.R"))

ddpcr <- read_ddpcr() %>%
  filter(type %in% c("biofilm", "subcuticle"), sp %in% species_11) %>%
  mutate(sp = factor(sp, levels = species_11)) %>%
  arrange(type, sp, sample)

write.csv(select(ddpcr, sample, code, type, class, sp, copiasnirs, copiasnosz,
                 log10nirs, log10nosz, rationirsnosz, log10ratio),
          file.path(tables_dir, "Fig3-5_ddpcr_per_sample.csv"), row.names = FALSE)

by_species <- function(variable) {
  ddpcr %>%
    group_by(type, class, sp) %>%
    group_modify(~ describe(.x[[variable]])) %>%
    ungroup()
}

fig3 <- left_join(by_species("copiasnirs"),
                  select(by_species("log10nirs"), type, sp, log10_median = median),
                  by = c("type", "sp"))
fig4 <- left_join(by_species("copiasnosz"),
                  select(by_species("log10nosz"), type, sp, log10_median = median),
                  by = c("type", "sp"))
fig5 <- bind_rows(
  mutate(by_species("log10ratio"), variable = "log10(nirS + 1) / log10(nosZ + 1)"),
  mutate(by_species("rationirsnosz"), variable = "nirS / nosZ")
) %>%
  select(variable, everything())

write.csv(fig3, file.path(tables_dir, "Fig3_nirS_by_species.csv"), row.names = FALSE)
write.csv(fig4, file.path(tables_dir, "Fig4_nosZ_by_species.csv"), row.names = FALSE)
write.csv(fig5, file.path(tables_dir, "Fig5_ratio_by_species.csv"), row.names = FALSE)
