# =============================================================================
# 05. PICRUSt2-predicted abundance of denitrification genes (Figure 2)
# -----------------------------------------------------------------------------
# Manuscript section: Results, "Estimation of denitrification potential with
# PICRUSt2".
#
# For each sample, the abundance of each gene is the sum of the contributions
# of all OTUs. The values shown in Figure 2 are the mean per host species and
# compartment, using the samples of the curated ddPCR dataset.
#
# Input:  data/picrust2/denitrification_KOs_contrib.xlsx, results/tables/ddpcr_curated.csv
# Output: results/tables/Fig2_picrust2_genes_by_species.csv
# =============================================================================

source(file.path("scripts", "00_settings.R"))

genes <- c("narG", "napA", "nirS", "nirK", "norB", "nosZ")

picrust2 <- read_picrust2() %>% filter(gene %in% genes)

# Only samples with PICRUSt2 data are kept (code "X" = not sequenced)
fig2 <- read_ddpcr() %>%
  filter(type %in% c("biofilm", "subcuticle")) %>%
  inner_join(picrust2, by = "code") %>%
  group_by(type, class, sp, gene) %>%
  summarise(n_samples       = n(),
            abundance_total = mean(abundance_total, na.rm = TRUE),
            abundance_rel   = mean(abundance_rel, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(class = factor(class, levels = class_names[1:5]),
         gene  = factor(gene, levels = genes)) %>%
  arrange(type, class, sp, gene)

write.csv(fig2, file.path(tables_dir, "Fig2_picrust2_genes_by_species.csv"), row.names = FALSE)
