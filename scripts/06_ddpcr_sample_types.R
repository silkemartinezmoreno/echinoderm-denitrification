# =============================================================================
# 06. nirS and nosZ abundance across sample types (Figure S1)
# -----------------------------------------------------------------------------
# Manuscript section: Results, "ddPCR quantification of nirS and nosZ genes".
#
# Biofilm, subcuticle, seawater and sediment samples are compared with
# Kruskal-Wallis tests followed by Dunn's post-hoc tests with Bonferroni
# correction (two-sided p-values). Both tests are rank-based, so the results
# are the same for copies/ng-DNA and for its log10 transformation.
#
# Input:  results/tables/ddpcr_curated.csv
# Output: results/tables/FigS1_ddpcr_by_sample_type.csv   (descriptive statistics)
#         results/tables/FigS1_kruskal_wallis.csv
#         results/tables/FigS1_dunn_bonferroni.csv
# =============================================================================

source(file.path("scripts", "00_settings.R"))
library(dunn.test)

type_names <- c("biofilm" = "Biofilm", "subcuticle" = "Subcuticle",
                "Agua" = "Seawater", "Sedimento" = "Sediment")

ddpcr <- read_ddpcr() %>%
  filter(type %in% names(type_names)) %>%
  mutate(sample_type = factor(type_names[type], levels = type_names))

genes <- c(nirS = "copiasnirs", nosZ = "copiasnosz")

descriptive <- list()
kruskal <- list()
dunn <- list()

for (gene in names(genes)) {
  d <- ddpcr[!is.na(ddpcr[[genes[gene]]]), ]
  x <- d[[genes[gene]]]

  for (st in levels(d$sample_type)) {
    descriptive[[paste(gene, st)]] <- data.frame(gene = gene, sample_type = st,
                                                 describe(x[d$sample_type == st]))
  }

  kw <- kruskal.test(x ~ d$sample_type)
  kruskal[[gene]] <- data.frame(gene = gene, n = length(x), chi_squared = unname(kw$statistic),
                                df = unname(kw$parameter), p_value = kw$p.value)

  dt <- dunn.test(x, d$sample_type, method = "bonferroni", altp = TRUE, kw = FALSE, table = FALSE)
  dunn[[gene]] <- data.frame(gene = gene, comparison = dt$comparisons, z = dt$Z,
                             p_value = dt$altP, p_bonferroni = dt$altP.adjusted)
}

write.csv(bind_rows(descriptive), file.path(tables_dir, "FigS1_ddpcr_by_sample_type.csv"), row.names = FALSE)
write.csv(bind_rows(kruskal), file.path(tables_dir, "FigS1_kruskal_wallis.csv"), row.names = FALSE)
write.csv(bind_rows(dunn), file.path(tables_dir, "FigS1_dunn_bonferroni.csv"), row.names = FALSE)
