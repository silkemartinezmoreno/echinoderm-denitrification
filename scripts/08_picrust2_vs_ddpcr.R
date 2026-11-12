# =============================================================================
# 08. PICRUSt2-predicted vs ddPCR-quantified nirS and nosZ (Figures 6 and 7)
# -----------------------------------------------------------------------------
# Manuscript section: Results, "Correlations between the PICRUSt2-predicted and
# ddPCR quantified gene abundances".
#
# Spearman correlations between log10(ddPCR copies/ng-DNA + 1) and
# log10(PICRUSt2 abundance + 1), separately for biofilm and subcuticle, using
# all samples with both measurements. The linear regression (slope) is the
# trend line drawn in the figures.
#
# Input:  results/tables/ddpcr_curated.csv, data/picrust2/denitrification_KOs_contrib.xlsx
# Output: results/tables/Fig6-7_picrust2_vs_ddpcr_per_sample.csv
#         results/tables/Fig6-7_spearman_correlations.csv
# =============================================================================

source(file.path("scripts", "00_settings.R"))

ddpcr_column <- c(nirS = "log10nirs", nosZ = "log10nosz")

picrust2 <- read_picrust2()
ddpcr <- read_ddpcr() %>% filter(type %in% c("biofilm", "subcuticle"))

per_sample <- list()
correlations <- list()

for (g in names(ddpcr_column)) {
  d <- ddpcr %>%
    inner_join(filter(picrust2, gene == g), by = "code") %>%
    mutate(log10_ddpcr   = .data[[ddpcr_column[g]]],
           log10_picrust2 = log10(abundance_total + 1)) %>%
    filter(!is.na(log10_ddpcr))
  per_sample[[g]] <- select(d, gene, sample, code, type, class, sp, log10_ddpcr, log10_picrust2)

  for (t in c("biofilm", "subcuticle")) {
    dt <- d[d$type == t, ]
    # Ties in the data: p-values use the asymptotic approximation
    ct  <- suppressWarnings(cor.test(dt$log10_ddpcr, dt$log10_picrust2, method = "spearman"))
    fit <- lm(log10_picrust2 ~ log10_ddpcr, data = dt)
    correlations[[paste(g, t)]] <- data.frame(
      gene = g, type = t, n = nrow(dt),
      spearman_rho = unname(ct$estimate), p_value = ct$p.value,
      lm_intercept = coef(fit)[[1]], lm_slope = coef(fit)[[2]]
    )
  }
}

write.csv(bind_rows(per_sample), file.path(tables_dir, "Fig6-7_picrust2_vs_ddpcr_per_sample.csv"), row.names = FALSE)
write.csv(bind_rows(correlations), file.path(tables_dir, "Fig6-7_spearman_correlations.csv"), row.names = FALSE)
