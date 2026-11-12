# =============================================================================
# Reproduces all the results of the manuscript from the data in data/
# -----------------------------------------------------------------------------
# Open echinoderm-denitrification.Rproj in RStudio and run this file.
# Result tables are written to results/tables/.
#
# Scripts 01-03 (sequence processing, PICRUSt2 and KO extraction) are not run
# here: they need the raw reads and many hours of computation, and their
# output is already provided in data/picrust2/. See README.md.
# =============================================================================

scripts <- c(
  "04_ddpcr_curation.R",                # ddPCR pseudoreplicates -> copies/ng-DNA per sample
  "05_picrust2_denitrification_genes.R", # Figure 2
  "06_ddpcr_sample_types.R",            # Figure S1
  "07_ddpcr_species.R",                 # Figures 3, 4 and 5
  "08_picrust2_vs_ddpcr.R",             # Figures 6 and 7
  "09_table_S1.R"                       # Table S1
)

for (s in scripts) {
  message("Running ", s)
  source(file.path("scripts", s), local = new.env())
}

writeLines(capture.output(sessionInfo()), file.path("results", "sessionInfo.txt"))
