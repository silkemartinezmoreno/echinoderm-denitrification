# =============================================================================
# 03. Extraction of denitrification KOs from the PICRUSt2 output
# -----------------------------------------------------------------------------
# Manuscript section: Methods, "Functional inference using PICRUSt2".
#
# pred_metagenome_contrib.tsv (script 02) is too large to load in memory, so
# it is read in chunks of 1e6 lines and only the KOs of interest are kept.
#
# Input:  results/picrust2/KO_metagenome_out/pred_metagenome_contrib.tsv
#         (not included in the repository because of its size)
# Output: data/picrust2/denitrification_KOs_contrib.xlsx
#         (included in the repository; used by scripts 05, 08 and 09)
#
# Skip this script if you do not re-run PICRUSt2: its output is already provided.
# =============================================================================

library(data.table)
library(openxlsx)

contrib_file <- file.path("results", "picrust2", "KO_metagenome_out", "pred_metagenome_contrib.tsv")
output_file  <- file.path("data", "picrust2", "denitrification_KOs_contrib.xlsx")

# KEGG Orthologs of each step of denitrification (Equation 1), plus anammox KOs
# https://www.genome.jp/kegg/ko.html
gen_map <- c(
  "ko:K00370" = "narG",  # nitrate reductase (membrane-bound)
  "ko:K02567" = "napA",  # nitrate reductase (periplasmic)
  "ko:K15864" = "nirS",  # nitrite reductase (cytochrome cd1)
  "ko:K00368" = "nirK",  # nitrite reductase (copper)
  "ko:K04561" = "norB",  # nitric oxide reductase
  "ko:K00376" = "nosZ",  # nitrous oxide reductase
  "ko:K10535" = "hao",   # hydroxylamine oxidoreductase (anammox)
  "ko:K20935" = "hdh",   # hydrazine dehydrogenase (anammox)
  "ko:K20934" = "hzs"    # hydrazine synthase (anammox)
)
kos <- names(gen_map)

con <- file(contrib_file, "r")
total_size <- file.info(contrib_file)$size

header <- readLines(con, n = 1)
colnames <- strsplit(header, "\t")[[1]]
bytes_read <- nchar(header) + 1

result <- data.table()
repeat {
  lines <- readLines(con, n = 1e6)
  if (length(lines) == 0) break

  bytes_read <- bytes_read + sum(nchar(lines)) + length(lines)
  cat(sprintf("Progress: %.2f%%\n", bytes_read / total_size * 100))

  chunk <- fread(text = paste(lines, collapse = "\n"), sep = "\t", header = FALSE)
  setnames(chunk, colnames)
  result <- rbind(result, chunk[chunk[["function"]] %in% kos])
}
close(con)

result[, gen_funcional := gen_map[result[["function"]]]]
write.xlsx(result, output_file)
