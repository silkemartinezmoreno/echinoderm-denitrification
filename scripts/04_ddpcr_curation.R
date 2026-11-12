# =============================================================================
# 04. Curation of ddPCR data (nirS and nosZ)
# -----------------------------------------------------------------------------
# Manuscript section: Methods, "Data processing and visualisation".
#
# Each sample was run in triplicate (three pseudoreplicates). Outlier
# pseudoreplicates were identified with a robust z-score (median and MAD,
# |z| > 3.5) and then checked for a relative deviation from the mean > 0.25.
# A pseudoreplicate was removed only when both criteria were met. The
# pseudoreplicates finally removed are listed explicitly below.
#
# Input:  data/ddpcr/<gene>_<plate>_manual.xlsx  (QX Manager export, one file per plate)
#         data/ddpcr_sample_metadata.csv          (species, class, type and
#                                                  sequencing code of each sample)
# Output: results/tables/ddpcr_replicates.csv     (all pseudoreplicates and outlier criteria)
#         results/tables/ddpcr_curated.csv        (mean copies/ng-DNA per sample)
# =============================================================================

source(file.path("scripts", "00_settings.R"))

# Read the six plates of a gene. Only the sample name, the ng of DNA and the
# concentration in copies/ng-DNA (last column) are used. `row` is the row of
# each pseudoreplicate in its Excel file and identifies it below.
read_plates <- function(gene) {
  bind_rows(lapply(1:6, function(plate) {
    x <- read_excel(file.path("data", "ddpcr", sprintf("%s_%d_manual.xlsx", gene, plate)))
    x <- x[, c(1, 2, ncol(x))]
    colnames(x) <- c("sample", "ng_sample", "copias_ng")
    x$gene  <- gene
    x$plate <- plate
    x$row   <- seq_len(nrow(x))
    x
  }))
}

# Sample names: "1HT2S" -> "1_HT2_S". A final "D" (diluted re-run) is removed,
# except in sediment samples ("SED").
clean_names <- function(d) {
  d$sample_sep <- sub("^(\\d+)([A-Za-z]+)(\\d+)([A-Za-z]+)$", "\\1_\\2\\3_\\4", d$sample)
  d$sample_noD <- sub("D$", "", d$sample_sep)
  d$sample_noD <- ifelse(grepl("SE$", d$sample_noD), paste0(d$sample_noD, "D"), d$sample_noD)
  d$sample <- gsub(" ", "", d$sample)
  d
}

# Outlier criteria (per sample, across its pseudoreplicates)
flag_outliers <- function(d) {
  d %>%
    group_by(sample) %>%
    mutate(
      med      = median(copias_ng, na.rm = TRUE),
      MAD      = mad(copias_ng, constant = 1, na.rm = TRUE),
      mz       = 0.6745 * (copias_ng - med) / MAD,
      mad_outlier = abs(mz) > 3.5,
      rel_deviation = abs(copias_ng - mean(copias_ng, na.rm = TRUE)) / mean(copias_ng, na.rm = TRUE)
    ) %>%
    ungroup()
}

# Pseudoreplicates removed as outliers (gene, plate, row in the Excel file)
removed <- read.csv(text = "
gene,plate,row,sample
nirS,1,51,3HT1H
nirS,1,69,7HT3S
nirS,1,70,8HT1S
nirS,1,72,c-
nirS,2,21,3AL1S
nirS,2,36,16AL2H
nirS,2,43,16HT2H
nirS,2,50,10HT1H
nirS,2,59,12AL1H
nirS,2,68,16HT4S
nirS,2,70,7AL1S
nirS,2,89,12HT2S
nirS,3,20,3PL2H
nirS,3,43,3PL1H
nirS,3,54,1PL2S
nirS,3,57,6PL1S
nirS,3,86,3MG1H
nirS,3,93,6PL2H
nirS,4,5,2AG
nirS,4,37,12PL3S
nirS,4,38,16PL3SR
nirS,4,42,7AG
nirS,4,48,18EE2S
nirS,4,63,10OL1S
nirS,4,68,16AG
nirS,4,80,16AG
nirS,4,90,8AG
nirS,5,5,1CT1S
nirS,5,6,12OF1H
nirS,5,17,15CL3S
nirS,5,36,15SEDR
nirS,5,43,4AM1H
nirS,5,51,17CC1HR
nirS,5,71,3SG1SR
nirS,5,82,6SEDR
nirS,5,92,11AM1H
nirS,5,94,6SEDR
nirS,6,11,16HA3H
nirS,6,33,9HA1H
nosZ,1,5,3HT3H
nosZ,1,33,7HT3H
nosZ,1,40,3HT2S
nosZ,1,78,6HT2S
nosZ,1,86,1HT3S
nosZ,2,11,10AL1H
nosZ,2,12,12AL1S
nosZ,2,43,16HT2H
nosZ,2,80,16HT4S
nosZ,2,91,16HT2S
nosZ,2,93,6AL1S
nosZ,3,60,12PL1S
nosZ,3,70,6PL3H
nosZ,4,30,7AG
nosZ,4,76,2AG
nosZ,5,1,3SG1H
nosZ,5,15,1OF1H
nosZ,5,30,11AM1SD
nosZ,5,40,3HP1HD
nosZ,5,63,12OF1SD
nosZ,5,92,3SG1SR
nosZ,6,22,16HA3H
nosZ,6,59,8HT2HR
nosZ,6,83,8HT2HR
")
removed$removed <- TRUE

replicates <- bind_rows(read_plates("nirS"), read_plates("nosZ")) %>%
  clean_names() %>%
  group_by(gene) %>%
  group_modify(~ flag_outliers(.x)) %>%
  ungroup() %>%
  left_join(select(removed, gene, plate, row, removed_sample = sample, removed),
            by = c("gene", "plate", "row")) %>%
  mutate(removed = !is.na(removed))

# Check that the listed rows are the expected samples and meet the MAD criterion
stopifnot(sum(replicates$removed) == nrow(removed),
          with(subset(replicates, removed), all(sample == removed_sample)),
          with(subset(replicates, removed), all(mad_outlier)))

write.csv(select(replicates, gene, plate, row, sample, sample_id = sample_noD, ng_sample,
                 copias_ng, mz, mad_outlier, rel_deviation, removed),
          file.path(tables_dir, "ddpcr_replicates.csv"), row.names = FALSE)

# Mean and SD of the remaining pseudoreplicates. Negative controls are dropped.
# One nirS pseudoreplicate (plate 6, row 22) was labelled "16HA2h" instead of
# "16HA2H" in the QX Manager export. Because of this labelling error it was
# not assigned to sample 16_HA2_H during data processing.
gene_means <- replicates %>%
  filter(!removed, !sample_noD %in% c("c-", "C-", "16_HA2_h")) %>%
  group_by(gene, sample = sample_noD) %>%
  summarise(copias = mean(copias_ng, na.rm = TRUE),
            sd     = sd(copias_ng, na.rm = TRUE),
            .groups = "drop")

nirs_mean <- gene_means %>% filter(gene == "nirS") %>% select(sample, copiasnirs = copias, sdnirs = sd)
nosz_mean <- gene_means %>% filter(gene == "nosZ") %>% select(sample, copiasnosz = copias, sdnosz = sd)

metadata <- read.csv(file.path("data", "ddpcr_sample_metadata.csv"))

curated <- full_join(nirs_mean, nosz_mean, by = "sample") %>%
  mutate(rationirsnosz = copiasnirs / copiasnosz,
         log10nirs = log10(copiasnirs + 1),
         log10nosz = log10(copiasnosz + 1),
         log10ratio = log10nirs / log10nosz) %>%   # nirS/nosZ ratio of Figure 5
  left_join(metadata, by = "sample") %>%
  mutate(class = unname(class_names[clases])) %>%
  arrange(sample)

write.csv(curated, file.path(tables_dir, "ddpcr_curated.csv"), row.names = FALSE)
