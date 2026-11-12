# =============================================================================
# 01. 16S rRNA amplicon processing with DADA2
# -----------------------------------------------------------------------------
# Manuscript section: Methods, "Sequencing and data analysis".
# Raw reads: NCBI SRA, BioProject PRJNA1454003 (download them into data/raw/).
# Outputs used in the next step: results/sequence_processing/otu_table_picrust2.tsv
# and otu_sequences.fasta (inputs of PICRUSt2, script 02).
#
# Samples were sequenced in two runs. Error rates are run-specific, so steps
# 1-3 are done separately for each run; the sequence tables of both runs are
# merged in step 4.
#
# From raw paired-end FASTQ files to a rarefied 97 % OTU table:
#   1. Quality inspection
#   2. Filtering and trimming
#   3. Error learning, denoising and read merging
#   4. Chimera removal and length filtering
#   5. Taxonomic assignment (SILVA)
#   6. Removal of non-bacterial ASVs (mitochondria, chloroplasts, archaea)
#   7. Rarefaction
#   8. Clustering of ASVs into 97 % OTUs
#   9. OTU tables
# =============================================================================

library(dada2)
library(tidyverse)
library(vegan)
library(DECIPHER)


# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------
# Paths are relative to the project folder (open the .Rproj file in RStudio).
input_dir   <- file.path("data", "raw")   # raw FASTQ files: *_1.fastq.gz and *_2.fastq.gz
output_dir  <- file.path("results", "sequence_processing")   # tables and R objects
figures_dir <- file.path(output_dir, "figures")
silva_db    <- file.path("data", "reference", "silva_nr99_v138.1_wSpecies_train_set.fa.gz")

# Sequencing run of each sample: columns "sample" and "run"
runs_file   <- file.path("data", "sequencing_runs.csv")

# Truncate reads where quality starts to drop (see the quality profiles)
trunc_len <- c(185, 175)   # forward, reverse
max_ee    <- c(2, 2)       # maximum expected errors: forward, reverse

# Expected amplicon length (see the ASV length distribution)
asv_length_range <- 251:256

# Rarefaction depth and samples that do not reach it (see the rarefaction curves)
rarefaction_depth <- 43100
samples_to_remove <- c("X1", "S8", "L13", "R7")
rarefaction_reps  <- 100

# OTU clustering threshold
otu_identity <- 0.97   # 97 % similarity

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# 1. Quality inspection
# =============================================================================
fnFs <- sort(list.files(input_dir, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(input_dir, pattern = "_2.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
names(fnFs) <- sample.names
names(fnRs) <- sample.names

# Sequencing run of each sample
runs <- read.csv(runs_file)
stopifnot(all(sample.names %in% runs$sample))
run_names <- unique(runs$run)
samples_in_run <- function(run) sample.names[sample.names %in% runs$sample[runs$run == run]]

for (run in run_names) {
  s <- samples_in_run(run)

  # Inspect where quality starts to drop to choose trunc_len (four samples per plot)
  for (i in seq(1, length(s), by = 4)) {
    idx <- s[i:min(i + 3, length(s))]
    print(plotQualityProfile(fnFs[idx]))
    print(plotQualityProfile(fnRs[idx]))
  }

  # Save the quality profiles of the first samples of the run
  first_samples <- s[seq_len(min(9, length(s)))]
  qplot_f <- plotQualityProfile(fnFs[first_samples]) + ggtitle(paste(run, "- Forward reads"))
  qplot_r <- plotQualityProfile(fnRs[first_samples]) + ggtitle(paste(run, "- Reverse reads"))
  ggsave(file.path(figures_dir, paste0("quality_forward_", run, ".pdf")), qplot_f, width = 410, height = 300, units = "mm")
  ggsave(file.path(figures_dir, paste0("quality_reverse_", run, ".pdf")), qplot_r, width = 410, height = 300, units = "mm")
}


# =============================================================================
# 2-3. Filtering, trimming, error learning, denoising and read merging (per run)
# =============================================================================
filtFs <- file.path(input_dir, "filtered", paste0(sample.names, "_F_filt.fq.gz"))
filtRs <- file.path(input_dir, "filtered", paste0(sample.names, "_R_filt.fq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

getN <- function(x) sum(getUniques(x))
seqtabs   <- list()   # sequence table of each run
track_run <- list()   # reads kept in each step, per run

for (run in run_names) {
  s <- samples_in_run(run)

  # 2. Filtering and trimming
  # truncLen: reads are truncated to this length; shorter reads are discarded
  # maxN:     DADA2 does not allow ambiguous bases (N)
  # maxEE:    reads with more expected errors are discarded, EE = sum(10^(-Q/10))
  # truncQ:   reads are truncated at the first base with quality <= truncQ
  # rm.phix:  reads matching the phiX genome are discarded
  out <- filterAndTrim(fnFs[s], filtFs[s], fnRs[s], filtRs[s],
                       truncLen = trunc_len, maxN = 0, maxEE = max_ee, truncQ = 2,
                       rm.phix = TRUE, compress = TRUE, multithread = TRUE)
  print(head(out))

  # 3. Error learning (run-specific)
  errF <- learnErrors(filtFs[s], multithread = TRUE, randomize = TRUE, verbose = TRUE)
  errR <- learnErrors(filtRs[s], multithread = TRUE, randomize = TRUE, verbose = TRUE)

  ggsave(file.path(figures_dir, paste0("errors_forward_", run, ".pdf")),
         plotErrors(errF, nominalQ = TRUE), width = 9, height = 8)
  ggsave(file.path(figures_dir, paste0("errors_reverse_", run, ".pdf")),
         plotErrors(errR, nominalQ = TRUE), width = 9, height = 8)

  # Dereplication collapses identical reads and averages their quality scores
  derepFs <- derepFastq(filtFs[s], verbose = TRUE)
  derepRs <- derepFastq(filtRs[s], verbose = TRUE)
  names(derepFs) <- s
  names(derepRs) <- s

  # Sample inference
  dadaFs <- dada(derepFs, err = errF, multithread = TRUE, verbose = TRUE)
  dadaRs <- dada(derepRs, err = errR, multithread = TRUE, verbose = TRUE)

  # Merge read pairs; pairs that do not overlap enough or have mismatches
  # in the overlap region are rejected
  mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

  seqtabs[[run]] <- makeSequenceTable(mergers)
  cat(run, "- samples:", nrow(seqtabs[[run]]), "- ASVs:", ncol(seqtabs[[run]]), "\n")
  saveRDS(seqtabs[[run]], file.path(output_dir, paste0("seqtab_with_chimeras_", run, ".rds")))

  track_run[[run]] <- data.frame(
    sample    = s,
    run       = run,
    input     = out[, "reads.in"],
    filtered  = out[, "reads.out"],
    denoisedF = sapply(dadaFs, getN),
    denoisedR = sapply(dadaRs, getN),
    merged    = sapply(mergers, getN)
  )
}


# =============================================================================
# 4. Chimera removal and length filtering
# =============================================================================
# The tables of both runs are merged and chimeras are removed from all of them
# together.
st.all <- if (length(seqtabs) > 1) mergeSequenceTables(tables = seqtabs) else seqtabs[[1]]

seqtab.nochim <- removeBimeraDenovo(st.all, method = "consensus", multithread = TRUE, verbose = TRUE)
cat(sprintf("%.2f%% of reads kept after chimera removal\n",
            100 * sum(seqtab.nochim) / sum(st.all)))

# Keep only ASVs within the expected amplicon length
length_dist <- table(nchar(getSequences(seqtab.nochim)))
length_dist
plot(length_dist, xlab = "ASV length (bp)", ylab = "Number of ASVs")

seqtab.len <- seqtab.nochim[, nchar(colnames(seqtab.nochim)) %in% asv_length_range]
table(nchar(getSequences(seqtab.len)))
cat(sprintf("%.2f%% of reads kept after length filtering\n",
            100 * sum(seqtab.len) / sum(seqtab.nochim)))

# Collapse ASVs that are identical except for length differences
seqtab.final <- collapseNoMismatch(seqtab.len, minOverlap = 50, verbose = TRUE)
cat("ASVs after collapsing:", ncol(seqtab.final), "\n")

# Track reads through the pipeline
track <- bind_rows(track_run) %>%
  mutate(nonchim         = rowSums(seqtab.nochim)[sample],
         length_filtered = rowSums(seqtab.final)[sample],
         kept_total      = round(length_filtered / input, 2))

head(track)
summary(track$kept_total)
write_excel_csv(track, file.path(output_dir, "read_tracking.csv"))

saveRDS(seqtab.final, file.path(output_dir, "seqtab_final.rds"))

# Export ASV sequences. Note: these IDs follow abundance order and are not the
# same as the ASV names assigned in step 6.
uniquesToFasta(seqtab.final, file.path(output_dir, "asv_sequences.fasta"),
               ids = paste0("asv", seq_len(ncol(seqtab.final)), ";size=", colSums(seqtab.final)))


# =============================================================================
# 5. Taxonomic assignment (SILVA)
# =============================================================================
set.seed(42)   # assignTaxonomy uses random bootstrapping
# tryRC = TRUE also tries the reverse complement of each sequence
taxa <- assignTaxonomy(seqtab.final, silva_db, multithread = TRUE, tryRC = TRUE, minBoot = 80)

saveRDS(taxa, file.path(output_dir, "taxonomy.rds"))
write_excel_csv(as.data.frame(taxa) %>% rownames_to_column("asv"),
                file.path(output_dir, "taxonomy.csv"))


# =============================================================================
# 6. Removal of non-bacterial ASVs
# =============================================================================
asv_counts <- as.data.frame(t(seqtab.final)) %>% rownames_to_column("asv")
tax_df     <- as.data.frame(taxa) %>% rownames_to_column("asv")

# merge() sorts rows by sequence; ASV names below follow this order
seq_tax_tab <- merge(asv_counts, tax_df, by = "asv")

clean_tab <- seq_tax_tab %>%
  filter(!is.na(Kingdom),
         Kingdom != "Archaea",
         !Order %in% "Chloroplast",
         !Family %in% "Mitochondria")

clean_tab$asvnames <- paste0("ASV", seq_len(nrow(clean_tab)))
rownames(clean_tab) <- clean_tab$asvnames

annotation_cols <- c("asv", colnames(taxa), "asvnames")
sample_cols     <- setdiff(colnames(clean_tab), annotation_cols)

tax_table <- clean_tab[, annotation_cols]              # sequence and taxonomy of each ASV
asv_table <- as.matrix(clean_tab[, sample_cols])       # read counts: ASVs x samples

write.table(clean_tab, file.path(output_dir, "asv_taxonomy_table_clean.txt"))
write.table(tax_table, file.path(output_dir, "taxonomy_table_clean.txt"))


# =============================================================================
# 7. Rarefaction
# =============================================================================
samples_x_asvs <- t(asv_table)

reads_per_sample <- sort(rowSums(samples_x_asvs), decreasing = TRUE)
reads_per_sample
write_excel_csv(enframe(reads_per_sample, name = "sample", value = "reads"),
                file.path(output_dir, "reads_per_sample.csv"))

# Rarefaction curves help choose rarefaction_depth and samples_to_remove
observed <- specnumber(samples_x_asvs)
raremax  <- min(rowSums(samples_x_asvs))
rarefied_richness <- rarefy(samples_x_asvs, raremax)

png(file.path(figures_dir, "observed_vs_rarefied_asvs.png"),
    width = 2644, height = 1566, units = "px", res = 300)
plot(observed, rarefied_richness,
     xlab = "Observed number of ASVs", ylab = "Rarefied number of ASVs")
abline(0, 1)
dev.off()

png(file.path(figures_dir, "rarefaction_curves.png"),
    width = 2644, height = 1566, units = "px", res = 300)
rarecurve(samples_x_asvs, step = 20, sample = raremax, col = "black", cex = 0.6,
          xlab = "Reads", ylab = "Observed number of ASVs")
dev.off()

# Remove samples below the rarefaction depth.
# data.frame() makes sample names syntactic (e.g. "1" becomes "X1").
tab <- data.frame(asv_table)
rar_input <- t(tab[, !colnames(tab) %in% samples_to_remove])
cat("Minimum reads in the remaining samples:", min(rowSums(rar_input)), "\n")

# Rarefy many times and average the results to reduce random subsampling effects
set.seed(41)
rarefied <- lapply(seq_len(rarefaction_reps), function(i) rrarefy(rar_input, rarefaction_depth))
mean_rarefied <- as.data.frame(round(apply(simplify2array(rarefied), c(1, 2), mean)))

# Because of rounding, samples are close to but not exactly at rarefaction_depth
rowSums(mean_rarefied)

# Drop ASVs that are absent from every sample after rarefaction
mean_rarefied <- mean_rarefied[, colSums(mean_rarefied) > 0]
dim(mean_rarefied)

write.table(mean_rarefied, file.path(output_dir, "asv_table_rarefied.txt"))


# =============================================================================
# 8. Clustering of ASVs into 97 % OTUs
# =============================================================================
# ASVs of the rarefied table are clustered with DECIPHER. Each OTU is
# represented by its most abundant ASV, which provides the OTU sequence and
# taxonomy.
asvs_x_samples <- t(mean_rarefied)
asv_seqs <- tax_table[rownames(asvs_x_samples), "asv"]

# method = "shortest": distances are relative to the shorter sequence of each
# pair. See the DECIPHER vignette "Clustering sequences" for details.
clusters <- Clusterize(DNAStringSet(asv_seqs), method = "shortest",
                       cutoff = 1 - otu_identity, processors = NULL)
cat(length(unique(clusters$cluster)), "OTUs created from", length(asv_seqs), "ASVs\n")

asv_otu <- tibble(asvnames = rownames(asvs_x_samples),
                  cluster  = clusters$cluster,
                  reads    = rowSums(asvs_x_samples)) %>%
  group_by(cluster) %>%
  mutate(representative = asvnames[which.max(reads)]) %>%   # first one if tied
  ungroup()

# Name OTUs by decreasing abundance (OTU1 is the most abundant)
otu_names <- asv_otu %>%
  group_by(representative) %>%
  summarise(reads = sum(reads)) %>%
  arrange(desc(reads)) %>%
  mutate(otunames = paste0("OTU", row_number()))

asv_otu <- left_join(asv_otu, dplyr::select(otu_names, representative, otunames),
                     by = "representative")


# =============================================================================
# 9. OTU tables
# =============================================================================
# Sum the reads of all ASVs in each OTU: OTUs x samples
otu_table <- rowsum(asvs_x_samples, asv_otu$otunames)[otu_names$otunames, ]

# Taxonomy and sequence of each OTU, taken from its representative ASV
otu_taxonomy_seqs <- otu_names %>%
  left_join(tax_table, by = c("representative" = "asvnames")) %>%
  dplyr::select(otunames, representative, sequence = asv, all_of(colnames(taxa)))

# Working table without sequences; sequences are kept in separate files
otu_taxonomy <- otu_taxonomy_seqs %>%
  dplyr::select(-representative, -sequence) %>%
  column_to_rownames("otunames")

write.table(otu_table, file.path(output_dir, "otu_table.txt"))
write.table(otu_taxonomy, file.path(output_dir, "otu_taxonomy.txt"))
write.table(otu_taxonomy_seqs, file.path(output_dir, "otu_taxonomy_sequences.txt"),
            row.names = FALSE)
write.table(dplyr::select(asv_otu, asvnames, otunames, reads),
            file.path(output_dir, "asv_otu_correspondence.txt"), row.names = FALSE)
writeXStringSet(DNAStringSet(setNames(otu_taxonomy_seqs$sequence, otu_taxonomy_seqs$otunames)),
                file.path(output_dir, "otu_sequences.fasta"))

# PICRUSt2 input: tab-separated OTU table (OTUs x samples) with the
# "#OTU ID" header expected by biom-format. Script 02 converts it to an
# HDF5 .biom file.
write.table(data.frame("#OTU ID" = rownames(otu_table), otu_table, check.names = FALSE),
            file.path(output_dir, "otu_table_picrust2.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
