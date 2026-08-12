#!/usr/bin/env Rscript

# =====================================================
# MAKE FINAL CONSERVATION TSV
# =====================================================

suppressPackageStartupMessages({
  library(data.table)
})

setDTthreads(8)

# =====================================================
# INPUT
# =====================================================

base_dir <- "/BLUES/eric/ONT_WGBS/Figure_3/Jitter_plot_conservation"

dfam_pc_file <- file.path(
  base_dir,
  "Dfam_TE_inf_PC_summary.tsv"
)

dfam_linc_file <- file.path(
  base_dir,
  "Dfam_TE_inf_lincRNA_summary.tsv"
)

rna_file <- "/BLUES/eric/ONT_WGBS/Heatmap/W3_TEsubfamily_RNAseq_summary_meanFPKM_no_rRNA_tRNA.tsv"

meth_file <- "/BLUES/eric/ONT_WGBS/Heatmap/W3_TEsubfamily_DNA_methylation_summary.tsv"

# =====================================================
# OUTPUT
# =====================================================

outfile <- file.path(
  base_dir,
  "TE_conservation_methylation_expression.tsv"
)

# =====================================================
# LOAD Dfam
# =====================================================

pc <- fread(dfam_pc_file)
linc <- fread(dfam_linc_file)

dfam <- rbind(pc, linc, fill = TRUE)

# remove duplicates
dfam <- unique(
  dfam,
  by = c("V7")
)

setnames(dfam,
         c("V7", "V8"),
         c("subfamily", "class"))

dfam <- dfam[, .(
  subfamily,
  class,
  MYA,
  Group
)]

cat("Dfam rows:", nrow(dfam), "\n")

# =====================================================
# LOAD RNA
# =====================================================

rna <- fread(rna_file)

rna <- dcast(
  rna,
  subfamily ~ sample,
  value.var = "mean_fpkm"
)

setnames(
  rna,
  c("Primed", "Naive", "TSC"),
  c(
    "Primed_expression",
    "Naive_expression",
    "TSC_expression"
  )
)

cat("RNA rows:", nrow(rna), "\n")

# =====================================================
# LOAD METHYLATION
# =====================================================

meth <- fread(meth_file)

meth <- dcast(
  meth,
  subfamily ~ sample,
  value.var = "wmean_value"
)

setnames(
  meth,
  c("Primed", "Naive", "TSC"),
  c(
    "Primed_methylation",
    "Naive_methylation",
    "TSC_methylation"
  )
)

cat("Methylation rows:", nrow(meth), "\n")

# =====================================================
# MERGE
# =====================================================

dt <- merge(
  dfam,
  rna,
  by = "subfamily",
  all = FALSE
)

dt <- merge(
  dt,
  meth,
  by = "subfamily",
  all = FALSE
)

# =====================================================
# CLEAN
# =====================================================

dt <- dt[
  is.finite(MYA)
]

# =====================================================
# SAVE
# =====================================================

fwrite(
  dt,
  outfile,
  sep = "\t",
  quote = FALSE
)

cat("Saved:\n")
cat(outfile, "\n")
cat("Final rows:", nrow(dt), "\n")
