#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
})

setDTthreads(8)

# =====================================================
# INPUT / OUTPUT
# =====================================================

indir <- "/BLUES/eric/ONT_WGBS/Heatmap/RNAseq_counts"
outdir <- "/BLUES/eric/ONT_WGBS/Figure_3/MA_Plot/deseq"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

files <- list(
  Primed_Rep1 = file.path(indir, "W3_Primed_Rep1_subFcounts.txt"),
  Primed_Rep2 = file.path(indir, "W3_Primed_Rep2_subFcounts.txt"),
  Naive_Rep1  = file.path(indir, "W3_Naive_Rep1_subFcounts.txt"),
  Naive_Rep2  = file.path(indir, "W3_Naive_Rep2_subFcounts.txt"),
  TSC_Rep1    = file.path(indir, "W3_TSC_Rep1_subFcounts.txt"),
  TSC_Rep2    = file.path(indir, "W3_TSC_Rep2_subFcounts.txt")
)

count_matrix_out <- file.path(outdir, "SQuIRE_tot_counts_matrix.tsv")
sample_info_out  <- file.path(outdir, "SQuIRE_sample_info.tsv")
deseq_out        <- file.path(outdir, "RNAseq_log2FC_DESeq2.tsv")

# =====================================================
# READ SQuIRE COUNTS
# =====================================================

read_squire_counts <- function(f, sample_name) {
  dt <- fread(f)

  if (!"Subfamily:Family:Class" %in% names(dt)) {
    stop("Missing Subfamily:Family:Class in ", f)
  }

  if (!"tot_counts" %in% names(dt)) {
    stop("Missing tot_counts in ", f)
  }

  dt[, subfamily := tstrsplit(
    `Subfamily:Family:Class`,
    ":",
    fixed = TRUE,
    keep = 1L
  )]

  dt[, count := as.numeric(tot_counts)]
  dt[, sample := sample_name]

  dt <- dt[
    !is.na(subfamily) &
      subfamily != "" &
      !is.na(count)
  ]

  dt <- dt[
    ,
    .(count = sum(count, na.rm = TRUE)),
    by = .(subfamily, sample)
  ]

  dt[, count := round(count)]

  dt
}

# =====================================================
# COMBINE COUNTS
# =====================================================

lst <- lapply(
  names(files),
  function(nm) {
    cat("Reading:", nm, "\n")
    read_squire_counts(files[[nm]], nm)
  }
)

counts_long <- rbindlist(
  lst,
  use.names = TRUE,
  fill = TRUE
)

cat("counts_long rows:", nrow(counts_long), "\n")
cat("unique subfamilies:", uniqueN(counts_long$subfamily), "\n")
cat("unique samples:", uniqueN(counts_long$sample), "\n")

# =====================================================
# COUNT MATRIX
# =====================================================

sample_cols <- names(files)

counts_wide <- dcast(
  counts_long,
  subfamily ~ sample,
  value.var = "count",
  fun.aggregate = sum,
  fill = 0
)

setcolorder(
  counts_wide,
  c("subfamily", sample_cols)
)

fwrite(
  counts_wide,
  count_matrix_out,
  sep = "\t",
  quote = FALSE
)

count_mat <- as.matrix(
  counts_wide[, ..sample_cols]
)

rownames(count_mat) <- counts_wide$subfamily
storage.mode(count_mat) <- "integer"

# =====================================================
# SAMPLE INFO
# =====================================================

sample_info <- data.table(
  sample = sample_cols,
  condition = c(
    "Primed", "Primed",
    "Naive", "Naive",
    "TSC", "TSC"
  )
)

sample_info[, condition := factor(
  condition,
  levels = c("Primed", "Naive", "TSC")
)]

fwrite(
  sample_info,
  sample_info_out,
  sep = "\t",
  quote = FALSE
)

coldata <- as.data.frame(sample_info)
rownames(coldata) <- coldata$sample

# =====================================================
# DESEQ2
# =====================================================

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = coldata,
  design = ~ condition
)

dds <- dds[rowSums(counts(dds)) >= 10, ]

dds <- DESeq(dds)

res_NP <- results(
  dds,
  contrast = c("condition", "Naive", "Primed")
)

res_TN <- results(
  dds,
  contrast = c("condition", "TSC", "Naive")
)

res_NP_dt <- as.data.table(res_NP, keep.rownames = "subfamily")
res_TN_dt <- as.data.table(res_TN, keep.rownames = "subfamily")

setnames(
  res_NP_dt,
  old = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"),
  new = c("baseMean_NP", "log2FC_NP", "lfcSE_NP", "stat_NP", "pvalue_NP", "padj_NP")
)

setnames(
  res_TN_dt,
  old = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"),
  new = c("baseMean_TN", "log2FC_TN", "lfcSE_TN", "stat_TN", "pvalue_TN", "padj_TN")
)

out <- merge(
  res_NP_dt,
  res_TN_dt,
  by = "subfamily",
  all = TRUE
)

fwrite(
  out,
  deseq_out,
  sep = "\t",
  quote = FALSE
)

cat("Saved:\n")
cat(count_matrix_out, "\n")
cat(sample_info_out, "\n")
cat(deseq_out, "\n")
cat("Rows:", nrow(out), "\n")