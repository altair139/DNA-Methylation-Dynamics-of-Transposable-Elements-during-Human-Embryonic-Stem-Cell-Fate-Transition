#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

outdir <- "/BLUES/eric/ONT_WGBS/Figure_2"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# methylC format: chr start end context meth_rate strand coverage
files <- c(
  "W3 naive"  = "/BLUES/eric/ONT/MAPQ_bam/ONT_Naive.MAPQ10.primary.noY.onC.methylC.gz",
  "W3 primed" = "/BLUES/eric/ONT/MAPQ_bam/ONT_Primed.MAPQ10.primary.noY.onC.methylC.gz",
  "W3 TSC"    = "/BLUES/eric/ONT/MAPQ_bam/ONT_TSC.MAPQ10.primary.noY.onC.methylC.gz",
  "H9 primed" = "/BLUES/eric/WGBS/CpG/H9-primed_onC.methylC.gz",
  "H9 naive"  = "/BLUES/eric/WGBS/CpG/H9-Naive_onC.methylC.gz",
  "H9 TSC"    = "/BLUES/eric/WGBS/CpG/H9-TSC_onC.methylC.gz"
)

# Order 
sample_order <- c("W3 naive", "W3 primed", "W3 TSC", "H9 primed", "H9 naive", "H9 TSC")

# Conditions
min_cov <- 1L
exclude_chr <- "chrY"  # whole genome except chrY

count_cpg_whole_genome_noY <- function(methylc_gz, min_cov = 1L, exclude_chr = "chrY") {
  if (!file.exists(methylc_gz)) stop("Missing methylC: ", methylc_gz)

  # Count rows that are CG, coverage>=min_cov, and chr != chrY.
  # NOTE: This counts "CpG records" in the methylC file (position-level).
  cmd <- sprintf(
    "zcat %s | awk '$1!=\"%s\" && $4==\"CG\" && $7>=%d {c++} END{print c+0}'",
    shQuote(methylc_gz), exclude_chr, as.integer(min_cov)
  )
  out <- system(cmd, intern = TRUE)
  as.numeric(out[1])
}

counts <- data.table(
  sample  = names(files),
  methylC = unname(files)
)

counts[, CpG_ge1_wholeGenome_noY := vapply(
  methylC, count_cpg_whole_genome_noY,
  FUN.VALUE = numeric(1),
  min_cov = min_cov,
  exclude_chr = exclude_chr
)]

counts[, tech := ifelse(grepl("^W3", sample), "LR (ONT)", "SR (WGBS)")]
counts[, sample := factor(sample, levels = sample_order)]

# Save table
fwrite(
  counts[, .(sample, tech, CpG_ge1_wholeGenome_noY)],
  file.path(outdir, "WholeGenome_noY_CpGcounts_ge1_Q10.tsv"),
  sep = "\t"
)

# Plot
p <- ggplot(counts, aes(x = sample, y = CpG_ge1_wholeGenome_noY, fill = tech)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.4) +
  geom_text(aes(label = comma(CpG_ge1_wholeGenome_noY)),
            vjust = -0.4, size = 3.6) +
  scale_y_continuous(
    labels = comma,
    breaks = c(0, 1e7, 2e7, 3e7),  
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    title = "Detected CpGs (≥1×) across the whole genome (chrY excluded)",
    x = NULL,
    y = "Detected CpG count (coverage ≥1×)",
    fill = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.6),
    axis.line.y = element_line(color = "black", linewidth = 0.6),
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey80", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "top"
  )

ggsave(
  filename = file.path(outdir, "Fig2A_left_WholeGenome_noY_CpGcounts_Q10.pdf"),
  plot = p, width = 10.5, height = 4.8, device = cairo_pdf
)

