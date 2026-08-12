#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

outdir <- "/BLUES/eric/ONT_WGBS/Figure_2"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Reference genome FASTA 
ref_fa <- "/BLUES/eric/ONT_WGBS/refs/hg38.fa"

# methylC format: chr start end context meth_rate strand coverage
files <- c(
  "W3 naive"  = "/BLUES/eric/ONT/MAPQ_bam/ONT_Naive.MAPQ10.primary.noY.onC.methylC.gz",
  "W3 primed" = "/BLUES/eric/ONT/MAPQ_bam/ONT_Primed.MAPQ10.primary.noY.onC.methylC.gz",
  "W3 TSC"    = "/BLUES/eric/ONT/MAPQ_bam/ONT_TSC.MAPQ10.primary.noY.onC.methylC.gz",
  "H9 primed" = "/BLUES/eric/WGBS/CpG/onC_noY/H9-prime_noY_onC.methylC.gz",
  "H9 naive"  = "/BLUES/eric/WGBS/CpG/onC_noY/H9-Naive_noY_onC.methylC.gz",
  "H9 TSC"    = "/BLUES/eric/WGBS/CpG/onC_noY/H9-TSC_noY_onC.methylC.gz"
)

sample_order <- c("W3 naive", "W3 primed", "W3 TSC", "H9 primed", "H9 naive", "H9 TSC")

# ------------------------------------------------------------
# Count total reference CpGs in hg38 (chrY excluded)
# ------------------------------------------------------------
count_reference_cpg_noY <- function(fasta) {
  cmd <- sprintf(
    paste(
      "awk '",
      "/^>/ {chr=substr($1,2); next}",
      "chr!=\"chrY\" {",
      "  seq=toupper($0);",
      "  while (match(seq, /CG/)) {",
      "    c++;",
      "    seq=substr(seq, RSTART+1)",
      "  }",
      "}",
      "END{print c+0}' %s",
      sep = " "
    ),
    shQuote(fasta)
  )
  as.numeric(system(cmd, intern = TRUE)[1])
}

# ------------------------------------------------------------
# Count CpG records in methylC file with coverage >= threshold
# chrY excluded
# ------------------------------------------------------------
count_cpg_ge <- function(methylc_gz, min_cov = 1L, exclude_chr = "chrY") {
  if (!file.exists(methylc_gz)) stop("Missing methylC: ", methylc_gz)

  cmd <- sprintf(
    "zcat %s | awk '$1!=\"%s\" && $1!=\"Y\" && $4==\"CG\" && $7>=%d {c++} END{print c+0}'",
    shQuote(methylc_gz), exclude_chr, as.integer(min_cov)
  )
  as.numeric(system(cmd, intern = TRUE)[1])
}

cat("[1] Counting total reference CpGs (chrY excluded)...\n")
total_ref_cpg <- count_reference_cpg_noY(ref_fa)
cat("Total reference CpGs (chrY excluded):", comma(total_ref_cpg), "\n")

# ------------------------------------------------------------
# Count detected CpGs per sample
# ------------------------------------------------------------
dt <- data.table(
  sample  = names(files),
  methylC = unname(files)
)

cat("[2] Counting CpGs per sample...\n")
dt[, CpG_ge1 := vapply(methylC, count_cpg_ge, FUN.VALUE = numeric(1), min_cov = 1L)]
dt[, CpG_ge10 := vapply(methylC, count_cpg_ge, FUN.VALUE = numeric(1), min_cov = 10L)]

dt[, CpG_1to9 := CpG_ge1 - CpG_ge10]
dt[, CpG_0 := total_ref_cpg - CpG_ge1]

# save wide table
fwrite(
  dt[, .(sample, tech, total_ref_cpg, CpG_ge1, CpG_ge10, CpG_1to9, CpG_0)],
  file.path(outdir, "CpG_coverage_6samples_counts.tsv"),
  sep = "\t"
)

# ------------------------------------------------------------
# Long format for plotting
# ------------------------------------------------------------
plot_dt <- rbindlist(list(
  dt[, .(sample, tech, category = "Uncovered (0×)",        count = CpG_0)],
  dt[, .(sample, tech, category = "Low coverage (1–9×)",   count = CpG_1to9)],
  dt[, .(sample, tech, category = "Covered (≥10×)",        count = CpG_ge10)]
))

plot_dt[, pct := 100 * count / total_ref_cpg]

# plotting order left -> right
plot_dt[, category := factor(
  category,
  levels = c("Uncovered (0×)", "Low coverage (1–9×)", "Covered (≥10×)")
)]

plot_dt[, label := sprintf("%.1f%%", pct)]
plot_dt[pct < 2, label := ""]

# save long table
fwrite(plot_dt, file.path(outdir, "CpG_coverage_6samples_long.tsv"), sep = "\t")

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
cols <- c(
  "Covered (≥10×)"      = "#F87171",
  "Low coverage (1–9×)" = "#E5E7EB",
  "Uncovered (0×)"      = "#9CA3AF"
)

p <- ggplot(plot_dt, aes(x = pct, y = sample, fill = category)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.25) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 3.5
  ) +
  scale_fill_manual(
    values = cols,
    breaks = c("Covered (≥10×)", "Low coverage (1–9×)", "Uncovered (0×)"),
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    title = "CpG coverage categories across the whole genome (chrY excluded)",
    subtitle = "Covered = reference CpGs with ≥10× coverage; low coverage = 1–9×; uncovered = 0×",
    x = "Fraction of reference CpGs (%)",
    y = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "top",
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 11)
  )

ggsave(
  filename = file.path(outdir, "CpG_coverage_6samples.pdf"),
  plot = p,
  width = 12,
  height = 4.8,
  device = cairo_pdf,
  bg = "white"
)

cat("[Done] wrote:\n")
