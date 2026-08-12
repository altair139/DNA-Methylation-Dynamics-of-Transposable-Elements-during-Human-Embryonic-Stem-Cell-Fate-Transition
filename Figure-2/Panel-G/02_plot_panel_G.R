#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: plot_Fig2D_PI_fractionCaptured.R /path/to/Panel_D")

panel_d <- args[1]
setDTthreads(threads = 8)

infile <- file.path(
  panel_d,
  "outputs",
  "Fig2D_top10_fractionCaptured.tsv"
)

toplist <- file.path(
  panel_d,
  "inputs",
  "top10_subfamilies_from_pie.txt"
)

outfile <- file.path(
  panel_d,
  "outputs",
  "Fig2D_dotplot_top10_fractionCaptured.pdf"
)

dt <- fread(infile, header = FALSE)
setnames(dt, c("subfamily","class","family","copies_in_LRonly","copies_genome","pct_captured"))

ord <- fread(toplist, header = FALSE)[[1]]
dt <- dt[subfamily %in% ord]
dt[, subfamily := factor(subfamily, levels = rev(ord))]

p <- ggplot(dt, aes(x = pct_captured, y = subfamily)) +
  geom_point(
    aes(size = copies_in_LRonly, color = family),
    alpha = 0.9
  ) +
  scale_size_continuous(
    name = "Copies in LR-only",
    breaks = c(1000, 3000, 5000, 10000),
    labels = scales::comma,
    range = c(3, 12)
  ) +
  labs(
    x = "Percent of genome-wide copies captured in LR-only region (%)",
    y = "TE subfamily (bp-ranked top 10)",
    title = "Top 10 TE subfamilies: fraction of genome-wide copies captured in LR-only region",
    subtitle = "X = (# copies in LR-only) / (total # copies in genome)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )
ggsave(outfile, p, width = 12, height = 6)

