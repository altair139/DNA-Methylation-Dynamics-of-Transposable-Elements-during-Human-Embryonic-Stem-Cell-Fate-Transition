#!/usr/bin/env Rscript
#Usage: Rscript plot_Fig2C_pies.R /path/to/Panel_C
suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

panel_c <- args[1]
setDTthreads(threads = 8)

wo_file <- file.path(
  panel_c,
  "outputs",
  "LRonly_x_TE.allchr.wo.tsv"
)

lr_file <- file.path(
  panel_c,
  "intermediate",
  "LR_only_windows_100bp.sorted.bed"
)

out_pdf <- file.path(
  panel_c,
  "outputs",
  "Fig2C_TE_pies_top10_legend.pdf"
)

# ---- Function: calculate unique genomic length ----
# This merges overlapping or adjacent intervals so that each base
# is counted only once.
union_bp <- function(intervals) {
  intervals <- intervals[
    !is.na(chr) &
    !is.na(start) &
    !is.na(end) &
    end > start
  ]

  intervals <- unique(intervals[, .(
    chr   = as.character(chr),
    start = as.numeric(start),
    end   = as.numeric(end)
  )])

  setorder(intervals, chr, start, end)

  intervals[, previous_max_end :=
    shift(cummax(end), fill = -Inf),
    by = chr
  ]

  intervals[, interval_group :=
    cumsum(start > previous_max_end),
    by = chr
  ]

  merged <- intervals[, .(
    start = min(start),
    end   = max(end)
  ), by = .(chr, interval_group)]

  merged[, sum(end - start)]
}

# ---- Load the complete LR-only window universe ----
lr <- fread(
  lr_file,
  sep = "\t",
  header = FALSE,
  select = 1:3,
  col.names = c("chr", "start", "end")
)

# This should be 122,042,900 bp
total_bp <- union_bp(lr)

# ---- Load LR-window × TE overlaps ----
dt <- fread(
  wo_file,
  sep = "\t",
  header = FALSE
)

setnames(dt, c(
  "A_chr", "A_start", "A_end",
  "B_chr", "B_start", "B_end", "B_strand",
  "subfamily", "class", "family", "overlap_bp"
))

dt[, `:=`(
  A_start   = as.numeric(A_start),
  A_end     = as.numeric(A_end),
  B_start   = as.numeric(B_start),
  B_end     = as.numeric(B_end),
  overlap_bp = as.numeric(overlap_bp)
)]

# Remove exact duplicate intersection records, if any
dt <- unique(dt)

# ---- Normalize TE labels ----
for (col in c("subfamily", "class", "family")) {
  dt[, (col) := fifelse(
    is.na(get(col)) |
      get(col) == "" |
      get(col) == "NA",
    "Unclassified",
    as.character(get(col))
  )]
}

# ---- Calculate unique TE-overlapping bp ----
# Clip each TE to the portion that lies inside its LR-only window.
te_intervals <- dt[, .(
  chr   = A_chr,
  start = pmax(A_start, B_start),
  end   = pmin(A_end, B_end)
)]

# Merge overlapping TE intervals before summing.
te_bp <- union_bp(te_intervals)

if (te_bp > total_bp) {
  stop(
    "Unique TE bp exceeds total LR-only bp. ",
    "The overlap file and LR-window file may represent different datasets."
  )
}

non_te_bp <- total_bp - te_bp

# Raw sum is reported only as a diagnostic.
raw_te_bp <- dt[, sum(overlap_bp, na.rm = TRUE)]

message("Total unique LR-only window bp: ", total_bp)
message("Raw summed TE overlap bp:       ", raw_te_bp)
message("Unique TE-overlapping bp:       ", te_bp)
message("Unique non-TE bp:               ", non_te_bp)
message(
  "TE percentage:                    ",
  sprintf("%.3f%%", 100 * te_bp / total_bp)
)
message(
  "Non-TE percentage:                ",
  sprintf("%.3f%%", 100 * non_te_bp / total_bp)
)

# ---- Top-N function ----
topN_by_bp <- function(DT, col, topN = 10) {
  x <- DT[
    ,
    .(bp = sum(overlap_bp, na.rm = TRUE)),
    by = c(col)
  ][order(-bp)]

  if (nrow(x) > topN) {
    other <- copy(x[1])

    other[, (col) := "Other"]
    other[, bp := sum(x[(topN + 1):.N, bp])]

    x <- rbindlist(
      list(
        x[1:topN],
        other
      ),
      use.names = TRUE
    )
  }

  x
}

make_palette <- function(n) {
  grDevices::hcl.colors(n, palette = "Dark 3")
}

plot_pie_with_legend <- function(values, labels, title) {
  pct <- 100 * values / sum(values)
  legend_text <- sprintf("%s (%.2f%%)", labels, pct)

  cols <- make_palette(length(values))
  cols[labels == "Other"] <- "grey50"

  par(mar = c(2, 2, 3.5, 10))

  pie(
    values,
    labels = NA,
    col = cols,
    border = "white",
    main = title
  )

  legend(
    "right",
    inset = c(-0.35, 0),
    xpd = TRUE,
    legend = legend_text,
    fill = cols,
    bty = "n",
    cex = 0.9
  )
}

# ---- Pie 1: unique TE bp versus unique non-TE bp ----
pie1_labels <- c("TE", "Non-TE")
pie1_vals <- c(te_bp, non_te_bp)
pie1_cols <- c("#1b9e77", "#d95f02")

# ---- Top-10 TE summaries ----
pie2 <- topN_by_bp(dt, "class")
pie3 <- topN_by_bp(dt, "family")
pie4 <- topN_by_bp(dt, "subfamily")

# ---- Plot ----
pdf(out_pdf, width = 13, height = 8.5)

par(
  mfrow = c(2, 2),
  cex = 1
)

par(mar = c(2, 2, 3.5, 2))

pie(
  pie1_vals,
  labels = sprintf(
    "%s\n%.2f%%",
    pie1_labels,
    100 * pie1_vals / sum(pie1_vals)
  ),
  col = pie1_cols,
  border = "white",
  main = "TE composition of LR-only 100-bp windows"
)

plot_pie_with_legend(
  pie2$bp,
  pie2$class,
  "Top 10 TE classes by overlap length"
)

plot_pie_with_legend(
  pie3$bp,
  pie3$family,
  "Top 10 TE families by overlap length"
)

plot_pie_with_legend(
  pie4$bp,
  pie4$subfamily,
  "Top 10 TE subfamilies by overlap length"
)

dev.off()

message("Wrote: ", out_pdf)