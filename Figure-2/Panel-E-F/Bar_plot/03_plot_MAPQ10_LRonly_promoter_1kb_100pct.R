#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# Input/output:
base_dir <- "/BLUES/eric/ONT_WGBS/Figure_2/Panel_E/outputs"

all100_file <- file.path(
  base_dir,
  "promoters_1000bp_100pct_LRonly_MAPQ10_regions.summary.tsv"
)

pc100_file <- file.path(
  base_dir,
  "protein_coding_promoters_1000bp_100pct_LRonly_MAPQ10_regions.summary.tsv"
)

out_pdf <- file.path(
  base_dir,
  "LRonly_MAPQ10_overlap_100pct_promoter1kb_barplot.pdf"
)

out_png <- file.path(
  base_dir,
  "LRonly_MAPQ10_overlap_100pct_promoter1kb_barplot.png"
)


# Read headerless gene tables

column_names <- c(
  "chr",
  "start",
  "end",
  "gene_id",
  "gene_name",
  "gene_type",
  "strand"
)

read_summary <- function(file) {
  if (!file.exists(file)) {
    stop("Missing input file: ", file)
  }

  dt <- fread(
    file,
    header = FALSE,
    sep = "\t"
  )

  if (ncol(dt) != 7) {
    stop(
      "Expected 7 columns but found ",
      ncol(dt),
      " in: ",
      file
    )
  }

  setnames(dt, column_names)
  dt
}

all_dt <- read_summary(all100_file)
pc_dt  <- read_summary(pc100_file)

# Count unique genes, not simply rows
all_val <- uniqueN(all_dt$gene_id)
pc_val  <- uniqueN(pc_dt$gene_id)

if (pc_val > all_val) {
  stop("Protein-coding count exceeds all-gene count.")
}

cat("All-gene promoter rows:            ", nrow(all_dt), "\n")
cat("Unique all-gene IDs:               ", all_val, "\n")
cat("Protein-coding promoter rows:      ", nrow(pc_dt), "\n")
cat("Unique protein-coding gene IDs:    ", pc_val, "\n")

# Check for duplicate gene IDs
all_duplicates <- nrow(all_dt) - all_val
pc_duplicates  <- nrow(pc_dt) - pc_val

cat("Duplicate all-gene records:        ", all_duplicates, "\n")
cat("Duplicate protein-coding records:  ", pc_duplicates, "\n")

# Plot table

plot_dt <- data.table(
  group = factor(
    c("All genes", "Protein-coding"),
    levels = c("All genes", "Protein-coding")
  ),
  count = c(all_val, pc_val)
)

plot_dt[, label := comma(count)]

# Colors

fill_cols <- c(
  "All genes" = "#377EB8",
  "Protein-coding" = "#E69F00"
)

# Plot

p <- ggplot(
  plot_dt,
  aes(x = group, y = count, fill = group)
) +
  geom_col(
    width = 0.60,
    color = "black",
    linewidth = 0.35
  ) +
  geom_text(
    aes(label = label),
    vjust = -0.30,
    size = 5
  ) +
  scale_fill_manual(values = fill_cols) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.10))
  ) +
  labs(
    title = "Genes overlapping MAPQ10 LR-only regions",
    subtitle = "Promoter ±1 kb • 100% contained within LR-only regions",
    x = NULL,
    y = "Number of unique genes"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(
      face = "bold",
      size = 17
    ),
    plot.subtitle = element_text(size = 12),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(
      face = "bold",
      size = 13,
      color = "black"
    ),
    axis.text.y = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

ggsave(
  out_pdf,
  p,
  width = 5.2,
  height = 6,
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  out_png,
  p,
  width = 5.2,
  height = 6,
  dpi = 300,
  bg = "white"
)

cat("\nWritten:\n")
cat(out_pdf, "\n")
cat(out_png, "\n")
