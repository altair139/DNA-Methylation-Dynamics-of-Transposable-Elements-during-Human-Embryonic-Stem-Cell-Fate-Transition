#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

setDTthreads(8)

# =====================================================
# INPUT / OUTPUT
# =====================================================

infile <- "/BLUES/eric/ONT_WGBS/Figure_3/Jitter_plot_conservation/TE_conservation_methylation_expression.tsv"

outdir <- "/BLUES/eric/ONT_WGBS/Figure_3/Jitter_plot_conservation"

outfile <- file.path(
  outdir,
  "TE_age_vs_methylation_difference.pdf"
)

# =====================================================
# LOAD
# =====================================================

dt <- fread(infile)

cat("Loaded rows:", nrow(dt), "\n")

# =====================================================
# CALCULATE METHYLATION DIFFERENCES
# =====================================================

dt[, NP_methylation_difference := Naive_methylation - Primed_methylation]
dt[, TN_methylation_difference := TSC_methylation - Naive_methylation]

# =====================================================
# COLORS / SHAPES
# =====================================================

group_colors <- c(
  "Amniota"  = "#F8766D",
  "Mammalia" = "#00BFC4",
  "Primate"  = "#C77CFF",
  "Homo"     = "#7CAE00"
)

shape_vals <- c(
  "DNA"  = 17,
  "LINE" = 15,
  "LTR"  = 16,
  "SINE" = 18
)

# =====================================================
# LABEL OUTLIERS
# =====================================================

np_high <- quantile(dt$NP_methylation_difference, 0.99, na.rm = TRUE)
np_low  <- quantile(dt$NP_methylation_difference, 0.01, na.rm = TRUE)

tn_high <- quantile(dt$TN_methylation_difference, 0.99, na.rm = TRUE)
tn_low  <- quantile(dt$TN_methylation_difference, 0.01, na.rm = TRUE)

label_np <- dt[
  Group == "Homo" |
    NP_methylation_difference >= np_high |
    NP_methylation_difference <= np_low
]

label_tn <- dt[
  Group == "Homo" |
    TN_methylation_difference >= tn_high |
    TN_methylation_difference <= tn_low
]

all_labels <- unique(
  rbind(label_np, label_tn, fill = TRUE),
  by = "subfamily"
)

cat("Labels:", nrow(all_labels), "\n")

# =====================================================
# THEME
# =====================================================

theme_te <- theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 16
    ),
    axis.title = element_text(
      face = "bold",
      size = 14
    ),
    axis.text = element_text(size = 12),
    legend.title = element_text(
      face = "bold",
      size = 13
    ),
    legend.text = element_text(size = 11)
  )

# =====================================================
# N - P PANEL
# =====================================================

p_np <- ggplot(
  dt,
  aes(
    x = (-1) * MYA,
    y = NP_methylation_difference,
    color = Group,
    shape = class
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey40"
  ) +
  geom_jitter(
    width = 2,
    height = 0,
    alpha = 0.8,
    size = 2
  ) +
  geom_text_repel(
    data = all_labels,
    aes(label = subfamily),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.2,
    segment.size = 0.3,
    min.segment.length = 0
  ) +
  scale_color_manual(values = group_colors) +
  scale_shape_manual(values = shape_vals) +
  labs(
    title = "Naive - Primed",
    x = "Insertion time of TE (million years)",
    y = "DNA methylation difference"
  ) +
  theme_te

# =====================================================
# T - N PANEL
# =====================================================

p_tn <- ggplot(
  dt,
  aes(
    x = (-1) * MYA,
    y = TN_methylation_difference,
    color = Group,
    shape = class
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey40"
  ) +
  geom_jitter(
    width = 2,
    height = 0,
    alpha = 0.8,
    size = 2
  ) +
  geom_text_repel(
    data = all_labels,
    aes(label = subfamily),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.2,
    segment.size = 0.3,
    min.segment.length = 0
  ) +
  scale_color_manual(values = group_colors) +
  scale_shape_manual(values = shape_vals) +
  labs(
    title = "TSC - Naive",
    x = "Insertion time of TE (million years)",
    y = "DNA methylation difference"
  ) +
  theme_te

combined <- p_np + p_tn +
  plot_annotation(
    title = "TE evolutionary age vs DNA methylation difference",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 18
      )
    )
  )

ggsave(
  outfile,
  combined,
  width = 14,
  height = 6
)

cat("Saved:\n")
cat(outfile, "\n")