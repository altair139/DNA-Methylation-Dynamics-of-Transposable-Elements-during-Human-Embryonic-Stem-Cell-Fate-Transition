#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

# Paths
input_tsv <- "/BLUES/eric/ONT_WGBS/Figure_2/Panel_B/B2/B2_mappability_violin.tsv"
output_pdf <- "/BLUES/eric/ONT_WGBS/Figure_2/Panel_B/B2/Fig2B_mappability_violin.pdf"

# Read
dt <- fread(input_tsv)

# Clean + factor order
dt2 <- dt %>%
  filter(!is.na(mappability)) %>%   # drop NA only (keep 0)
  mutate(
    # standardize category text 
    category = gsub(">=1x", "≥1x", category, fixed = TRUE),
    category = gsub(">=1x", "≥1x", category),

    category = factor(
      category,
      levels = c("Detected (≥1x)", "Undetected (<1x)")
    ),
    tech = factor(tech, levels = c("WGBS", "ONT")),
    state = factor(state, levels = c("Primed", "Naive", "TSC"))
  )

# Colors by tech
tech_cols <- c("WGBS" = "#2B6CB0", "ONT" = "#ED8936")

# Plot:
# - violins by tech inside each Detected/Undetected
# - boxplots inside
pd <- position_dodge(width = 0.85)

p <- ggplot(dt2, aes(x = category, y = mappability, fill = tech)) +
  geom_violin(
    position = pd,
    trim = TRUE,
    scale = "width",
    color = "black",
    linewidth = 0.25
  ) +
  geom_boxplot(
    position = pd,
    width = 0.15,
    outlier.shape = NA,
    alpha = 0.55,
    color = "black",
    linewidth = 0.25
  ) +
  facet_wrap(~ state, nrow = 1) +
  scale_fill_manual(values = tech_cols) +
  labs(
    x = NULL,
    y = "Mappability",
    title = "Mappability of detected vs undetected regions",
    fill = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )

# Save
pdf(output_pdf, width = 12, height = 4)
print(p)
dev.off()

