library(tidyverse)
library(ggplot2)

MIN_COV <- 10

files <- c(
  # H9 (WGBS methylC)
  "/BLUES/eric/WGBS/CpG/H9-primed_onC.methylC.gz",
  "/BLUES/eric/WGBS/CpG/H9-Naive_onC.methylC.gz",
  "/BLUES/eric/WGBS/CpG/H9-TSC_onC.methylC.gz",

  # W3MECP2 (ONT methylC)
  "/BLUES/eric/ONT/Primed/Primed_W3MECP2_aln.methylC.gz",
  "/BLUES/eric/ONT/Naive/091724_Naive_W3MECPC2_aln.methylC.gz",
  "/BLUES/eric/ONT/TSC/TSC_GFP_W3MECP2_merged_aln.methylC.gz"
)

labels <- c(
  "H9 primed",
  "H9 naive",
  "H9 TSC",
  "W3MECP2 primed",
  "W3MECP2 naive",
  "W3MECP2 TSC"
)

# Order
desired_order <- c(
  "H9 primed",
  "W3MECP2 primed",
  "H9 naive",
  "W3MECP2 naive",
  "H9 TSC",
  "W3MECP2 TSC"
)

read_methylC <- function(file, label, min_cov = 10) {
  read_tsv(
    file,
    comment = "#",
    col_names = c("chr","start","end","context","meth_rate","strand","coverage"),
    show_col_types = FALSE
  ) %>%
    filter(context == "CG") %>%
    filter(!is.na(meth_rate), !is.na(coverage)) %>%
    filter(coverage >= min_cov) %>%
    mutate(
      category = case_when(
        meth_rate < 0.2 ~ "0-0.2",
        meth_rate < 0.8 ~ "0.2-0.8",
        meth_rate <= 1  ~ "0.8-1",
        TRUE            ~ NA_character_
      ),
      sample = label
    ) %>%
    filter(!is.na(category))
}

all_data <- map2_dfr(files, labels, ~read_methylC(.x, .y, min_cov = MIN_COV))

plot_df <- all_data %>%
  count(sample, category) %>%
  group_by(sample) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup() %>%
  mutate(
    sample = factor(sample, levels = desired_order),
    is_ONT = grepl("W3MECP2", sample),

    # KEY FIX (color direction):
    # Use this order so AFTER coord_flip(), bars read left->right:
    # 0-0.2 (blue) -> 0.2-0.8 (green) -> 0.8-1 (orange)
    category = factor(category, levels = c("0.8-1", "0.2-0.8", "0-0.2"))
  )

col_map <- c(
  "0-0.2"   = scales::alpha("#8ecae6", 0.85),
  "0.2-0.8" = scales::alpha("#a6e3a1", 0.85),
  "0.8-1"   = scales::alpha("#f4a261", 0.85)
)

# Outside label only for naive + 0.8-1
outside_naive <- plot_df %>%
  filter(grepl("naive", as.character(sample), ignore.case = TRUE),
         category == "0.8-1") %>%
  mutate(
    # 0.8-1 will appear on the RIGHT in the flipped plot with the category order above
    seg_mid  = 100 - (percent / 2),
    line_end = 104,
    text_pos = 106
  )

inside_labels <- plot_df %>%
  anti_join(outside_naive, by = c("sample", "category"))

p <- ggplot(plot_df, aes(x = sample, y = percent, fill = category)) +

  geom_col(
    data = plot_df %>% filter(!is_ONT),
    width = 0.8
  ) +
  geom_col(
    data = plot_df %>% filter(is_ONT),
    width = 0.8,
    color = "black",
    linewidth = 1.2
  ) +

  # Labels 
  geom_text(
    data = inside_labels,
    aes(label = sprintf("%.1f%%", percent)),
    position = position_stack(vjust = 0.5),
    size = 3.6
  ) +

  # Leader line + outside label
  geom_segment(
    data = outside_naive,
    aes(x = sample, xend = sample, y = seg_mid, yend = line_end),
    linewidth = 0.6
  ) +
  geom_text(
    data = outside_naive,
    aes(x = sample, y = text_pos, label = sprintf("%.1f%%", percent)),
    hjust = 0,
    size = 3.6
  ) +

  coord_flip(clip = "off") +

  # Forcing the order
  scale_x_discrete(limits = rev(desired_order)) +

  scale_y_continuous(
    limits = c(0, 112),
    breaks = c(0, 25, 50, 75, 100),
    expand = c(0, 0)
  ) +

  scale_fill_manual(
    values = col_map,
    name   = "Methylation range",
    breaks = c("0.8-1", "0.2-0.8", "0-0.2")
  ) +

  labs(
    title = "Percent methylation composition",
    x = NULL,
    y = "Percentage"
  ) +

  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    plot.margin = margin(t = 10, r = 30, b = 10, l = 10)
  )

print(p)

ggsave(
  "Combined_methylC_barplot_cov10.pdf",
  p,
  width = 9,
  height = 6,
  device = cairo_pdf,
  bg = "white"
)
