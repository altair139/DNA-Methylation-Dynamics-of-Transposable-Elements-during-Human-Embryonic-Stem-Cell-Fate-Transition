suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

in_tsv <- "/BLUES/eric/ONT_WGBS/Figure_2/Panel_A/Fig2A_right_genome_fraction_3cat_MAPQ10_primary.tsv"
out_pdf <- "/BLUES/eric/ONT_WGBS/Figure_2/Panel_A/Fig2A_genome_fraction_MAPQ.pdf"

dt <- fread(in_tsv)
dt <- dt[, .(state, technology, category, bases=as.numeric(bases), pct=as.numeric(percent_genome))]


# Sample label
dt[, sample := factor(paste(state, technology), levels = rev(c(
  "Primed WGBS","Primed ONT",
  "Naive WGBS","Naive ONT",
  "TSC WGBS", "TSC ONT"
)))]

# STACK ORDER
stack_levels <- c(
  "Detected (≥1×, non-N)",
  "Undetected (<1×, non-N)",
  "Non-detectable (N-only)"
)
dt[, category := factor(category, levels = stack_levels)]

# Labels: 
dt[, lab := sprintf("%.1f%%", pct)]
dt[pct < 2, lab := ""]   

# Colors mapped to category names
cols <- c(
  "Detected (≥1×, non-N)"      = "#F87171",
  "Undetected (<1×, non-N)"    = "#E5E7EB",
  "Non-detectable (N-only)"    = "#9CA3AF"
)

p <- ggplot(dt, aes(x=pct, y=sample, fill=category)) +
  geom_col(width=0.72, color="black", linewidth=0.25) +
  geom_text(aes(label=lab),
            position=position_stack(vjust=0.5),
            size=3.6) +
  scale_fill_manual(values=cols, drop=FALSE) +
  scale_x_continuous(limits=c(0,100),
                     breaks=c(0,25,50,75,100),
                     expand=expansion(mult=c(0,0.01))) +
  labs(
    title="Genome coverage categories (chrY excluded)",
    subtitle="Detected = ≥1× in callable (non-N) bases; Undetected = <1× in callable bases; Non-detectable = N-only bases",
    x="Genome fraction (%)", y=NULL, fill=NULL
  ) +
  theme_minimal(base_size=13) +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position="top",
    axis.text.y = element_text(size=11, color="black"),
    axis.text.x = element_text(size=11, color="black"),
    axis.title.x = element_text(size=12, color="black"),
    plot.title = element_text(size=16, face="bold"),
    plot.subtitle = element_text(size=11)
  )

ggsave(out_pdf, p, width=12, height=4.5, device=cairo_pdf, bg="white")
cat("[ok] wrote:", out_pdf, "\n")
