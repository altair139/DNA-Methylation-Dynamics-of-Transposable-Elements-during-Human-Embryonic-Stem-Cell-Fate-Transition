#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(dplyr)
})


# Settings

setDTthreads(threads = 8)

output_dir <- "/BLUES/eric/ONT_WGBS/Figure_1/Panel_B"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

DENSITY_ADJUST <- 1.8  
EPS <- 1e-4
COV_MIN <- 10
MAX_N <- 100000


# Sample definitions 

samples <- list(
  # WGBS methylC
  "H9-Naive"  = list(file="/BLUES/eric/WGBS/CpG/H9-Naive_onC.methylC.gz",   label="H9 Naive",   type="Naive",  platform="WGBS"),
  "H9-TSC"    = list(file="/BLUES/eric/WGBS/CpG/H9-TSC_onC.methylC.gz",     label="H9 TSC",     type="TSC",    platform="WGBS"),
  "H9-Primed" = list(file="/BLUES/eric/WGBS/CpG/H9-primed_onC.methylC.gz",  label="H9 Primed",  type="Primed", platform="WGBS"),

  # ONT methylC
  "W3-Naive"  = list(file="/BLUES/eric/ONT/MAPQ_bam/ONT_Naive.MAPQ10.primary.noY.onC.methylC.gz",
                     label="W3MECP2 Naive", type="Naive",  platform="ONT"),
  "W3-TSC"    = list(file="/BLUES/eric/ONT/MAPQ_bam/ONT_TSC.MAPQ10.primary.noY.onC.methylC.gz",
                     label="W3MECP2 TSC",   type="TSC",    platform="ONT"),
  "W3-Primed" = list(file="/BLUES/eric/ONT/MAPQ_bam/ONT_Primed.MAPQ10.primary.noY.onC.methylC.gz",
                     label="W3MECP2 Primed", type="Primed", platform="ONT")
)


# Reader for methylC.gz 

read_methylC <- function(file, label, type, platform) {
  if (!file.exists(file)) stop(paste("File does not exist:", file))

  df <- fread(file, header = FALSE, showProgress = FALSE)
  if (ncol(df) < 7) stop(paste("Unexpected methylC columns:", file, "ncol=", ncol(df)))

  # V1 chr, V2 start, V3 end, V4 context, V5 meth_rate, V6 strand, V7 coverage
  dt <- df[, .(
    chr = as.character(V1),
    start = as.integer(V2),
    end = as.integer(V3),
    context = as.character(V4),
    meth_rate = as.numeric(V5),
    strand = as.character(V6),
    cov = as.numeric(V7)
  )]

  # Filter to CpG only + coverage threshold
  dt <- dt[context == "CG" & cov >= COV_MIN]

  # Subsample for speed/plot readability
  if (nrow(dt) > MAX_N) dt <- dt[sample(.N, MAX_N)]

  # Clip away from exact boundaries to reduce KDE spikes at 0 or 1
  dt[, meth_rate := pmin(pmax(meth_rate, EPS), 1 - EPS)]

  dt[, sample := label]
  dt[, type := type]
  dt[, platform := platform]

  return(dt[, .(meth_rate, sample, type, platform)])
}


# Load all data

all_data <- rbindlist(lapply(samples, function(s) {
  read_methylC(s$file, s$label, s$type, s$platform)
}))

# Theme
theme_clean <- theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank())


# Combined plot 

combined_plot <- ggplot(all_data, aes(x = meth_rate, color = sample)) +
  geom_density(linewidth = 1, adjust = DENSITY_ADJUST) +
  labs(
    title = "Methylation Rate Density - WGBS and ONT (methylC, cov>=10)",
    x = "Methylation Rate",
    y = "Density"
  ) +
  coord_cartesian(ylim = c(0, 8)) +
  theme_clean

pdf(file.path(output_dir, "Methylation_Rate_Density_WGBS_ONT_methylC_cov10.pdf"), width = 8, height = 6)
print(combined_plot)
dev.off()


# Per-cell-type plots (with fill)

for (cell_type in unique(all_data$type)) {
  sub <- filter(all_data, type == cell_type)

  p <- ggplot(sub, aes(x = meth_rate, fill = sample)) +
    geom_density(alpha = 0.4, color = "black", adjust = DENSITY_ADJUST) +
    labs(
      title = paste("Methylation Density -", cell_type, "(methylC, cov>=10)"),
      x = "Methylation Rate",
      y = "Density"
    ) +
    coord_cartesian(ylim = c(0, 8)) +
    theme_clean

  pdf(file.path(output_dir, paste0("Methylation_Density_", cell_type, "_methylC_cov10.pdf")),
      width = 8, height = 6)
  print(p)
  dev.off()
}

cat("Done. Outputs written to: ", output_dir, "\n", sep = "")
