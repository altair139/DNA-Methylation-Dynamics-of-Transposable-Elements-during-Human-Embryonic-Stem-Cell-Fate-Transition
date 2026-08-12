#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

setDTthreads(8)

# =====================================================
# INPUT / OUTPUT
# =====================================================

infile <- "/BLUES/eric/ONT_WGBS/Heatmap/DNA_methylation_matrix_FULL_no_rRNA_tRNA.tsv"
te_bed <- "/BLUES/eric/ONT/hg38_TE_noY.bed.gz"

outdir <- "/BLUES/eric/ONT_WGBS/Figure_3/Violin_plot"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

out_pdf <- file.path(
  outdir,
  "TE_repeatClass_DNA_methylation_violin_allClasses.pdf"
)

out_png <- file.path(
  outdir,
  "TE_repeatClass_DNA_methylation_violin_allClasses.png"
)

out_class_summary <- file.path(
  outdir,
  "TE_repeatClass_DNA_methylation_violin_class_summary.tsv"
)

out_unannotated <- file.path(
  outdir,
  "TE_repeatClass_DNA_methylation_violin_unannotated_subfamilies.tsv"
)

# =====================================================
# LOAD TE ANNOTATION FROM BED
# =====================================================

load_te_annotation <- function(te_bed) {

  te <- fread(cmd = paste("zcat", shQuote(te_bed)), header = FALSE)

  # Expected:
  # chr start end strand attributes
  # repeat_id "SVA_E,chr9,15838676,15838746";repeat_class "Retroposon";repeat_family "SVA";

  te[, repeat_id := sub('.*repeat_id "([^"]+)".*', '\\1', V5)]
  te[, repeat_class := sub('.*repeat_class "([^"]+)".*', '\\1', V5)]
  te[, repeat_family := sub('.*repeat_family "([^"]+)".*', '\\1', V5)]

  te[, subfamily := tstrsplit(repeat_id, ",", fixed = TRUE, keep = 1L)]

  anno <- unique(te[, .(subfamily, repeat_class, repeat_family)])

  # If one subfamily maps to multiple classes/families, keep one sorted record
  ambig <- anno[, .(
    n_class = uniqueN(repeat_class),
    n_family = uniqueN(repeat_family)
  ), by = subfamily][n_class > 1 | n_family > 1]

  if (nrow(ambig) > 0) {
    warning(sprintf(
      "%d subfamilies map to >1 repeat_class/repeat_family. Keeping first sorted annotation.",
      nrow(ambig)
    ))

    fwrite(
      ambig,
      file.path(outdir, "ambiguous_subfamily_repeatClass_repeatFamily.tsv"),
      sep = "\t"
    )
  }

  setorder(anno, subfamily, repeat_class, repeat_family)
  anno <- anno[, .SD[1], by = subfamily]

  return(anno)
}

anno <- load_te_annotation(te_bed)

# =====================================================
# LOAD METHYLATION MATRIX
# =====================================================

dt <- fread(infile)

dt <- dt[!grepl("^5S$|^7SK$|^7SLRNA$", subfamily)]

# Keep rows with complete methylation values
dt <- dt[complete.cases(dt[, .(Primed, Naive, TSC)])]

# Join true RepeatMasker class/family
dt <- merge(dt, anno, by = "subfamily", all.x = TRUE)

# Use true repeat_class directly
dt[, plot_class := repeat_class]
dt[is.na(plot_class) | plot_class == "", plot_class := "Unannotated"]

dt[, repeat_family := fifelse(
  is.na(repeat_family) | repeat_family == "",
  "Unannotated",
  repeat_family
)]

# =====================================================
# CLASS ORDER
# =====================================================

preferred_class_order <- c(
  "LINE",
  "SINE",
  "LTR",
  "Retroposon",
  "DNA",
  "RC",
  "Satellite",
  "Simple_repeat",
  "Low_complexity",
  "RNA",
  "Unknown",
  "DNA?",
  "LTR?",
  "RC?",
  "scRNA",
  "snRNA",
  "srpRNA",
  "tRNA",
  "rRNA",
  "Unannotated"
)

observed_classes <- unique(dt$plot_class)

class_order <- c(
  preferred_class_order[preferred_class_order %in% observed_classes],
  sort(setdiff(observed_classes, preferred_class_order))
)

dt[, plot_class := factor(plot_class, levels = class_order)]

# =====================================================
# CLASS SUMMARY 
# =====================================================

class_summary <- dt[, .N, by = plot_class][order(plot_class)]
fwrite(class_summary, out_class_summary, sep = "\t")

unannotated <- dt[plot_class == "Unannotated", .(
  subfamily,
  Primed,
  Naive,
  TSC
)]

if (nrow(unannotated) > 0) {
  fwrite(unannotated, out_unannotated, sep = "\t")
  cat("\nUnannotated subfamilies saved to:\n")
  cat(out_unannotated, "\n")
}

# =====================================================
# LONG FORMAT
# =====================================================

df <- melt(
  dt,
  id.vars = c("subfamily", "repeat_class", "repeat_family", "plot_class"),
  measure.vars = c("Primed", "Naive", "TSC"),
  variable.name = "celltype",
  value.name = "beta"
)

df[, celltype := factor(
  celltype,
  levels = c("Primed", "Naive", "TSC")
)]

df[, plot_class := factor(plot_class, levels = class_order)]

# =====================================================
# PLOT
# =====================================================

p <- ggplot(
  df,
  aes(
    x = celltype,
    y = beta,
    fill = celltype
  )
) +

  geom_violin(
    scale = "width",
    trim = TRUE,
    color = "black",
    linewidth = 0.2
  ) +

  geom_boxplot(
    width = 0.12,
    outlier.size = 0.2,
    fill = "white",
    linewidth = 0.25
  ) +

  facet_wrap(
    ~ plot_class,
    nrow = 3,
    scales = "fixed"
  ) +

  scale_fill_manual(
    values = c(
      "Primed" = "#d7301f",
      "Naive"  = "#fee08b",
      "TSC"    = "#4575b4"
    )
  ) +

  labs(
    title = "DNA methylation distribution across RepeatMasker classes",
    subtitle = "Violin plot of TE subfamily methylation beta values grouped by true repeat_class",
    x = NULL,
    y = "DNA methylation (beta)"
  ) +

  coord_cartesian(
    ylim = c(0, 1)
  ) +

  theme_bw(base_size = 12) +

  theme(
    strip.background = element_rect(
      fill = "grey95",
      color = "black"
    ),

    strip.text = element_text(
      face = "bold",
      size = 11
    ),

    axis.text.x = element_text(
      size = 9,
      face = "bold",
      angle = 45,
      hjust = 1
    ),

    axis.text.y = element_text(
      size = 10
    ),

    plot.title = element_text(
      size = 15,
      face = "bold",
      hjust = 0.5
    ),

    plot.subtitle = element_text(
      size = 11,
      hjust = 0.5
    ),

    legend.position = "none"
  )

# =====================================================
# SAVE
# =====================================================

ggsave(
  out_pdf,
  p,
  width = 13,
  height = 9,
  useDingbats = FALSE
)

ggsave(
  out_png,
  p,
  width = 13,
  height = 9,
  dpi = 300
)

cat("\nSaved:\n")
cat(out_pdf, "\n")
cat(out_png, "\n")
cat(out_class_summary, "\n\n")
