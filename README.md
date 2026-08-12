# DNA Methylation Dynamics of Transposable Elements During Human Embryonic Stem Cell Fate Transition

This repository contains the scripts used to generate analyses and plots for figures. The workflows compare whole-genome bisulfite sequencing (WGBS) with Oxford Nanopore Technologies (ONT) methylation data across primed, naive, and trophoblast stem cell (TSC) states, and integrate DNA methylation with transposable-element (TE) annotation, RNA expression, genomic coverage, and evolutionary conservation.

## Workflow conventions

- Scripts with numeric prefixes (`01_`, `02_`, `03_`) form a sequential workflow and must be run in ascending numeric order.
- Unnumbered scripts in the same folder are independent analyses that generate different figures; they are not sequential steps.
- Commands below assume they are run from the repository root.
- The repository contains analysis scripts only. Most large input and intermediate data files are not included.

## Repository structure

```text
Figure-1/
  Panel-B/   Methylation-rate density plots
  Panel-C/   WGBS-versus-ONT CpG concordance plots
  Panel-D/   Gene-body metagene methylation profile
  Panel-E/   Methylation-category composition plot
Figure-2/
  Panel-B/   Two independent whole-genome CpG coverage figures
  Panel-C/   Callable-genome coverage summary and plot
  Panel-D/   Mappability violin plot
  Panel-E-F/
    Bar_plot/  MAPQ10 LR-only promoter workflow
    Pie/       MAPQ10 LR-only TE-composition workflow
  Panel-G/   MAPQ10 LR-only TE-subfamily capture workflow
Figure-3/
  Panel-C/   TE-class DNA methylation violin plot
  Panel_D/   TE RNA differential-expression and classification workflow
  Panel_E/   TE conservation and methylation workflow
```

## Software requirements

### Command-line software

- Bash
- GNU `awk`, `sort`, `wc`, `gzip`/`zcat`, and coreutils
- BEDTools
- Python 3
- R
- Cairo support for R graphics where `cairo_pdf` is requested
- Optional: `nice` and `ionice` for the low-priority Figure 2 Panel C workflow

### Python packages

- `matplotlib`

```bash
python3 -m pip install matplotlib
```

### R packages

CRAN packages:

- `data.table`
- `dplyr`
- `ggplot2`
- `ggrepel`
- `ggtext`
- `patchwork`
- `cowplot`
- `scales`
- `tidyverse`

Bioconductor packages:

- `DESeq2`

```r
install.packages(c(
  "data.table", "dplyr", "ggplot2", "ggrepel", "ggtext",
  "patchwork", "cowplot", "scales", "tidyverse"
))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
```

`grid` is distributed with R and does not require a separate installation.

## Input formats used by the workflows

- **methylC:** seven tab-separated fields: chromosome, start, end, context, methylation rate, strand, and coverage.
- **Coverage bedGraph:** chromosome, start, end, and depth/value.
- **TE annotation BED:** genomic coordinates plus repeat attributes containing subfamily, class, and family information.
- **GENCODE GTF:** compressed gene annotation used to construct gene and promoter BED files.
- **SQuIRE subfamily counts:** tables containing `Subfamily:Family:Class`, `tot_counts`, or `fpkm`, depending on the analysis.
- **DNA methylation matrices/summaries:** TE subfamilies in rows with primed, naive, and TSC methylation values.
- **Dfam summaries:** TE subfamily, class, estimated age (`MYA`), and evolutionary group.

## Figure 1

All Figure 1 scripts are unnumbered and are run independently.

### Panel B — methylation-rate density

Script: `Figure-1/Panel-B/methylation_density_plot.R`

Purpose:

- Reads H9 WGBS and W3MECP2 ONT methylC files.
- Uses CpG (`CG`) records with coverage >= 10.
- The ONT inputs explicitly reference MAPQ10-primary, chrY-excluded methylC files.
- Randomly subsamples to at most 100,000 CpGs per sample and generates combined and cell-state-specific density plots.

Run:

```bash
Rscript Figure-1/Panel-B/methylation_density_plot.R
```

Outputs under `/BLUES/eric/ONT_WGBS/Figure_1/Panel_B/`:

- `Methylation_Rate_Density_WGBS_ONT_methylC_cov10.pdf`
- `Methylation_Density_Primed_methylC_cov10.pdf`
- `Methylation_Density_Naive_methylC_cov10.pdf`
- `Methylation_Density_TSC_methylC_cov10.pdf`

### Panel C — WGBS-versus-ONT CpG concordance

Script: `Figure-1/Panel-C/Density_Scatter_plot` (an R script without an `.R` suffix)

Purpose:

- Compares exact chromosome/start-matched WGBS and ONT CpGs on `chr10`.
- Retains `CG` records with coverage >= 10.
- Computes Pearson correlation and plots two-dimensional bins with width `0.03`.
- Generates separate comparisons for TSC H9, TSC AN, naive H9, naive AN, and primed samples.

Run:

```bash
Rscript Figure-1/Panel-C/Density_Scatter_plot
```

Outputs under `/BLUES/eric/ONT_WGBS/Scatter_plot/`:

- `binned_<comparison>_chr10_exact_cov10.pdf`

### Panel D — metagene methylation profile

Script: `Figure-1/Panel-D/metagene_5mc.py`

Purpose:

- Reads protein-coding genes from GENCODE v45.
- Excludes `chrM`, `chrY`, and genes shorter than 100 bp.
- Summarizes coverage-filtered methylation bedGraphs across 5 kb upstream, a scaled gene body, and 5 kb downstream.
- Uses 50 upstream, 100 gene-body, and 50 downstream bins while orienting genes 5-prime to 3-prime.

Run:

```bash
python3 Figure-1/Panel-D/metagene_5mc.py
```

Inputs:

- `/BLUES/eric/ONT_WGBS/TSS_TES_plot/gencode.v45.annotation.gtf.gz`
- Six `*.CG.cov10.bg` files under `/BLUES/eric/ONT_WGBS/TSS_TES_plot/bedgraph/`

Outputs under `/BLUES/eric/ONT_WGBS/TSS_TES_plot/`:

- `Out_mean_<sample>.txt` for each sample
- `Metagene_methylation_pm5kb_cov10.pdf`

### Panel E — methylation-category composition

Script: `Figure-1/Panel-E/Bar_plot.R`

Purpose:

- Reads six WGBS/ONT methylC files and retains CpGs with coverage >= 10.
- Classifies methylation rates into `0–0.2`, `0.2–0.8`, and `0.8–1`.
- Produces a stacked percentage bar plot, with ONT bars outlined in black.

Run:

```bash
Rscript Figure-1/Panel-E/Bar_plot.R
```

## Figure 2

### Panel B — independent whole-genome CpG figures

The two unnumbered scripts are independent and may be run in either order.

#### CpG coverage categories

Script: `Figure-2/Panel-B/CpG_Coverage.R`

Purpose:

- Counts reference CpGs in hg38 with chrY excluded.
- Counts observed CpGs at coverage >= 1 and >= 10 for six samples.
- Divides reference CpGs into uncovered, 1–9x, and >=10x categories.
- ONT inputs explicitly reference MAPQ10-primary, chrY-excluded methylC files.

Run:

```bash
Rscript Figure-2/Panel-B/CpG_Coverage.R
```

Inputs:

- `/BLUES/eric/ONT_WGBS/refs/hg38.fa`
- Six methylC files listed in the script

Intended outputs under `/BLUES/eric/ONT_WGBS/Figure_2/`:

- `CpG_coverage_6samples_counts.tsv`
- `CpG_coverage_6samples_long.tsv`
- `CpG_coverage_6samples.pdf`


#### Detected whole-genome CpGs

Script: `Figure-2/Panel-B/Detected_CpG_Whole_Genome.R`

Purpose:

- Counts position-level `CG` records detected at coverage >= 1 with chrY excluded.
- Compares ONT and WGBS across primed, naive, and TSC states.
- ONT inputs explicitly reference MAPQ10-primary, chrY-excluded methylC files.

Run:

```bash
Rscript Figure-2/Panel-B/Detected_CpG_Whole_Genome.R
```

Outputs under `/BLUES/eric/ONT_WGBS/Figure_2/`:

- `WholeGenome_noY_CpGcounts_ge1_Q10.tsv`
- `Fig2A_left_WholeGenome_noY_CpGcounts_Q10.pdf`

### Panel C — genome coverage categories

Run in numeric order:

```bash
bash Figure-2/Panel-C/01_make_summary.sh
Rscript Figure-2/Panel-C/02_plot_genome_coverage_detected_bases.R
```

#### Step 01: build the coverage summary

`01_make_summary.sh`:

- Builds chrY-excluded chromosome sizes and genome intervals.
- Converts the hg38 gap table into merged N-only intervals.
- Defines callable sequence as the chrY-excluded genome minus N-only intervals.
- Counts non-N bases detected at coverage >= 1 from WGBS and MAPQ10-primary ONT coverage bedGraphs.

Principal inputs:

- `/BLUES/eric/WGBS/hg38.chrom.sizes`
- `/BLUES/eric/ONT_WGBS/Figure_2/Panel_A/hg38_gap.txt.gz`
- Three WGBS and three `ONT_*_MAPQ10_primary.bedGraph.gz` files

Outputs:

- `/BLUES/eric/ONT_WGBS/Figure_2/hg38.nochrY.chrom.sizes`
- `/BLUES/eric/ONT_WGBS/Figure_2/N_only_noY.raw.bed`
- `/BLUES/eric/ONT_WGBS/Figure_2/N_only_noY.merged.bed`
- `/BLUES/eric/ONT_WGBS/Figure_2/genome_noY.bed`
- `/BLUES/eric/ONT_WGBS/Figure_2/callable_noY.bed`
- `/BLUES/eric/ONT_WGBS/Figure_2/Panel_A/Fig2A_right_genome_fraction_3cat_MAPQ10_primary.tsv`

#### Step 02: plot the summary

`02_plot_genome_coverage_detected_bases.R` reads the summary TSV from Step 01 and writes:

- `/BLUES/eric/ONT_WGBS/Figure_2/Panel_A/Fig2A_genome_fraction_MAPQ.pdf`

Although these scripts are stored under `Figure-2/Panel-C/` in the repository, their hard-coded working outputs are under `Figure_2/Panel_A/`.

### Panel D — mappability violin plot

Script: `Figure-2/Panel-D/mappability_violin_plot.R`

Purpose:

- Reads a precomputed table of mappability values for detected and undetected intervals.
- Draws WGBS-versus-ONT violin and internal box plots faceted by cell state.

Run:

```bash
Rscript Figure-2/Panel-D/mappability_violin_plot.R
```

Input:

- `/BLUES/eric/ONT_WGBS/Figure_2/Panel_B/B2/B2_mappability_violin.tsv`

Output:

- `/BLUES/eric/ONT_WGBS/Figure_2/Panel_B/B2/Fig2B_mappability_violin.pdf`

The script does not generate the input TSV; that preprocessing step is not present in this repository.

### Panels E–F — LR-only promoter bar plot

Run in numeric order:

```bash
bash Figure-2/Panel-E-F/Bar_plot/01_prepare_gencode_gene_promoter_inputs.sh
bash Figure-2/Panel-E-F/Bar_plot/02_count_MAPQ10_LRonly_promoter_1kb_100pct.sh
Rscript Figure-2/Panel-E-F/Bar_plot/03_plot_MAPQ10_LRonly_promoter_1kb_100pct.R
```

#### Step 01: prepare annotation and promoter inputs

- Parses GENCODE v47 gene records into BED6 format.
- Retains primary chromosomes and creates +/-1 kb, +/-2 kb, and +/-5 kb promoter BED files.
- Also merges the legacy `LR_only_windows_100bp.bed`; this legacy LR-only intermediate is not consumed by Steps 02–03.

Main input:

- `/BLUES/eric/refs/D/gencode.v47.primary_assembly.annotation.gtf.gz`

Key outputs under `/BLUES/eric/ONT_WGBS/Figure_2/Panel_E/`:

- `inputs/hg38.primary.chrom.sizes`
- `intermediate/gencode_v47_genes.primary.sorted.bed`
- `intermediate/gencode_v47_promoters_{1000,2000,5000}bp.primary.sorted.bed`

#### Step 02: identify promoters fully covered by MAPQ10 LR-only regions

- Uses `Primed_LR_only_MAPQ10_windows_100bp.bed`.
- Merges MAPQ10 LR-only windows.
- Creates protein-coding-specific gene and +/-1 kb promoter BEDs.
- Uses `bedtools coverage` to retain promoters for which covered bases equal promoter length.

Principal input:

- `/BLUES/eric/ONT_WGBS/Figure_2/LR_only_MAPQ/Primed_LR_only_MAPQ10_windows_100bp.bed`

Key outputs under `Figure_2/Panel_E/outputs/`:

- `promoters_1000bp_LRonly_MAPQ10_coverage.tsv`
- `promoters_1000bp_100pct_LRonly_MAPQ10_regions.bed`
- `promoters_1000bp_100pct_LRonly_MAPQ10_regions.summary.tsv`
- `protein_coding_promoters_1000bp_LRonly_MAPQ10_coverage.tsv`
- `protein_coding_promoters_1000bp_100pct_LRonly_MAPQ10_regions.bed`
- `protein_coding_promoters_1000bp_100pct_LRonly_MAPQ10_regions.summary.tsv`

#### Step 03: plot promoter counts

- Counts unique all-gene and protein-coding gene IDs from the Step 02 summaries.
- Writes:
  - `LRonly_MAPQ10_overlap_100pct_promoter1kb_barplot.pdf`
  - `LRonly_MAPQ10_overlap_100pct_promoter1kb_barplot.png`

### Panels E–F — LR-only TE-composition pie plots

This is a separate sequential workflow from the bar-plot folder:

```bash
bash Figure-2/Panel-E-F/Pie/01_make_summary_pie_LR_only.sh
Rscript Figure-2/Panel-E-F/Pie/02_plot_pie.R \
  /BLUES/eric/ONT_WGBS/Figure_2/Panel_C
```

#### Step 01: prepare LR-only and TE intervals

- Sorts the primed MAPQ10 LR-only 100 bp windows.
- Parses subfamily, class, and family from the chrY-excluded TE annotation.
- Identifies whole LR-only windows that have no TE overlap and merges adjacent no-TE windows.

Inputs:

- `/BLUES/eric/ONT_WGBS/Figure_2/LR_only_MAPQ/Primed_LR_only_MAPQ10_windows_100bp.bed`
- `/BLUES/eric/ONT/hg38_TE_noY.bed.gz`

Key outputs under `/BLUES/eric/ONT_WGBS/Figure_2/Panel_C/`:

- `intermediate/LR_only_windows_100bp.sorted.bed`
- `intermediate/TE.parsed.bed`
- `intermediate/TE.parsed.sorted.bed`
- `outputs/LR_only_windows_100bp.noTE.bed`
- `outputs/LR_only_regions.noTE.merged.bed`

#### Step 02: calculate overlap composition and plot

- Calculates unique LR-only and TE-overlapping base-pair lengths.
- Produces TE versus non-TE composition and top-10 class, family, and subfamily pie charts.
- Writes `outputs/Fig2C_TE_pies_top10_legend.pdf`.

### Panel G — top TE-subfamily fraction captured

Run in numeric order:

```bash
bash Figure-2/Panel-G/01_Make_inputs.sh
Rscript Figure-2/Panel-G/02_plot_panel_G.R \
  /BLUES/eric/ONT_WGBS/Figure_2/Panel_D
```

#### Step 01: calculate capture fractions

- Merges primed MAPQ10 LR-only windows.
- Parses and sorts the chrY-excluded TE annotation.
- Defines the ten TE subfamilies listed in the script.
- Counts total genome-wide TE records and TE records overlapping LR-only regions.
- Calculates the percentage of genome-wide copies captured in LR-only regions.

Outputs under `/BLUES/eric/ONT_WGBS/Figure_2/Panel_D/`:

- `inputs/top10_subfamilies_from_pie.txt`
- `intermediate/top10_total_genome_copies.tsv`
- `intermediate/TE_overlapping_LRonly.bed`
- `intermediate/top10_LRonly_copies.tsv`
- `outputs/Fig2D_top10_fractionCaptured.tsv`

#### Step 02: plot capture fractions

- Reads the fraction table and the ordered top-10 list from Step 01.
- Plots percent captured against TE subfamily; point size represents copies in LR-only regions and color represents TE family.
- Writes `outputs/Fig2D_dotplot_top10_fractionCaptured.pdf`.

## Figure 3

### Panel C — TE-class methylation distributions

Script: `Figure-3/Panel-C/Plot_Violin_TE_class_DNA_methylation.R`

Purpose:

- Joins a full TE-subfamily DNA methylation matrix to subfamily/class/family annotations parsed from the TE BED.
- Removes selected RNA-derived labels, requires complete primed/naive/TSC values, and plots subfamily methylation distributions by RepeatMasker class.

Run:

```bash
Rscript Figure-3/Panel-C/Plot_Violin_TE_class_DNA_methylation.R
```

Inputs:

- `/BLUES/eric/ONT_WGBS/Heatmap/DNA_methylation_matrix_FULL_no_rRNA_tRNA.tsv`
- `/BLUES/eric/ONT/hg38_TE_noY.bed.gz`

Outputs under `/BLUES/eric/ONT_WGBS/Figure_3/Violin_plot/`:

- `TE_repeatClass_DNA_methylation_violin_allClasses.pdf`
- `TE_repeatClass_DNA_methylation_violin_allClasses.png`
- `TE_repeatClass_DNA_methylation_violin_class_summary.tsv`
- `TE_repeatClass_DNA_methylation_violin_unannotated_subfamilies.tsv`, when applicable
- `ambiguous_subfamily_repeatClass_repeatFamily.tsv`, when applicable

The methylation-matrix filename does not establish MAPQ10-primary filtering; verify this input before a strict Q10 reproduction.

### Panel D — differential TE expression and grid classification

Run in numeric order:

```bash
Rscript Figure-3/Panel_D/01_DEseq2_TE_subfam.R
Rscript Figure-3/Panel_D/02_plot_grid_based_classification_deseq
```

#### Step 01: DESeq2 analysis of SQuIRE subfamily counts

- Reads two SQuIRE replicates for each of the primed, naive, and TSC conditions.
- Constructs a rounded integer count matrix from `tot_counts`.
- Retains subfamilies with at least 10 total counts.
- Calculates Naive-versus-Primed and TSC-versus-Naive DESeq2 contrasts.

Outputs under `/BLUES/eric/ONT_WGBS/Figure_3/MA_Plot/deseq/`:

- `SQuIRE_tot_counts_matrix.tsv`
- `SQuIRE_sample_info.tsv`
- `RNAseq_log2FC_DESeq2.tsv`

#### Step 02: integrate RNA and MAPQ10-primary DNA methylation

- Reads `DNA_methylation_matrix_FULL_MAPQ10_primary_noY_no_rRNA_tRNA.tsv` and the Step 01 DESeq2 result.
- Calculates Naive-minus-Primed and TSC-minus-Naive methylation changes.
- Applies the committed RNA (`2`) and DNA methylation (`0.1`) classification cutoffs.
- Writes:
  - `Grid_based_classification_DESeq2_log2FC.tsv`
  - `Grid_based_classification_DESeq2_log2FC.pdf`

### Panel E — TE age versus methylation difference

Run in numeric order:

```bash
Rscript Figure-3/Panel_E/01_make_summary_input.R
Rscript Figure-3/Panel_E/02_Plot_conservation.R
```

#### Step 01: combine conservation, expression, and methylation summaries

- Combines protein-coding and lincRNA Dfam TE-age summaries.
- Joins mean SQuIRE FPKM and TE-subfamily DNA methylation summaries by subfamily.
- Writes `/BLUES/eric/ONT_WGBS/Figure_3/Jitter_plot_conservation/TE_conservation_methylation_expression.tsv`.

Inputs:

- `Dfam_TE_inf_PC_summary.tsv`
- `Dfam_TE_inf_lincRNA_summary.tsv`
- `/BLUES/eric/ONT_WGBS/Heatmap/W3_TEsubfamily_RNAseq_summary_meanFPKM_no_rRNA_tRNA.tsv`
- `/BLUES/eric/ONT_WGBS/Heatmap/W3_TEsubfamily_DNA_methylation_summary.tsv`

The methylation-summary filename does not establish MAPQ10-primary filtering; verify this input before a strict Q10 reproduction.

#### Step 02: plot TE age against methylation differences

- Calculates Naive-minus-Primed and TSC-minus-Naive methylation differences.
- Plots TE insertion time against methylation difference, colored by evolutionary group and shaped by TE class.
- Labels Homo-group TEs and methylation-difference outliers.
- Writes `TE_age_vs_methylation_difference.pdf`.

## Absolute paths and portability

All scripts were written for a specific HPC filesystem and contain hard-coded paths beginning with:

```text
/BLUES/eric/
```

Other users must change these paths to match their local or cluster filesystem before running the workflows. This includes input paths, output directories, reference files, and paths embedded inside inline Python blocks. Simply cloning the repository will not provide the large external datasets.

Some repository folder names also differ from the hard-coded output locations—for example, scripts stored under `Figure-2/Panel-C/` write into `Figure_2/Panel_A/`. Follow the paths in each script or update them consistently for the target environment.

