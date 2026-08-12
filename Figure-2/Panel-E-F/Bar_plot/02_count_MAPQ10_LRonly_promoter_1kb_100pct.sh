#!/usr/bin/env bash
set -euo pipefail

OUTDIR="/BLUES/eric/ONT_WGBS/Figure_2/Panel_E"

# MAPQ-filtered LR-only 100bp windows
LR_MAPQ_WIN="/BLUES/eric/ONT_WGBS/Figure_2/LR_only_MAPQ/Primed_LR_only_MAPQ10_windows_100bp.bed"

# Existing Panel E files
GENOME="${OUTDIR}/inputs/hg38.primary.chrom.sizes"
ALL_GENES="${OUTDIR}/intermediate/gencode_v47_genes.primary.sorted.bed"

# All-gene promoter file
ALL_PROM1KB="${OUTDIR}/intermediate/gencode_v47_promoters_1000bp.primary.sorted.bed"

# Protein-coding-specific files
PCGENES="${OUTDIR}/intermediate/gencode_v47_protein_coding_genes.primary.sorted.bed"
PC_PROM1KB="${OUTDIR}/intermediate/gencode_v47_protein_coding_promoters_1000bp.primary.sorted.bed"

# MAPQ-specific LR-only merged file
LR_MAPQ_MERGED="${OUTDIR}/intermediate/LR_only_MAPQ10_merged.primary.sorted.bed"

# All-gene outputs
ALL_COV="${OUTDIR}/outputs/promoters_1000bp_LRonly_MAPQ10_coverage.tsv"
ALL_OUTBED="${OUTDIR}/outputs/promoters_1000bp_100pct_LRonly_MAPQ10_regions.bed"
ALL_SUMMARY="${OUTDIR}/outputs/promoters_1000bp_100pct_LRonly_MAPQ10_regions.summary.tsv"

# Protein-coding outputs
PC_COV="${OUTDIR}/outputs/protein_coding_promoters_1000bp_LRonly_MAPQ10_coverage.tsv"
PC_OUTBED="${OUTDIR}/outputs/protein_coding_promoters_1000bp_100pct_LRonly_MAPQ10_regions.bed"
PC_SUMMARY="${OUTDIR}/outputs/protein_coding_promoters_1000bp_100pct_LRonly_MAPQ10_regions.summary.tsv"

mkdir -p "${OUTDIR}/intermediate" "${OUTDIR}/outputs"

echo "LR_MAPQ_WIN = $LR_MAPQ_WIN"
echo "ALL_GENES   = $ALL_GENES"
echo "GENOME      = $GENOME"

[[ -s "$LR_MAPQ_WIN" ]] || { echo "ERROR: missing $LR_MAPQ_WIN" >&2; exit 1; }
[[ -s "$ALL_GENES" ]] || { echo "ERROR: missing $ALL_GENES. Run 01_prepare_gencode_gene_promoter_inputs.sh first." >&2; exit 1; }
[[ -s "$GENOME" ]] || { echo "ERROR: missing $GENOME" >&2; exit 1; }

# ------------------------------------------------------------
# 1. Make merged MAPQ10 LR-only regions
# ------------------------------------------------------------
echo "Making merged primary MAPQ10 LR-only regions..."

awk 'BEGIN{OFS="\t"} $1 ~ /^chr([0-9]+|X|Y)$/ {print $1,$2,$3}' "$LR_MAPQ_WIN" \
  | bedtools sort -g "$GENOME" \
  | bedtools merge -i - \
  | bedtools sort -g "$GENOME" \
  > "$LR_MAPQ_MERGED"

echo "MAPQ10 LR-only merged regions:"
wc -l "$LR_MAPQ_MERGED"

# ------------------------------------------------------------
# 2. Make all-gene promoter ±1kb if missing
# ------------------------------------------------------------
if [[ ! -s "$ALL_PROM1KB" ]]; then
  echo "Creating all-gene promoter ±1kb BED..."

  awk 'BEGIN{FS=OFS="\t"}
  {
    chr=$1
    start=$2
    end=$3
    name=$4
    strand=$6

    if(strand=="+") {
      tss=start
    } else if(strand=="-") {
      tss=end
    } else {
      next
    }

    ps=tss-1000
    pe=tss+1000
    if(ps<0) ps=0

    print chr,ps,pe,name,0,strand
  }' "$ALL_GENES" \
    | bedtools sort -g "$GENOME" \
    > "$ALL_PROM1KB"
fi

# ------------------------------------------------------------
# 3. Make protein-coding genes and protein-coding promoter ±1kb
# ------------------------------------------------------------
echo "Creating protein-coding gene BED..."

awk 'BEGIN{FS=OFS="\t"}
{
  split($4,a,"|")
  if(a[3]=="protein_coding") print
}' "$ALL_GENES" \
  | bedtools sort -g "$GENOME" \
  > "$PCGENES"

echo "Protein-coding genes:"
wc -l "$PCGENES"

echo "Creating protein-coding promoter ±1kb BED..."

awk 'BEGIN{FS=OFS="\t"}
{
  chr=$1
  start=$2
  end=$3
  name=$4
  strand=$6

  if(strand=="+") {
    tss=start
  } else if(strand=="-") {
    tss=end
  } else {
    next
  }

  ps=tss-1000
  pe=tss+1000
  if(ps<0) ps=0

  print chr,ps,pe,name,0,strand
}' "$PCGENES" \
  | bedtools sort -g "$GENOME" \
  > "$PC_PROM1KB"

echo "[info] Protein-coding promoters ±1kb:"
wc -l "$PC_PROM1KB"

# ------------------------------------------------------------
# 4. Function: keep promoters 100% covered by MAPQ10 LR-only
# ------------------------------------------------------------
run_100pct_promoter_overlap () {
  local label="$1"
  local prom_bed="$2"
  local cov_file="$3"
  local out_bed="$4"
  local out_summary="$5"

  echo
  echo "[info] Processing $label"

  bedtools coverage \
    -a "$prom_bed" \
    -b "$LR_MAPQ_MERGED" \
    > "$cov_file"

  # bedtools coverage output for BED6:
  # columns 1-6 = promoter
  # column 7 = overlap_count
  # column 8 = covered_bases
  # column 9 = promoter_length
  # column 10 = fraction_covered
  #
  # Require covered_bases == promoter_length
  awk 'BEGIN{FS=OFS="\t"} $(NF-2) == $(NF-1) {print $1,$2,$3,$4,$5,$6}' "$cov_file" \
    > "$out_bed"

  awk 'BEGIN{FS=OFS="\t"}
  {
    split($4,a,"|")
    gene_id=a[1]
    gene_name=a[2]
    gene_type=a[3]
    print $1,$2,$3,gene_id,gene_name,gene_type,$6
  }' "$out_bed" \
    > "$out_summary"

  echo "[ok] BED:"
  echo "$out_bed"
  wc -l "$out_bed"

  echo "[ok] summary:"
  echo "$out_summary"
  wc -l "$out_summary"

  echo "gene type counts:"
  awk 'BEGIN{FS=OFS="\t"} {count[$6]++} END{for(k in count) print k,count[k]}' "$out_summary" \
    | sort -k1,1
}

# ------------------------------------------------------------
# 5. Run both all-gene and protein-coding-only outputs
# ------------------------------------------------------------
run_100pct_promoter_overlap \
  "all_gene_promoters_1kb_100pct_MAPQ10_LRonly" \
  "$ALL_PROM1KB" \
  "$ALL_COV" \
  "$ALL_OUTBED" \
  "$ALL_SUMMARY"

run_100pct_promoter_overlap \
  "protein_coding_promoters_1kb_100pct_MAPQ10_LRonly" \
  "$PC_PROM1KB" \
  "$PC_COV" \
  "$PC_OUTBED" \
  "$PC_SUMMARY"

echo
echo "All-gene outputs:"
echo "$ALL_OUTBED"
echo "$ALL_SUMMARY"
echo
echo "Protein-coding-only outputs:"
echo "$PC_OUTBED"
echo "$PC_SUMMARY"
