#!/usr/bin/env bash
set -euo pipefail

OUTDIR="/BLUES/eric/ONT_WGBS/Figure_2/Panel_A"
CHROMSIZES="/BLUES/eric/WGBS/hg38.chrom.sizes"
GAP_TXT_GZ="$OUTDIR/hg38_gap.txt.gz"

# bedGraph.gz inputs: chr start end depth
WGBS_PRIMED="/BLUES/eric/ONT_WGBS/bedGraph/H9-primed.bedGraph.gz"
WGBS_NAIVE="/BLUES/eric/ONT_WGBS/bedGraph/H9-naive.bedGraph.gz"
WGBS_TSC="/BLUES/eric/ONT_WGBS/bedGraph/H9-naive-TSC.bedGraph.gz"

ONT_PRIMED="/BLUES/eric/ONT_WGBS/bedGraph/ONT_Primed_MAPQ10_primary.bedGraph.gz"
ONT_NAIVE="/BLUES/eric/ONT_WGBS/bedGraph/ONT_Naive_MAPQ10_primary.bedGraph.gz"
ONT_TSC="/BLUES/eric/ONT_WGBS/bedGraph/ONT_TSC_MAPQ10_primary.bedGraph.gz"

mkdir -p "$OUTDIR"

FAIDX_NOY="/BLUES/eric/ONT_WGBS/Figure_2/hg38.nochrY.chrom.sizes"
N_BED_RAW="/BLUES/eric/ONT_WGBS/Figure_2/N_only_noY.raw.bed"
N_BED="/BLUES/eric/ONT_WGBS/Figure_2/N_only_noY.merged.bed"
GENOME_BED="/BLUES/eric/ONT_WGBS/Figure_2/genome_noY.bed"
CALLABLE_BED="/BLUES/eric/ONT_WGBS/Figure_2/callable_noY.bed"

OUTTSV="$OUTDIR/Fig2A_right_genome_fraction_3cat_MAPQ10_primary.tsv"

run_lowprio() {
  if command -v ionice >/dev/null 2>&1; then
    ionice -c 3 nice -n 15 "$@"
  else
    nice -n 15 "$@"
  fi
}

echo "Building static/reference files if needed..."

if [[ ! -s "$FAIDX_NOY" ]]; then
  awk '$1!="chrY"' "$CHROMSIZES" > "$FAIDX_NOY"
fi

if [[ ! -s "$N_BED" ]]; then
  zcat "$GAP_TXT_GZ" \
    | awk -v OFS="\t" '
        NR==FNR {ok[$1]=1; next}
        ($2 in ok) && ($2!="chrY") {print $2,$3,$4}
      ' "$FAIDX_NOY" - \
    | bedtools sort -faidx "$FAIDX_NOY" -i - \
    > "$N_BED_RAW"

  bedtools merge -i "$N_BED_RAW" > "$N_BED"
fi

if [[ ! -s "$GENOME_BED" ]]; then
  awk -v OFS="\t" '{print $1,0,$2}' "$FAIDX_NOY" > "$GENOME_BED"
fi

if [[ ! -s "$CALLABLE_BED" ]]; then
  bedtools subtract -a "$GENOME_BED" -b "$N_BED" \
    | bedtools sort -faidx "$FAIDX_NOY" -i - \
    > "$CALLABLE_BED"
fi

TOTAL_GENOME_BASES=$(awk '{s+=$2} END{print s+0}' "$FAIDX_NOY")
TOTAL_CALLABLE_BASES=$(awk '{s+=($3-$2)} END{print s+0}' "$CALLABLE_BED")
TOTAL_N_BASES=$(( TOTAL_GENOME_BASES - TOTAL_CALLABLE_BASES ))

# Initialize output table
echo -e "state\ttechnology\tcategory\tbases\tpercent_genome" > "$OUTTSV"

detected_callable_bases_ge1 () {
  local BG_GZ="$1"

  run_lowprio bash -c '
    set -euo pipefail

    FAIDX="$1"
    CALLABLE="$2"
    BG="$3"

    zcat "$BG" \
      | awk -v OFS="\t" "
          NR==FNR {ok[\$1]=1; next}
          /^track/ || /^browser/ || /^#/ || NF<4 {next}
          !(\$1 in ok) {next}
          \$1==\"chrY\" {next}
          \$2 !~ /^[0-9]+$/ {next}
          \$3 !~ /^[0-9]+$/ {next}
          \$3 <= \$2 {next}
          (\$4+0) < 1 {next}
          {print \$1,\$2,\$3}
        " "$FAIDX" - \
      | bedtools intersect -a "$CALLABLE" -b - -wo \
      | awk "{sum += \$NF} END{printf \"%.0f\\n\", sum+0}"
  ' _ "$FAIDX_NOY" "$CALLABLE_BED" "$BG_GZ"
}

emit_rows () {
  local STATE="$1"
  local TECH="$2"
  local BG="$3"

  echo "Processing: $STATE $TECH"
  local DET
  DET=$(detected_callable_bases_ge1 "$BG")

  local UND=$(( TOTAL_CALLABLE_BASES - DET ))
  local ND="$TOTAL_N_BASES"

  awk -v OFS="\t" \
      -v st="$STATE" \
      -v te="$TECH" \
      -v tot="$TOTAL_GENOME_BASES" \
      -v nd="$ND" \
      -v und="$UND" \
      -v det="$DET" '
    BEGIN{
      printf "%s\t%s\tNon-detectable (N-only)\t%d\t%.6f\n", st,te,nd,100*nd/tot;
      printf "%s\t%s\tUndetected (<1×, non-N)\t%d\t%.6f\n", st,te,und,100*und/tot;
      printf "%s\t%s\tDetected (≥1×, non-N)\t%d\t%.6f\n", st,te,det,100*det/tot;
    }' >> "$OUTTSV"
}

echo "Running Panel A summary..."

emit_rows Primed WGBS "$WGBS_PRIMED"
emit_rows Primed ONT  "$ONT_PRIMED"
emit_rows Naive  WGBS "$WGBS_NAIVE"
emit_rows Naive  ONT  "$ONT_NAIVE"
emit_rows TSC    WGBS "$WGBS_TSC"
emit_rows TSC    ONT  "$ONT_TSC"

echo "Done:"
echo "$OUTTSV"