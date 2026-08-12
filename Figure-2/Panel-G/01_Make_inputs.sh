#!/usr/bin/env bash
set -euo pipefail

THREADS=8

FIGDIR="/BLUES/eric/ONT_WGBS/Figure_2"
OUTDIR="${FIGDIR}/Panel_D"

LRWIN="${FIGDIR}/LR_only_MAPQ/Primed_LR_only_MAPQ10_windows_100bp.bed"
TEGZ="/BLUES/eric/ONT/hg38_TE_noY.bed.gz"

TOP10="${OUTDIR}/inputs/top10_subfamilies_from_pie.txt"

mkdir -p "$OUTDIR"/{inputs,intermediate,outputs,scripts,perchr}

# Inputs
cp -f "$LRWIN" "$OUTDIR/inputs/LR_only_windows_100bp.bed"

# Merge LR-only windows into regions
sort -k1,1 -k2,2n "$OUTDIR/inputs/LR_only_windows_100bp.bed" \
  | bedtools merge -i - \
  > "$OUTDIR/intermediate/LR_only_merged.bed"

# Parse TE annotation: chr start end strand subfamily class family
zcat "$TEGZ" \
| awk 'BEGIN{OFS="\t"}
{
  chr=$1; st=$2; en=$3; strand=$4;
  name="NA"; cls="NA"; fam="NA";

  if (match($0, /repeat_id "[^"]+"/)) {
    rid=substr($0, RSTART, RLENGTH);
    gsub(/repeat_id "|"/, "", rid);
    split(rid, a, ",");
    name=a[1];
  }
  if (match($0, /repeat_class "[^"]+"/)) {
    rc=substr($0, RSTART, RLENGTH);
    gsub(/repeat_class "|"/, "", rc);
    cls=rc;
  }
  if (match($0, /repeat_family "[^"]+"/)) {
    rf=substr($0, RSTART, RLENGTH);
    gsub(/repeat_family "|"/, "", rf);
    fam=rf;
  }

  print chr, st, en, strand, name, cls, fam;
}' > "$OUTDIR/intermediate/hg38_TE_noY.parsed.bed"

# Sort parsed TE (bedtools expects sorted for some operations)
sort -k1,1 -k2,2n "$OUTDIR/intermediate/hg38_TE_noY.parsed.bed" \
> "$OUTDIR/intermediate/hg38_TE_noY.parsed.sorted.bed"

# Updated top 10 from MAPQ10 Fig. 2C
cat > "$TOP10" << 'EOF'
L1PA2
L1PA3
L1HS
AluSx
AluY
AluSx1
L1PA4
AluSz
AluJr
L2a
EOF

# Denominator: total copies genome-wide per subfamily (restricted to TOP10) 
awk 'BEGIN{OFS="\t"}
FNR==NR{wanted[$1]=1; next}
{
  sf=$5; cls=$6; fam=$7;
  if (wanted[sf]) {
    total[sf]+=1;
    scls[sf]=cls; sfam[sf]=fam;
  }
}
END{
  for (sf in total) print sf, scls[sf], sfam[sf], total[sf];
}' "$TOP10" "$OUTDIR/intermediate/hg38_TE_noY.parsed.sorted.bed" \
> "$OUTDIR/intermediate/top10_total_genome_copies.tsv"
# cols: subfamily class family total_copies_genome

# Numerator: copies overlapping LR-only region (count each TE record once) 
# Intersect TE (as A) with LR-only (as B) and keep unique A records
bedtools intersect -a "$OUTDIR/intermediate/hg38_TE_noY.parsed.sorted.bed" \
  -b "$OUTDIR/intermediate/LR_only_merged.bed" -u \
> "$OUTDIR/intermediate/TE_overlapping_LRonly.bed"

awk 'BEGIN{OFS="\t"}
FNR==NR{wanted[$1]=1; next}
{
  sf=$5;
  if (wanted[sf]) inlr[sf]+=1;
}
END{
  for (sf in inlr) print sf, inlr[sf];
}' "$TOP10" "$OUTDIR/intermediate/TE_overlapping_LRonly.bed" \
> "$OUTDIR/intermediate/top10_LRonly_copies.tsv"
# cols: subfamily copies_in_LRonly

# Join + compute fraction captured 
awk 'BEGIN{OFS="\t"}
FNR==NR{inlr[$1]=$2; next}
{
  sf=$1; cls=$2; fam=$3; tot=$4;
  lr = (sf in inlr) ? inlr[sf] : 0;
  pct = (tot>0) ? (100.0*lr/tot) : 0;
  print sf, cls, fam, lr, tot, pct;
}' "$OUTDIR/intermediate/top10_LRonly_copies.tsv" \
   "$OUTDIR/intermediate/top10_total_genome_copies.tsv" \
> "$OUTDIR/outputs/Fig2D_top10_fractionCaptured.tsv"

echo "Saved: $OUTDIR/outputs/Fig2D_top10_fractionCaptured.tsv"
