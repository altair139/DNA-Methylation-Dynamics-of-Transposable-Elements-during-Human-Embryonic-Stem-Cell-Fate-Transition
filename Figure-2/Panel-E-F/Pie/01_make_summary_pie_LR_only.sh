#!/usr/bin/env bash
set -euo pipefail

THREADS=8 

PANEL_C="/BLUES/eric/ONT_WGBS/Figure_2/Panel_C"
LR_BED="/BLUES/eric/ONT_WGBS/Figure_2/LR_only_MAPQ/Primed_LR_only_MAPQ10_windows_100bp.bed"
TE_GZ="/BLUES/eric/ONT/hg38_TE_noY.bed.gz"

mkdir -p "$PANEL_C"/{inputs,intermediate,outputs,scripts}

command -v bedtools >/dev/null || { echo "ERROR: bedtools not in PATH" >&2; exit 1; }
command -v zcat >/dev/null || { echo "ERROR: zcat not in PATH" >&2; exit 1; }

[[ -s "$LR_BED" ]] || { echo "ERROR: missing LR-only bed: $LR_BED" >&2; exit 1; }
[[ -s "$TE_GZ"  ]] || { echo "ERROR: missing TE gz: $TE_GZ" >&2; exit 1; }

echo "[info] PANEL_C = $PANEL_C"
echo "[info] LR_BED  = $LR_BED"
echo "[info] TE_GZ   = $TE_GZ"

# ---- 1) Sort LR-only windows ----
LR_SORTED="$PANEL_C/intermediate/LR_only_windows_100bp.sorted.bed"
bedtools sort -i "$LR_BED" > "$LR_SORTED"
echo "[ok] LR sorted: $LR_SORTED"

# ---- 2) Parse TE annotation ----
# Expected input:
# chr start end strand attrs
#
# attrs contains:
# repeat_id "L1MC5a,chr1,11485,11676";
# repeat_class "LINE";
# repeat_family "L1";
#
# Output:
# chr start end strand subfamily class family

TE_PARSED="$PANEL_C/intermediate/TE.parsed.bed"

echo "[info] Parsing TE annotation -> $TE_PARSED"

zcat "$TE_GZ" |
awk '
BEGIN {
    FS = "\t"
    OFS = "\t"
}

# Skip blank lines and comments
NF == 0 || $0 ~ /^#/ {
    next
}

{
    chr    = $1
    start  = $2
    end    = $3
    strand = $4

    subfam = "NA"
    teclass = "NA"
    fam = "NA"

    if (NF >= 7 && $5 !~ /repeat_id/) {
        subfam  = $5
        teclass = $6
        fam     = $7
    } else {
        attrs = $5

        # repeat_id "L1MC5a,chr1,11485,11676"
        if (match(attrs, /repeat_id "[^"]+"/)) {
            value = substr(attrs, RSTART, RLENGTH)
            sub(/^repeat_id "/, "", value)
            sub(/"$/, "", value)

            split(value, fields, ",")
            subfam = fields[1]
        }

        # repeat_class "LINE"
        if (match(attrs, /repeat_class "[^"]+"/)) {
            value = substr(attrs, RSTART, RLENGTH)
            sub(/^repeat_class "/, "", value)
            sub(/"$/, "", value)
            teclass = value
        }

        # repeat_family "L1"
        if (match(attrs, /repeat_family "[^"]+"/)) {
            value = substr(attrs, RSTART, RLENGTH)
            sub(/^repeat_family "/, "", value)
            sub(/"$/, "", value)
            fam = value
        }
    }

    print chr, start, end, strand, subfam, teclass, fam
}
' > "$TE_PARSED"

TE_SORTED="$PANEL_C/intermediate/TE.parsed.sorted.bed"
bedtools sort -i "$TE_PARSED" > "$TE_SORTED"

# ---- 3) LR-only windows with NO TE overlap ----
# -v means "report A entries that have no overlap with B"
LR_NO_TE="$PANEL_C/outputs/LR_only_windows_100bp.noTE.bed"
bedtools intersect -v -a "$LR_SORTED" -b "$TE_SORTED" > "$LR_NO_TE"
echo "[ok] Wrote: $LR_NO_TE"

# ---- 4) Merge adjacent/overlapping windows into larger regions ----
LR_NO_TE_MERGED="$PANEL_C/outputs/LR_only_regions.noTE.merged.bed"
bedtools merge -i "$LR_NO_TE" > "$LR_NO_TE_MERGED"
echo "[ok] Wrote: $LR_NO_TE_MERGED"

# ---- 5) Paste results (quick check) ----
TOTAL_WIN=$(wc -l < "$LR_SORTED" | awk "{print \$1}")
NOTE_WIN=$(wc -l < "$LR_NO_TE" | awk "{print \$1}")
echo "[stats] LR-only windows total: $TOTAL_WIN"
echo "[stats] LR-only windows with NO TE: $NOTE_WIN"
echo "[stats] fraction no-TE: $(awk -v a="$NOTE_WIN" -v t="$TOTAL_WIN" 'BEGIN{print (t? a/t : 0)}')"

echo "[done]"

