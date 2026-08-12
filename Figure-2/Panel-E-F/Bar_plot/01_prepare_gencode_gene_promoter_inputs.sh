#!/usr/bin/env bash
set -euo pipefail

FIGDIR="/BLUES/eric/ONT_WGBS/Figure_2"
OUTDIR="${FIGDIR}/Panel_E"

LRWIN="${FIGDIR}/LR_only_windows_100bp.bed"
GTF="/BLUES/eric/refs/D/gencode.v47.primary_assembly.annotation.gtf.gz"

mkdir -p "$OUTDIR"/{inputs,intermediate,outputs,scripts}

# ------------------------------------------------------------
# 1. Copy LR-only 100 bp windows into Panel_E
# ------------------------------------------------------------
cp -f "$LRWIN" "$OUTDIR/inputs/LR_only_windows_100bp.bed"

# ------------------------------------------------------------
# 2. Merge adjacent/overlapping LR-only 100 bp windows
# ------------------------------------------------------------
sort -k1,1 -k2,2n "$OUTDIR/inputs/LR_only_windows_100bp.bed" \
  | bedtools merge -i - \
  > "$OUTDIR/intermediate/LR_only_merged.bed"

echo "LR-only windows:"
wc -l "$OUTDIR/inputs/LR_only_windows_100bp.bed"

echo "Merged LR-only regions:"
wc -l "$OUTDIR/intermediate/LR_only_merged.bed"

# ------------------------------------------------------------
# 3. Make primary chromosome genome order file
# ------------------------------------------------------------
cat > "$OUTDIR/inputs/hg38.primary.chrom.sizes" <<'EOF'
chr1	248956422
chr2	242193529
chr3	198295559
chr4	190214555
chr5	181538259
chr6	170805979
chr7	159345973
chr8	145138636
chr9	138394717
chr10	133797422
chr11	135086622
chr12	133275309
chr13	114364328
chr14	107043718
chr15	101991189
chr16	90338345
chr17	83257441
chr18	80373285
chr19	58617616
chr20	64444167
chr21	46709983
chr22	50818468
chrX	156040895
chrY	57227415
EOF

# ------------------------------------------------------------
# 4. Keep only primary chromosomes in LR-only merged file
# ------------------------------------------------------------
awk 'BEGIN{OFS="\t"} $1 ~ /^chr([0-9]+|X|Y)$/ {print}' \
  "$OUTDIR/intermediate/LR_only_merged.bed" \
  | bedtools sort -g "$OUTDIR/inputs/hg38.primary.chrom.sizes" \
  > "$OUTDIR/intermediate/LR_only_merged.primary.sorted.bed"

echo "Merged LR-only primary-chromosome regions:"
wc -l "$OUTDIR/intermediate/LR_only_merged.primary.sorted.bed"

# ------------------------------------------------------------
# 5. Create clean GENCODE v47 gene BED
# ------------------------------------------------------------
python3 - <<'PY'
import gzip

gtf = "/BLUES/eric/refs/D/gencode.v47.primary_assembly.annotation.gtf.gz"
out = "/BLUES/eric/ONT_WGBS/Figure_2/Panel_E/intermediate/gencode_v47_genes.bed"

def parse_attrs(s):
    d = {}
    for item in s.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        parts = item.split(" ", 1)
        if len(parts) == 2:
            d[parts[0]] = parts[1].strip().strip('"')
    return d

with gzip.open(gtf, "rt") as f, open(out, "w") as o:
    for line in f:
        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue

        chrom, source, feature, start, end, score, strand, frame, attrs_raw = fields

        if feature != "gene":
            continue

        attrs = parse_attrs(attrs_raw)

        gene_id = attrs.get("gene_id", "NA")
        gene_name = attrs.get("gene_name", "NA")
        gene_type = attrs.get("gene_type", attrs.get("gene_biotype", "NA"))

        bed_start = int(start) - 1
        bed_end = int(end)

        name = f"{gene_id}|{gene_name}|{gene_type}"

        o.write(f"{chrom}\t{bed_start}\t{bed_end}\t{name}\t0\t{strand}\n")
PY

# Keep primary chromosomes and sort
awk 'BEGIN{OFS="\t"} $1 ~ /^chr([0-9]+|X|Y)$/ {print}' \
  "$OUTDIR/intermediate/gencode_v47_genes.bed" \
  | bedtools sort -g "$OUTDIR/inputs/hg38.primary.chrom.sizes" \
  > "$OUTDIR/intermediate/gencode_v47_genes.primary.sorted.bed"

echo "GENCODE v47 primary chromosome genes:"
wc -l "$OUTDIR/intermediate/gencode_v47_genes.primary.sorted.bed"

# ------------------------------------------------------------
# 6. Make promoter BEDs: +/-1 kb, +/-2 kb, +/-5 kb around TSS
# ------------------------------------------------------------
for W in 1000 2000 5000; do
  awk -v W="$W" 'BEGIN{FS=OFS="\t"}
  {
    chr=$1
    start=$2
    end=$3
    name=$4
    strand=$6

    if (strand == "+") {
      tss=start
    } else if (strand == "-") {
      tss=end
    } else {
      next
    }

    ps=tss-W
    pe=tss+W
    if (ps < 0) ps=0

    print chr,ps,pe,name,0,strand
  }' "$OUTDIR/intermediate/gencode_v47_genes.primary.sorted.bed" \
    | bedtools sort -g "$OUTDIR/inputs/hg38.primary.chrom.sizes" \
    > "$OUTDIR/intermediate/gencode_v47_promoters_${W}bp.primary.sorted.bed"

  echo "Promoters +/- ${W} bp:"
  wc -l "$OUTDIR/intermediate/gencode_v47_promoters_${W}bp.primary.sorted.bed"
done

echo
echo "Panel E files are in:"
echo "$OUTDIR"
