#!/usr/bin/env python3
import gzip, re, sys
from bisect import bisect_left
import matplotlib.pyplot as plt

# settings 
UP = 5000
DOWN = 5000
NBIN_UP = 50
NBIN_BODY = 100
NBIN_DOWN = 50
TOTAL = 200
BIN_UP = UP // NBIN_UP      # 100bp
BIN_DOWN = DOWN // NBIN_DOWN # 100bp

BASE = "/BLUES/eric/ONT_WGBS/TSS_TES_plot"
GTF  = BASE + "/gencode.v45.annotation.gtf.gz"
BGDIR = BASE + "/bedgraph"

samples = [
  ("H9 primed",   BGDIR + "/H9_primed.CG.cov10.bg"),
  ("W3 primed",   BGDIR + "/W3_primed.CG.cov10.bg"),
  ("H9 naive",    BGDIR + "/H9_naive.CG.cov10.bg"),
  ("W3 naive",    BGDIR + "/W3_naive.CG.cov10.bg"),
  ("H9 TSC",      BGDIR + "/H9_TSC.CG.cov10.bg"),
  ("W3 TSC",      BGDIR + "/W3_TSC.CG.cov10.bg"),
]

# read genes (protein_coding; skip chrM/Y)
genes = []
attr_re = re.compile(r'(\S+)\s+"([^"]+)"')
with gzip.open(GTF, "rt") as f:
    for line in f:
        if line.startswith("#"): continue
        t = line.rstrip("\n").split("\t")
        if len(t) < 9: continue
        chrn, feat, start, end, strand, attrs = t[0], t[2], int(t[3])-1, int(t[4]), t[6], t[8]
        if feat != "gene": continue
        if chrn in ("chrM","chrY"): continue
        d = dict(attr_re.findall(attrs))
        gt = d.get("gene_type") or d.get("gene_biotype")
        if gt != "protein_coding": continue
        if end - start < 100: continue
        genes.append((chrn, start, end, strand))

print("genes:", len(genes))

# helpers 
def load_bg(path):
    # per chrom: positions[] and values[]
    pos = {}
    val = {}
    with open(path, "r") as f:
        for line in f:
            c, s, e, m = line.rstrip("\n").split("\t")
            s = int(s); m = float(m)
            if c not in pos:
                pos[c] = []
                val[c] = []
            pos[c].append(s)
            val[c].append(m)
    return pos, val

def mean_interval(p_list, v_list, a, b):
    # mean for CpGs with pos in [a,b); return None if empty
    i0 = bisect_left(p_list, a)
    i1 = bisect_left(p_list, b)
    if i1 <= i0: return None
    s = 0.0
    n = 0
    for i in range(i0, i1):
        s += v_list[i]
        n += 1
    return (s / n) if n else None

def compute_profile(pos, val):
    mD = [0.0] * TOTAL
    cnt = [0] * TOTAL

    for (chrn, start, end, strand) in genes:
        if chrn not in pos: continue
        P = pos[chrn]
        V = val[chrn]

        body_len = end - start
        k = body_len / float(NBIN_BODY)

        if strand == "+":
            # upstream
            for i in range(NBIN_UP):
                a = start - UP + i*BIN_UP
                b = a + BIN_UP
                if b <= 0: continue
                if a < 0: a = 0
                mv = mean_interval(P, V, int(a), int(b))
                if mv is None: continue
                mD[i] += mv; cnt[i] += 1

            # body (scaled)
            for i in range(NBIN_BODY):
                a = start + i*k
                b = start + (i+1)*k
                mv = mean_interval(P, V, int(a), int(b))
                if mv is None: continue
                idx = NBIN_UP + i
                mD[idx] += mv; cnt[idx] += 1

            # downstream
            for i in range(NBIN_DOWN):
                a = end + i*BIN_DOWN
                b = a + BIN_DOWN
                mv = mean_interval(P, V, int(a), int(b))
                if mv is None: continue
                idx = NBIN_UP + NBIN_BODY + i
                mD[idx] += mv; cnt[idx] += 1

        else:
            # "-" strand: flip orientation so output is still 5'->3'
            # upstream (5') is end -> end+UP
            for i in range(NBIN_UP):
                a = end + i*BIN_UP
                b = a + BIN_UP
                mv = mean_interval(P, V, int(a), int(b))
                if mv is None: continue
                mD[i] += mv; cnt[i] += 1

            # body reversed: bins go end->start
            for i in range(NBIN_BODY):
                a = end - (i+1)*k
                b = end - i*k
                mv = mean_interval(P, V, int(a), int(b))
                if mv is None: continue
                idx = NBIN_UP + i
                mD[idx] += mv; cnt[idx] += 1

            # downstream (3') is start-UP -> start
            for i in range(NBIN_DOWN):
                a = start - DOWN + i*BIN_DOWN
                b = a + BIN_DOWN
                if b <= 0: continue
                if a < 0: a = 0
                mv = mean_interval(P, V, int(a), int(b))
                if mv is None: continue
                idx = NBIN_UP + NBIN_BODY + i
                mD[idx] += mv; cnt[idx] += 1

    out = []
    for i in range(TOTAL):
        out.append(mD[i]/cnt[i] if cnt[i] else float("nan"))
    return out

#run all samples, write outputs, plot
profiles = {}
for (name, bg) in samples:
    print("loading:", name)
    pos, val = load_bg(bg)
    prof = compute_profile(pos, val)
    profiles[name] = prof
    with open(BASE + "/Out_mean_" + name.replace(" ","_") + ".txt", "w") as o:
        for v in prof:
            o.write(str(v) + "\n")

# plot
plt.figure(figsize=(10,5))
x = list(range(TOTAL))
for (name, _) in samples:
    plt.plot(x, profiles[name], linewidth=2, label=name)

plt.axvline(NBIN_UP, color="gray", linewidth=1)
plt.axvline(NBIN_UP + NBIN_BODY, color="gray", linewidth=1)

plt.ylim(0, 1.02)
plt.ylabel("DNA methylation", fontsize=18)
plt.xlabel("Gene", fontsize=24)
plt.xticks([0, NBIN_UP, NBIN_UP+NBIN_BODY, TOTAL-1], ["-5kb","TSS","TES","5kb"], fontsize=14)
plt.yticks(fontsize=12)
plt.legend(
    ncol=1,
    fontsize=11,
    frameon=False,
    loc="center left",
    bbox_to_anchor=(1.02, 0.5)
)
plt.title("Metagene methylation (±5kb, cov≥10)", fontsize=18, pad=10)
plt.tight_layout()
plt.savefig(BASE + "/Metagene_methylation_pm5kb_cov10.pdf")
plt.close()

print("done:", BASE + "/Metagene_methylation_pm5kb_cov10.pdf")
