#!/usr/bin/env python3
"""Circular map of an organelle assembly, plus a one-row QC summary.

Tracks, from the outside in:
  1. genes coloured by functional class (forward strand outside, reverse inside)
  2. read depth (from mosdepth, folded back from the wrapped reference)
  3. GC content, as deviation from the genome mean
  4. GC skew (G-C)/(G+C)

Usage:
  plot_organelle.py --fasta asm.fa --gff asm.gff3 --depth wrapped.per-base.bed.gz \
      --junction junction.tsv --sample S1 --type chloro --genome-range 120000-200000 \
      --prefix S1
"""
import argparse
import gzip
import statistics
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from pycirclize import Circos  # noqa: E402

# Colour per gene class (OGDRAW-like where it has a convention)
CLASS_COLOURS = {
    "Complex I (NADH dehydrogenase)": "#e6c229",
    "Complex II (succinate dehydrogenase)": "#9b7653",
    "Complex III (cytochrome b)": "#7cae7a",
    "Complex IV (cytochrome c oxidase)": "#e88f8f",
    "ATP synthase": "#8fb3e8",
    "Cytochrome c biogenesis": "#c9a0dc",
    "Ribosomal proteins": "#c08552",
    "Photosystem I": "#2e7d32",
    "Photosystem II": "#66bb6a",
    "Cytochrome b6/f": "#5c8a8a",
    "RuBisCO": "#1b5e20",
    "RNA polymerase": "#b71c1c",
    "Maturases": "#ff8a65",
    "Conserved ORFs (ycf)": "#80cbc4",
    "tRNA": "#212121",
    "rRNA": "#d32f2f",
    "Other": "#9e9e9e",
}


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line.upper())
    return {k: "".join(v) for k, v in seqs.items()}


def read_gff(path):
    feats = defaultdict(list)
    if not path:
        return feats
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            attrs = dict(kv.split("=", 1) for kv in c[8].split(";") if "=" in kv)
            feats[c[0]].append({
                "start": int(c[3]) - 1, "end": int(c[4]), "strand": c[6],
                "name": attrs.get("Name", ""), "cls": attrs.get("gene_class", "Other"),
            })
    return feats


def read_depth(path, lengths):
    """Per-base depth on the original coordinates.

    The reads were mapped to each contig with its first W bases appended, so
    depth at wrapped position L+i belongs to position i.
    """
    depth = {k: np.zeros(n, dtype=np.int64) for k, n in lengths.items()}
    if not path:
        return depth
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            chrom, s, e, d = line.split("\t")
            if chrom not in depth:
                continue
            n = lengths[chrom]
            s, e, d = int(s), int(e), int(d)
            if d == 0:
                continue
            for a, b in ((s, min(e, n)), (max(s, n) - n, e - n)):
                if b > a:
                    depth[chrom][a:b] += d
    return depth


def read_junction(path):
    out = {}
    if path:
        with open(path) as fh:
            for line in fh:
                if line.startswith("contig"):
                    continue
                c = line.rstrip("\n").split("\t")
                out[c[0]] = int(c[2])
    return out


def windows(seq, size):
    xs, gc, skew = [], [], []
    for i in range(0, len(seq), size):
        w = seq[i : i + size]
        g, c = w.count("G"), w.count("C")
        acgt = sum(w.count(b) for b in "ACGT") or 1
        xs.append(i + len(w) / 2)
        gc.append((g + c) / acgt)
        skew.append((g - c) / (g + c) if g + c else 0.0)
    return np.array(xs), np.array(gc), np.array(skew)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--gff")
    ap.add_argument("--depth")
    ap.add_argument("--junction")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--type", default="")
    ap.add_argument("--genome-range", default="")
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--window", type=int, default=0, help="GC window (default: length/360, min 50)")
    a = ap.parse_args()

    seqs = read_fasta(a.fasta)
    lengths = {k: len(v) for k, v in seqs.items()}
    total = sum(lengths.values())
    feats = read_gff(a.gff)
    depth = read_depth(a.depth, lengths)
    junction = read_junction(a.junction)
    window = a.window or max(50, total // 360)

    all_seq = "".join(seqs.values())
    gc_total = sum(all_seq.count(b) for b in "GC") / max(1, sum(all_seq.count(b) for b in "ACGT"))
    all_depth = np.concatenate(list(depth.values())) if depth else np.array([0])
    mean_depth = float(all_depth.mean()) if all_depth.size else 0.0

    # ---------------- plot ----------------
    circos = Circos(lengths, space=0 if len(lengths) == 1 else 2)
    tick = 10 ** max(3, len(str(total)) - 1)
    if total / tick < 4:
        tick //= 2
    classes_seen = set()
    depth_max = max(1, float(np.percentile(all_depth, 99.5))) if all_depth.size else 1

    for sector in circos.sectors:
        name, n = sector.name, lengths[sector.name]
        ruler = sector.add_track((99, 100))
        ruler.axis(fc="black")
        ruler.xticks_by_interval(
            tick, label_size=6, label_orientation="vertical",
            label_formatter=lambda v: f"{v / 1000:.0f} kb",
        )

        fwd = sector.add_track((94, 98.5))
        rev = sector.add_track((89.5, 94))
        for f in feats.get(name, []):
            colour = CLASS_COLOURS.get(f["cls"], CLASS_COLOURS["Other"])
            classes_seen.add(f["cls"])
            track = fwd if f["strand"] != "-" else rev
            track.rect(f["start"], min(f["end"], n), fc=colour, ec="none")
        # Gene names, placed just outside the ruler with leader lines
        labelled = sorted(feats.get(name, []), key=lambda f: f["start"])
        if labelled and len(labelled) <= 250:
            for f in labelled:
                mid = (f["start"] + min(f["end"], n)) / 2
                ruler.annotate(mid, f["name"], label_size=5, min_r=101, max_r=112)

        cov = sector.add_track((74, 88), r_pad_ratio=0.05)
        cov.axis(fc="#f5f5f5", ec="none")
        step = max(1, n // 2000)
        xs = np.arange(0, n, step)
        ys = depth[name][::step].astype(float)
        cov.fill_between(xs, np.minimum(ys, depth_max), vmin=0, vmax=depth_max, fc="#5c6bc0", ec="none")

        xw, gc, skew = windows(seqs[name], window)
        gct = sector.add_track((60, 72))
        dev = gc - gc_total
        lim = max(0.01, float(np.abs(dev).max()))
        gct.fill_between(xw, np.where(dev > 0, dev, 0), y2=0, vmin=-lim, vmax=lim, fc="#555555", ec="none")
        gct.fill_between(xw, np.where(dev < 0, dev, 0), y2=0, vmin=-lim, vmax=lim, fc="#aaaaaa", ec="none")

        skt = sector.add_track((46, 58))
        slim = max(0.01, float(np.abs(skew).max()))
        skt.fill_between(xw, np.where(skew > 0, skew, 0), y2=0, vmin=-slim, vmax=slim, fc="#43a047", ec="none")
        skt.fill_between(xw, np.where(skew < 0, skew, 0), y2=0, vmin=-slim, vmax=slim, fc="#e53935", ec="none")

    centre = (
        f"{a.sample}\n{total:,} bp\nGC {gc_total * 100:.1f}%\n"
        f"mean depth {mean_depth:.0f}x"
    )
    circos.text(centre, size=11)
    fig = circos.plotfig(figsize=(10, 10))
    handles = [Patch(fc=CLASS_COLOURS[c], label=c) for c in CLASS_COLOURS if c in classes_seen]
    handles += [
        Patch(fc="#5c6bc0", label=f"Read depth (0-{depth_max:.0f}x)"),
        Patch(fc="#555555", label="GC above mean"),
        Patch(fc="#aaaaaa", label="GC below mean"),
        Patch(fc="#43a047", label="GC skew +"),
        Patch(fc="#e53935", label="GC skew -"),
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.02), fontsize=7, frameon=False, ncol=4)
    fig.savefig(f"{a.prefix}.circos.png", dpi=200, bbox_inches="tight")
    fig.savefig(f"{a.prefix}.circos.svg", bbox_inches="tight")
    plt.close(fig)

    # ---------------- summary row ----------------
    in_range = "NA"
    if a.genome_range:
        lo, hi = (int(x) for x in a.genome_range.split("-"))
        in_range = "yes" if lo <= total <= hi else "no"
    min_win = 0.0
    if all_depth.size and total >= 1000:
        k = 500
        cs = np.convolve(all_depth, np.ones(k) / k, mode="valid")
        min_win = float(cs.min())
    n_by = defaultdict(int)
    for fl in feats.values():
        for f in fl:
            n_by["tRNA" if f["cls"] == "tRNA" else "rRNA" if f["cls"] == "rRNA" else "protein"] += 1
    span = [junction.get(k, 0) for k in lengths]
    row = {
        "Sample": a.sample,
        "Type": a.type,
        "Contigs": len(lengths),
        "Length (bp)": total,
        "In expected range": in_range,
        "GC (%)": round(gc_total * 100, 2),
        "Mean depth": round(mean_depth, 1),
        "Median depth": round(float(statistics.median(all_depth)) if all_depth.size else 0, 1),
        "Min depth (500 bp window)": round(min_win, 1),
        "Reads spanning start/end": min(span) if span else 0,
        "Protein genes": n_by["protein"],
        "tRNAs": n_by["tRNA"],
        "rRNAs": n_by["rRNA"],
    }
    with open(f"{a.prefix}.summary.tsv", "w") as fh:
        fh.write("\t".join(row) + "\n")
        fh.write("\t".join(str(v) for v in row.values()) + "\n")


if __name__ == "__main__":
    main()
