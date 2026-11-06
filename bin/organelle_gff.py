#!/usr/bin/env python3
"""Turn organelle annotations into one gene-level GFF3 with a ``gene_class`` attribute.

Two modes:

* ``normalise``: read the GFF written by an annotator (MITOS2, Chloe) and keep
  one top-level feature per gene.
* ``homology``: merge miniprot protein alignments, BLASTN rRNA hits and
  tRNAscan-SE tRNAs into non-overlapping best-hit gene models.

The output is meant for plotting and counting. The annotators' own files are
published alongside it and remain the reference annotation.
"""
import argparse
import re
import sys
from urllib.parse import unquote

# --------------------------------------------------------------------------
# Gene classes (order matters: first match wins)
# --------------------------------------------------------------------------
CLASSES = [
    ("tRNA", r"^(trn|trna|tRNA-)"),
    ("rRNA", r"^(rrn|rnl|rns|rrl|rrs|\d+(\.\d+)?s(_|$)|(large|small)_subunit|12s|16s|18s|23s|26s|28s|5s|4\.5s)"),
    ("Complex I (NADH dehydrogenase)", r"^(nad|nd\d|ndh)"),
    ("Complex II (succinate dehydrogenase)", r"^sdh"),
    ("Complex III (cytochrome b)", r"^(cob|cytb)"),
    ("Complex IV (cytochrome c oxidase)", r"^(cox|co\d)"),
    ("ATP synthase", r"^(atp|atp\d)"),
    ("Cytochrome c biogenesis", r"^(ccm|cc[mb])"),
    ("Ribosomal proteins", r"^(rps|rpl)"),
    ("Photosystem I", r"^psa"),
    ("Photosystem II", r"^psb"),
    ("Cytochrome b6/f", r"^pet"),
    ("RuBisCO", r"^rbc"),
    ("RNA polymerase", r"^rpo"),
    ("Maturases", r"^mat"),
    ("Conserved ORFs (ycf)", r"^ycf"),
]


def gene_class(name):
    low = name.lower()
    for cls, pattern in CLASSES:
        if re.match(pattern, low):
            return cls
    return "Other"


AA3 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
    "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
    "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "SeC": "U",
    "Sup": "X", "Undet": "X", "fMet": "M", "iMet": "M", "Ile2": "I",
}


# --------------------------------------------------------------------------
# GFF helpers
# --------------------------------------------------------------------------
def parse_attrs(text):
    """Attributes keyed in lower case: annotators disagree on Parent/parent, Name/name."""
    attrs = {}
    for part in text.strip().strip(";").split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            attrs[k.strip().lower()] = unquote(v.strip())
    return attrs


def read_gff(path):
    with open(path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                continue
            yield {
                "seqid": cols[0], "source": cols[1], "type": cols[2],
                "start": int(cols[3]), "end": int(cols[4]),
                "score": float(cols[5]) if cols[5] not in (".", "") else 0.0,
                "strand": cols[6], "attrs": parse_attrs(cols[8]),
            }


def overlaps(a, b, min_frac=0.3):
    if a["seqid"] != b["seqid"]:
        return False
    ov = min(a["end"], b["end"]) - max(a["start"], b["start"]) + 1
    if ov <= 0:
        return False
    shorter = min(a["end"] - a["start"], b["end"] - b["start"]) + 1
    return ov / shorter >= min_frac


def best_non_overlapping(features):
    """Greedy: keep the highest-scoring feature, drop what overlaps it on the same strand."""
    kept = []
    for f in sorted(features, key=lambda x: -x["score"]):
        if not any(k["strand"] == f["strand"] and overlaps(k, f) for k in kept):
            kept.append(f)
    return kept


def clean_name(name):
    name = name.split("|")[0]
    # Drop RefSeq copy/part suffixes such as nad6-2 or atp1-a1, keep tRNA anticodons
    if not name.lower().startswith("trn"):
        name = re.sub(r"-(a?\d+|ct)$", "", name)
    return name


def write_gff(features, source, out):
    out.write("##gff-version 3\n")
    counts = {}
    for f in sorted(features, key=lambda x: (x["seqid"], x["start"])):
        name = f["name"]
        counts[name] = counts.get(name, 0) + 1
        fid = f"{name}" if counts[name] == 1 else f"{name}_{counts[name]}"
        fid = re.sub(r"[^A-Za-z0-9_.()-]", "_", fid)
        cls = gene_class(name)
        attrs = f"ID={fid};Name={name};gene_class={cls.replace(';', ',')}"
        if f.get("note"):
            attrs += f";Note={f['note']}"
        out.write("\t".join([
            f["seqid"], source, f["type"], str(f["start"]), str(f["end"]),
            ".", f["strand"], ".", attrs,
        ]) + "\n")


# --------------------------------------------------------------------------
# Mode: normalise an annotator's GFF
# --------------------------------------------------------------------------
TOP_TYPES = {"gene", "tRNA", "rRNA", "ncRNA", "rRNA_gene", "tRNA_gene", "ncRNA_gene"}
CHILD_FALLBACK = {"CDS", "exon"}


def normalise(path):
    feats = list(read_gff(path))
    top = [f for f in feats if f["type"] in TOP_TYPES and "parent" not in f["attrs"]]
    if not top:  # annotators that only emit CDS/exon lines
        top = [f for f in feats if f["type"] in CHILD_FALLBACK]
    out = []
    for f in top:
        name = f["attrs"].get("name") or f["attrs"].get("gene") or f["attrs"].get("id", "unknown")
        name = clean_name(name)
        cls = gene_class(name)
        ftype = "tRNA" if cls == "tRNA" else "rRNA" if cls == "rRNA" else "gene"
        out.append({**f, "name": name, "type": ftype})
    # Collapse duplicates (e.g. gene + same-named CDS fallback)
    uniq = {}
    for f in out:
        key = (f["seqid"], f["name"], f["start"], f["end"], f["strand"])
        uniq.setdefault(key, f)
    return list(uniq.values())


# --------------------------------------------------------------------------
# Mode: homology merge
# --------------------------------------------------------------------------
def miniprot_genes(path, min_identity):
    genes = []
    for f in read_gff(path):
        if f["type"] != "mRNA":
            continue
        a = f["attrs"]
        if float(a.get("identity", 0)) < min_identity:
            continue
        name = clean_name(a.get("target", "unknown").split()[0])
        note = f"miniprot identity={a.get('identity')} target={a.get('target', '').split()[0]}"
        genes.append({**f, "name": name, "type": "gene", "note": note})
    return best_non_overlapping(genes)


def blast_rrnas(path, max_evalue, min_cov):
    hits = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            qseqid, sseqid, pident, length, qstart, qend = c[0], c[1], float(c[2]), int(c[3]), int(c[4]), int(c[5])
            sstart, send, evalue, bitscore, slen = int(c[6]), int(c[7]), float(c[8]), float(c[9]), int(c[10])
            if evalue > max_evalue or length < min_cov * slen:
                continue
            strand = "+" if send >= sstart else "-"
            name = clean_name(sseqid)
            name = {"large_subunit_ribosomal_RNA": "rrnL", "small_subunit_ribosomal_RNA": "rrnS"}.get(name, name)
            m = re.match(r"^(\d+(\.\d+)?)S_ribosomal_RNA$", name)
            if m:
                name = f"rrn{m.group(1)}"
            hits.append({
                "seqid": qseqid, "start": min(qstart, qend), "end": max(qstart, qend),
                "strand": strand, "score": bitscore, "name": name, "type": "rRNA",
                "note": f"blastn pident={pident:.1f} hit={sseqid}",
            })
    return best_non_overlapping(hits)


def trnascan_trnas(path):
    out = []
    for f in read_gff(path):
        if f["type"] not in ("tRNA", "pseudogene"):
            continue
        a = f["attrs"]
        iso = a.get("isotype", "Undet")
        anti = a.get("anticodon", "nnn").lower()
        name = f"trn{AA3.get(iso, 'X')}({anti})"
        note = "pseudogene" if f["type"] == "pseudogene" or a.get("pseudo") else None
        out.append({**f, "name": name, "type": "tRNA", "note": note})
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="mode", required=True)
    n = sub.add_parser("normalise")
    n.add_argument("--gff", required=True)
    n.add_argument("--source", required=True)
    h = sub.add_parser("homology")
    h.add_argument("--miniprot", required=True)
    h.add_argument("--blast", required=True)
    h.add_argument("--trnascan", required=True)
    h.add_argument("--min-identity", type=float, default=0.4)
    h.add_argument("--max-evalue", type=float, default=1e-10)
    h.add_argument("--min-rrna-cov", type=float, default=0.3)
    for p in (n, h):
        p.add_argument("--out", required=True)
    a = ap.parse_args()

    if a.mode == "normalise":
        feats, source = normalise(a.gff), a.source
    else:
        feats = (
            miniprot_genes(a.miniprot, a.min_identity)
            + blast_rrnas(a.blast, a.max_evalue, a.min_rrna_cov)
            + trnascan_trnas(a.trnascan)
        )
        source = "nanoplasty-homology"
    with open(a.out, "w") as out:
        write_gff(feats, source, out)
    print(f"{len(feats)} features written to {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
