#!/usr/bin/env python3
"""Build the organelle reference sets used by the homology annotator.

Downloads RefSeq plastid and plant mitochondrial genomes from NCBI and writes,
for each organelle, a protein FASTA (for miniprot) and an rRNA FASTA (for
BLASTN). Headers are rewritten to ``<gene>|<accession>`` so the merger can
name features from the best hit.

Usage: python3 build_refs.py [outdir]    (default: directory of this script)
"""
import re
import sys
import time
import urllib.request
from pathlib import Path

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={acc}&rettype={rt}&retmode=text"

SETS = {
    "chloro": [
        "NC_000932.1",  # Arabidopsis thaliana
        "NC_001879.2",  # Nicotiana tabacum
        "NC_001666.2",  # Zea mays
        "NC_031333.1",  # Oryza sativa
        "NC_005086.1",  # Amborella trichopoda
        "NC_001631.1",  # Pinus thunbergii
        "NC_005087.2",  # Physcomitrium patens
    ],
    "plant_mito": [
        "NC_037304.1",  # Arabidopsis thaliana
        "NC_006581.1",  # Nicotiana tabacum
        "NC_007982.1",  # Zea mays
        "NC_021152.1",  # Liriodendron tulipifera
        "NC_010303.1",  # Cycas taitungensis
        "NC_037508.1",  # Marchantia polymorpha
    ],
}

# Lineage-specific ORFs add noise rather than annotation
SKIP_GENE = re.compile(r"^(orf|ORF|hypothetical)", re.I)


def fetch(acc, rettype):
    for attempt in range(5):
        try:
            with urllib.request.urlopen(EUTILS.format(acc=acc, rt=rettype), timeout=120) as r:
                return r.read().decode()
        except Exception:  # noqa: BLE001 - NCBI is flaky; retry
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"Could not fetch {acc} ({rettype})")


def parse_fasta(text):
    name, seq = None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(seq)
            name, seq = line[1:], []
        elif line.strip():
            seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def revcomp(s):
    return s[::-1].translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))


def proteins(acc):
    for header, seq in parse_fasta(fetch(acc, "fasta_cds_aa")):
        m = re.search(r"\[gene=([^\]]+)\]", header)
        gene = m.group(1) if m else None
        if gene is None or SKIP_GENE.match(gene) or "pseudo=true" in header:
            continue
        yield f"{gene}|{acc}", seq.rstrip("*")


def rrnas(acc):
    genome = "".join(s for _, s in parse_fasta(fetch(acc, "fasta")))
    feats, cur = [], None
    for line in fetch(acc, "ft").splitlines():
        cols = line.split("\t")
        if len(cols) >= 3 and cols[0] and cols[1] and cols[2]:
            cur = {"type": cols[2], "spans": [(cols[0], cols[1])], "name": None}
            feats.append(cur)
        elif len(cols) >= 2 and cols[0] and cols[1] and cur is not None and len(cols) < 3:
            cur["spans"].append((cols[0], cols[1]))
        elif len(cols) >= 5 and cur is not None and cols[3] in ("gene", "product") and cur["name"] is None:
            cur["name"] = cols[4]
    for f in feats:
        if f["type"] != "rRNA" or not f["name"]:
            continue
        seq = ""
        for a, b in f["spans"]:
            a, b = int(a.lstrip("<>")), int(b.lstrip("<>"))
            seq += genome[a - 1 : b] if a <= b else revcomp(genome[b - 1 : a])
        name = re.sub(r"\s+", "_", f["name"])
        yield f"{name}|{acc}", seq


def write(path, records):
    seen = set()
    with open(path, "w") as fh:
        for name, seq in records:
            key = (name, seq)
            if key in seen or not seq:
                continue
            seen.add(key)
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i : i + 80] + "\n")


def main():
    outdir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent
    outdir.mkdir(parents=True, exist_ok=True)
    for organelle, accs in SETS.items():
        prot, rrna = [], []
        for acc in accs:
            prot.extend(proteins(acc))
            rrna.extend(rrnas(acc))
            time.sleep(0.5)
        write(outdir / f"{organelle}_proteins.faa", prot)
        write(outdir / f"{organelle}_rrna.fna", rrna)
        print(f"{organelle}: {len(prot)} proteins, {len(rrna)} rRNAs", file=sys.stderr)


if __name__ == "__main__":
    main()
