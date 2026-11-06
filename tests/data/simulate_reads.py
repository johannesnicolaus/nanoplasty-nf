#!/usr/bin/env python3
"""Simulate long reads from a circular genome (test data only).

Reads are random arcs of the circle, from either strand, with uniform
substitution/insertion/deletion errors. This is not a realistic error model;
it is just enough for NANOPlasty and minimap2 to have something to work with.

Usage: simulate_reads.py genome.fa out.fastq.gz --depth 60 --mean-len 6000 --error 0.04 --seed 1
"""
import argparse
import gzip
import random


def read_fasta(path):
    with open(path) as fh:
        return "".join(line.strip() for line in fh if not line.startswith(">")).upper()


def revcomp(s):
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def mutate(seq, rate, rng):
    out = []
    for base in seq:
        r = rng.random()
        if r < rate / 3:
            out.append(rng.choice("ACGT".replace(base, "") or "A"))
        elif r < 2 * rate / 3:
            out.append(base + rng.choice("ACGT"))
        elif r < rate:
            continue
        else:
            out.append(base)
    return "".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("genome")
    ap.add_argument("out")
    ap.add_argument("--depth", type=float, default=60)
    ap.add_argument("--mean-len", type=int, default=6000)
    ap.add_argument("--error", type=float, default=0.04)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()

    rng = random.Random(a.seed)
    genome = read_fasta(a.genome)
    n = len(genome)
    doubled = genome + genome
    total, i = 0, 0
    with gzip.open(a.out, "wt") as fh:
        while total < a.depth * n:
            length = min(int(rng.lognormvariate(0, 0.5) * a.mean_len), n - 1)
            if length < 500:
                continue
            start = rng.randrange(n)
            seq = doubled[start : start + length]
            if rng.random() < 0.5:
                seq = revcomp(seq)
            seq = mutate(seq, a.error, rng)
            fh.write(f"@sim_{i}\n{seq}\n+\n{'5' * len(seq)}\n")
            total += length
            i += 1


if __name__ == "__main__":
    main()
