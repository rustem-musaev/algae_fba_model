#!/usr/bin/env python3
"""
Split contigs into eukaryote / prokaryote / other FASTA files using the
DeepMicroClass score table. Each contig goes to the class with the highest score.

Writes:
  <prefix>.eukaryote.fa   (class Eukaryote)
  <prefix>.prokaryote.fa  (class Prokaryote)
  <prefix>.other.fa       (Plasmid, EukaryoteVirus, ProkaryoteVirus, or not scored)
  <prefix>.class_summary.tsv
"""
import argparse
import csv
import sys
from collections import Counter

CLASSES = ["Eukaryote", "EukaryoteVirus", "Plasmid", "Prokaryote", "ProkaryoteVirus"]


def read_scores(path):
    best = {}
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        cols = {name: i for i, name in enumerate(header) if name in CLASSES}
        if not cols:
            sys.exit(f"No class columns ({', '.join(CLASSES)}) found in {path}. Header was: {header}")
        for row in reader:
            if not row:
                continue
            name = row[0].strip()
            scores = {c: float(row[i]) for c, i in cols.items()}
            winner = max(scores, key=scores.get)
            best[name] = winner
            best.setdefault(name.split()[0], winner)  # also match on the first word
    return best


def fasta_records(path):
    header, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, seq
                header, seq = line[1:], []
            else:
                seq.append(line)
    if header is not None:
        yield header, seq


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--prefix", required=True)
    a = ap.parse_args()

    best = read_scores(a.tsv)
    out = {
        "Eukaryote": open(f"{a.prefix}.eukaryote.fa", "w"),
        "Prokaryote": open(f"{a.prefix}.prokaryote.fa", "w"),
        "other": open(f"{a.prefix}.other.fa", "w"),
    }
    counts = Counter()
    for header, seq in fasta_records(a.fasta):
        cls = best.get(header) or best.get(header.split()[0]) or "not_scored"
        counts[cls] += 1
        fh = out.get(cls, out["other"])
        fh.write(">" + header + "\n" + "\n".join(seq) + "\n")
    for fh in out.values():
        fh.close()

    with open(f"{a.prefix}.class_summary.tsv", "w") as fh:
        fh.write("class\tcontigs\n")
        for cls, n in sorted(counts.items()):
            fh.write(f"{cls}\t{n}\n")


if __name__ == "__main__":
    main()
