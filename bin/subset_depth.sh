#!/usr/bin/env bash
# Keep only the depth-table rows for contigs present in a FASTA file
# (MetaBAT2 needs the same contigs in the same order in both files).
# Usage: subset_depth.sh contigs.fa depth.txt out_depth.txt
set -euo pipefail
fa="$1"; depth="$2"; out="$3"

awk '
    NR == FNR { if ($0 ~ /^>/) { name = substr($1, 2); keep[name] = 1 } ; next }
    FNR == 1  { print; next }
    ($1 in keep) { print }
' "$fa" "$depth" > "$out"
