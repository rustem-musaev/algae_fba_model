#!/usr/bin/env bash
# Turn a GeneMark-ES GTF into an AUGUSTUS hints file (same rules as the original
# protocol): gene/mRNA -> genicpart, any UTR -> UTR, start_codon -> start,
# stop_codon -> stop, and column 9 becomes "source=M".
# Usage: fix_genemark_gtf.sh genemark.gtf hints.gtf
set -euo pipefail
in="$1"; out="$2"

awk 'BEGIN { FS = OFS = "\t" }
    /^#/ || NF < 9 { next }
    {
        if ($3 == "gene" || $3 == "mRNA") $3 = "genicpart"
        if (index($3, "UTR") > 0)         $3 = "UTR"
        if ($3 == "start_codon")          $3 = "start"
        if ($3 == "stop_codon")           $3 = "stop"
        $9 = "source=M"
        print
    }' "$in" > "$out"
