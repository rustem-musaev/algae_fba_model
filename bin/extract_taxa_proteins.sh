#!/usr/bin/env bash
# For each taxon name in a list (one per line), pull out the proteins that
# EUKulele assigned to it. Same idea as the loop in the original protocol,
# but without needing seqkit.
# Usage: extract_taxa_proteins.sh estimated-taxonomy.tsv proteins.faa taxa.txt
# Output: <Taxon_name>_proteins.faa (spaces become underscores)
set -euo pipefail
tax="$1"; prot="$2"; taxa="$3"

# find the column with the protein names (called transcript_name by EUKulele; fall back to column 2)
col=$(head -n 1 "$tax" | awk -F'\t' '{ for (i = 1; i <= NF; i++) if ($i == "transcript_name") { print i; exit } }')
col=${col:-2}

while IFS= read -r taxon || [ -n "$taxon" ]; do
    taxon="${taxon%$'\r'}"
    [ -z "$taxon" ] && continue
    outname="${taxon// /_}_proteins.faa"
    echo "Processing: $taxon" >&2

    grep -F -- "$taxon" "$tax" | awk -F'\t' -v c="$col" '{ print $c }' | sort -u > ids.txt || true
    if [ ! -s ids.txt ]; then
        echo "  no proteins assigned to '$taxon'" >&2
        continue
    fi

    awk 'NR == FNR { want[$1] = 1; next }
         /^>/ { name = substr($1, 2); keep = (name in want) }
         keep' ids.txt "$prot" > "$outname"
    if [ ! -s "$outname" ]; then
        echo "  ids found but no matching proteins for '$taxon'" >&2
        rm -f "$outname"
        continue
    fi
    echo "  $(grep -c '^>' "$outname") proteins -> $outname" >&2
done < "$taxa"
rm -f ids.txt
