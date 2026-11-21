#!/usr/bin/env bash
# Give every contig a short unique name (c000001, c000002, ...) so RepeatMasker
# accepts it (it refuses headers longer than 50 characters).
# The old and new names are written to a two-column table.
# Usage: rename_headers.sh in.fa out.fa map.tsv
set -euo pipefail
in="$1"; out="$2"; map="$3"

awk -v map="$map" '
    BEGIN { print "new_name\toriginal_header" > map }
    /^>/ {
        n++
        new = sprintf("c%06d", n)
        print new "\t" substr($0, 2) > map
        print ">" new
        next
    }
    { print }
' "$in" > "$out"
