#!/usr/bin/env bash
# Keep only contigs at least MIN_LEN bases long.
# Usage: filter_contigs.sh in.fa out.fa MIN_LEN
set -euo pipefail
in="$1"; out="$2"; min="$3"

awk -v min="$min" '
    function flush() { if (h != "" && length(s) >= min) { print h; print s } }
    /^>/ { flush(); h = $0; s = ""; next }
    { s = s $0 }
    END { flush() }
' "$in" > "$out"
