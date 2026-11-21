#!/usr/bin/env bash
# Put the bin name in front of every protein name (g1.t1 -> BIN__g1.t1)
# and remove a trailing stop "*" (InterProScan does not accept it).
# Usage: tag_proteins.sh in.aa out.faa BIN
set -euo pipefail
in="$1"; out="$2"; bin="$3"

awk -v bin="$bin" '
    /^>/ { sub(/^>/, ">" bin "__"); print; next }
    { sub(/\*$/, ""); print }
' "$in" > "$out"
