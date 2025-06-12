#!/usr/bin/env bash
set -euo pipefail

# Configuration
# locate project root and resources
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
echo $ROOT
OUTPUT_FILE="$ROOT/manual/manual.pdf"
FILES=( "$ROOT/README.md" )
FILES+=( "$ROOT/INSTALL.md" )
FILES+=( "$ROOT/manual/"*.md )
INPUT_DIR="$PWD"
#BIB_FILE=( "$ROOT/manual/bib/references.bib" ) # Commented out
CSL_FILE=( "$ROOT/manual/bib/bioinformatics.csl" )

# Find all .bib files in the bib directory
BIB_FILES_FOUND=( "$ROOT"/manual/bib/*.bib )

# Prepare pandoc bibliography arguments
PANDOC_BIB_ARGS=()
if [ -e "${BIB_FILES_FOUND[0]}" ]; then # Check if any .bib files were actually found
  for bib_file in "${BIB_FILES_FOUND[@]}"; do
    PANDOC_BIB_ARGS+=( "--bibliography=$bib_file" )
  done
fi

#PDF_ENGINE="lualatex"
#PDF_ENGINE="xelatex"
PDF_ENGINE="pdflatex"

export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

# gfm for github style md
# highlight style theme for syntax
pandoc \
    --from=markdown \
    --to=pdf \
    --standalone \
    --embed-resources \
    --pdf-engine="$PDF_ENGINE" \
    --lua-filter="$ROOT/manual/rewrite.links.lua" \
    --include-in-header="$ROOT/manual/preamble.tex" \
    --toc \
    --toc-depth=2 \
    -V title="EPIHAPMPI manual" \
    -V author="Dzianis Prakapenka" \
    -V date="$(date '+%B %d, %Y')" \
    -V documentclass=book \
    -V classoption=oneside \
    -V classoption=openany \
    -V colorlinks=true \
    --top-level-division=chapter \
    "${PANDOC_BIB_ARGS[@]}" \
    --csl="$CSL_FILE" \
    --citeproc \
    --metadata link-citations=true \
    --highlight-style=zenburn \
    -V mainfont="Source Sans Pro" \
    -V monofont="Source Code Pro" \
    -V geometry:margin=1in \
    --output="$OUTPUT_FILE" \
    "${FILES[@]}"
#    "$ROOT/README.md" "$ROOT/INSTALL.md" "$ROOT/manual/03.running.md"

echo "PDF generated at $OUTPUT_FILE"

# TODO: maybe add custom template:
# pandoc ... --template=custom.tex ...
