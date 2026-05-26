#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ARTICLE_DIR="$ROOT_DIR/article"
MAIN_TEX="$ARTICLE_DIR/main.tex"
JOB_NAME="main"

if [[ ! -f "$MAIN_TEX" ]]; then
  echo "error: expected article source at $MAIN_TEX" >&2
  exit 1
fi

require_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "error: required command not found: $1" >&2
    exit 1
  fi
}

require_cmd pdflatex
require_cmd bibtex

cd "$ARTICLE_DIR"

pdflatex -interaction=nonstopmode -halt-on-error "$JOB_NAME.tex"
bibtex "$JOB_NAME"
pdflatex -interaction=nonstopmode -halt-on-error "$JOB_NAME.tex"
pdflatex -interaction=nonstopmode -halt-on-error "$JOB_NAME.tex"

warning_pattern='LaTeX Warning|Package .* Warning|Warning--|undefined|Undefined|Error|Citation|Overfull|Underfull|TODO'

if command -v rg >/dev/null 2>&1; then
  if rg -n "$warning_pattern" "$JOB_NAME.log" "$JOB_NAME.blg" "$JOB_NAME.tex" sections references.bib; then
    echo "error: compile finished, but warnings or issues were found" >&2
    exit 1
  fi
else
  if grep -RInE "$warning_pattern" "$JOB_NAME.log" "$JOB_NAME.blg" "$JOB_NAME.tex" sections references.bib; then
    echo "error: compile finished, but warnings or issues were found" >&2
    exit 1
  fi
fi

echo "ok: built $ARTICLE_DIR/$JOB_NAME.pdf"
