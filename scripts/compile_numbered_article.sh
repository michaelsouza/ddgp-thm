#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ARTICLE_DIR="$ROOT_DIR/article"
BASE_JOB="main_numbered_base"
FINAL_JOB="main_numbered"
OVERLAY_TEX="$ARTICLE_DIR/main_numbered_overlay.tex"

require_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "error: required command not found: $1" >&2
    exit 1
  fi
}

require_cmd python3
require_cmd pdflatex
require_cmd bibtex
require_cmd synctex
require_cmd pdfinfo

cd "$ARTICLE_DIR"
rm -f "$BASE_JOB".{aux,bbl,blg,log,out,pdf,synctex.gz} "$FINAL_JOB".{aux,log,out,pdf} "$OVERLAY_TEX"

pdflatex -jobname="$BASE_JOB" -synctex=1 -interaction=nonstopmode -halt-on-error main.tex
bibtex "$BASE_JOB"
pdflatex -jobname="$BASE_JOB" -synctex=1 -interaction=nonstopmode -halt-on-error main.tex
pdflatex -jobname="$BASE_JOB" -synctex=1 -interaction=nonstopmode -halt-on-error main.tex

cd "$ROOT_DIR"
python3 scripts/build_numbered_overlay.py \
  --article-dir "$ARTICLE_DIR" \
  --base-pdf "$BASE_JOB.pdf" \
  --out-tex "$OVERLAY_TEX"

cd "$ARTICLE_DIR"
pdflatex -jobname="$FINAL_JOB" -interaction=nonstopmode -halt-on-error "$(basename "$OVERLAY_TEX")"

echo "ok: built $ARTICLE_DIR/$FINAL_JOB.pdf"
