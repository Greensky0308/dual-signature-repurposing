#!/bin/bash
# Fetch the validation-cohort inputs.
#
# These are plain HTTP files on the GEO mirror, so they are scripted. Everything
# else is placed by hand in $HCC_HF_RAW -- see "Inputs" in README.md: the L1000
# release files and the discovery GEO matrices are too large or too indirect to
# script, and the pipeline reads them straight from that directory.
#
# Files land in $HCC_HF_VAL (default: ./downloads), which is where config.py looks.
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST="${HCC_HF_VAL:-$HERE/downloads}"
PROXY="${HCC_HF_PROXY:-}"
mkdir -p "$DEST"
cd "$DEST" || exit 1

fetch() {
  local url="$1" out="$2"
  [ -s "$out" ] && { echo "[skip] $out already present"; return 0; }
  if curl -sL ${PROXY:+-x "$PROXY"} --max-time 1800 -o "$out.part" "$url" \
     && [ -s "$out.part" ] && ! head -c 200 "$out.part" | grep -qi '<html'; then
    mv "$out.part" "$out"; echo "[ok] $out $(du -h "$out" | cut -f1)"; return 0
  fi
  echo "[FAIL] $url -- partial left at $out.part, re-run to retry"; return 1
}

S="https://ftp.ncbi.nlm.nih.gov/geo/series"
P="https://ftp.ncbi.nlm.nih.gov/geo/platforms"

# HCC validation cohort (Illumina HumanHT-12 v4)
fetch "$S/GSE76nnn/GSE76427/matrix/GSE76427_series_matrix.txt.gz" GSE76427_series_matrix.txt.gz
fetch "$P/GPL10nnn/GPL10558/annot/GPL10558.annot.gz" GPL10558.annot.gz

# Heart-failure validation cohort (GSE141910): per-sample RNA-seq matrices
# inside a supplementary tar; the RPKM table supplies its ENSG annotation
fetch "$S/GSE116nnn/GSE116250/suppl/GSE116250_rpkm.txt.gz" GSE116250_rpkm.txt.gz   # ENSG annotation, used by GSE141910
fetch "$S/GSE141nnn/GSE141910/matrix/GSE141910_series_matrix.txt.gz" GSE141910_series_matrix.txt.gz
fetch "$S/GSE141nnn/GSE141910/suppl/GSE141910_RAW.tar" GSE141910_RAW.tar

# the tar holds one gzipped matrix per sample
if [ -s GSE141910_RAW.tar ] && [ ! -d GSE141910_raw ]; then
  mkdir -p GSE141910_raw && tar -xf GSE141910_RAW.tar -C GSE141910_raw && echo "[ok] extracted GSE141910_raw/"
fi

echo "=== done ==="
ls -la
