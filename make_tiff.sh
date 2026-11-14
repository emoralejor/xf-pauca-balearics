#!/usr/bin/env bash
# Convierte los PDF de figures/ a TIFF 600 dpi (LZW) para envio a revista.
set -uo pipefail
cd "$(dirname "$0")"
for f in figures/*.pdf; do
  [ -e "$f" ] || continue; out="${f%.pdf}.tiff"
  if command -v magick >/dev/null 2>&1; then magick -density 600 "$f" -flatten -compress lzw "$out" && echo "   $out"
  elif command -v convert >/dev/null 2>&1; then convert -density 600 "$f" -flatten -compress lzw "$out" && echo "   $out"
  else echo "   instala imagemagick+ghostscript para generar TIFF"; break; fi
done
