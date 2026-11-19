#!/usr/bin/env bash
# Reconstruye todas las figuras desde la raiz del repo.
set -uo pipefail
cd "$(dirname "$0")"
for s in scripts/figure1_pauca_tree.R scripts/figure2_ST53_panel.R \
         scripts/figure3_tree_fastgear.R scripts/figureS1_midpoint.R \
         scripts/figureS2_ST53_tree.R scripts/figure4_beast_dated.R; do
  if [ -f "$s" ]; then echo ">> $s"; conda run -n rggtree Rscript "$s" || echo "   [error] $s"; fi
done
echo "Figuras en figures/"
