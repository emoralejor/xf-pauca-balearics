#!/usr/bin/env bash
# Coloca Mag1_15620328 y Mag1_10513216 en la filogenia pauca (Linaje I vs II).
set -uo pipefail
SNP=~/pauca/snippy
OUT=~/pauca/mag/run
mkdir -p "$OUT"

find_bin(){ find ~/miniconda3/envs -maxdepth 3 -name "$1" -type f 2>/dev/null | head -1; }
envof(){ echo "$1" | sed -E 's#.*/envs/([^/]+)/bin/.*#\1#'; }
SC=$(find_bin snippy-core); SS=$(find_bin snp-sites); IQ=$(find_bin iqtree2); [ -z "$IQ" ] && IQ=$(find_bin iqtree)
echo "snippy-core: $SC"; echo "snp-sites: $SS"; echo "iqtree: $IQ"
[ -z "$SC" ] && { echo "no encuentro snippy-core"; exit 1; }

cd "$SNP"
# todas las muestras pauca + Fillmore (outgroup) + los 2 Mag; se excluyen otras subespecies
DIRS=$(ls -d */ 2>/dev/null | sed 's#/##' | grep -vxE 'ATCC35871|Ann1|Riv11|XYL1752_17|XYL1968_18')
echo ">> muestras incluidas: $(echo "$DIRS" | wc -w)"
echo "$DIRS" | grep -i mag || echo "   (ojo: no veo dirs Mag en la lista)"

conda run -n "$(envof "$SC")" snippy-core --ref 9a5c/ref.fa --prefix "$OUT/core_mag" $DIRS
cd "$OUT"

# SNP alignment conservando datos faltantes (N) de los Mag
if [ -n "$SS" ]; then conda run -n "$(envof "$SS")" snp-sites -o core_mag.snps.aln core_mag.full.aln
else echo "no snp-sites; uso core.aln de snippy"; cp core_mag.aln core_mag.snps.aln; fi

# cobertura real de cada Mag (fraccion de posiciones ACGT en el alineamiento completo)
conda run -n rggtree Rscript -e '
a<-ape::read.dna("core_mag.full.aln",format="fasta"); m<-as.character(as.matrix(a))
f<-function(x) mean(x %in% c("a","c","g","t"))
for(s in grep("Mag1",rownames(m),value=TRUE)) cat(sprintf("%-16s coverage %.1f%%\n",s,100*f(m[s,])))'

# arbol rapido (GTR+ASC)
conda run -n "$(envof "$IQ")" "$(basename "$IQ")" -s core_mag.snps.aln -m GTR+ASC -B 1000 -T AUTO --prefix mag_tree -redo >/dev/null 2>&1 || \
conda run -n "$(envof "$IQ")" "$(basename "$IQ")" -s core_mag.snps.aln -m GTR -B 1000 -T AUTO --prefix mag_tree -redo

# clasificacion por vecino mas cercano (independiente del enraizado) + sister clade
conda run -n rggtree Rscript -e '
library(ape); t<-read.tree("mag_tree.treefile")
LI<-c("9a5c","CVC0251","CVC0256","J1a12","XRB","50066","50080","11399","B111","U24D","Hib4","Pr8x","6c","COF0324","RAAR16_Co13","RAAR15_Co33","Fb7","ALM4","AM2-Angelina","OLI17A","32","3124","XF39","XF43","XF64.2","XF71_BR","XF-73-BR","XF143","XF217-2","XF221","XF45-A3")
LII<-c("De_Donno","CoDiRO","APL2","APL9","APL12","APL18","APL31","APL32","APL37","APL42","APL50","APL57","APL64","APL69","APL71","APL72","APL75","APL76","MentOliv","MentPoly","COF0407","OLS0478","OLS0479","CFBP8072","CFBP8498","PD7211","PD7202","LSV47.09","LSV47.33","M5","M7","Salento-1","Salento-2")
D<-cophenetic(t)
for(mg in grep("Mag1",t$tip.label,value=TRUE)){
  d<-D[mg,]; d<-d[names(d)!=mg]; nn<-names(sort(d))[1:6]
  nI<-sum(nn%in%LI); nII<-sum(nn%in%LII)
  cl<-if(nII>nI)"Lineage II (ST53/Central-American)" else if(nI>nII)"Lineage I (South American)" else "ambiguo"
  cat("\n",mg,"\n  6 vecinos mas cercanos:",paste(nn,collapse=", "),"\n  -> ",cl," (",nI," LI /",nII," LII)\n")
}'
echo "arbol en: $OUT/mag_tree.treefile"
