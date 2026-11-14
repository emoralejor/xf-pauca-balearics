# Xylella fastidiosa subsp. pauca in the Balearic Islands — figures & code

Code and processed data to reproduce the figures of the manuscript
*"Two independent introductions of Xylella fastidiosa subsp. pauca into the
Balearic Islands revealed by whole-genome sequencing"* (Moralejo et al.).

## Reproduce
```bash
conda env create -f environment.yml      # creates env 'rggtree' (or reuse yours)
bash run_all.sh                          # regenerates every figure into figures/
bash make_tiff.sh                        # optional: 600-dpi TIFF for submission
```
All scripts use paths relative to the repository root; run from the repo root.

## Figure → script
| Figure | Script | Main input(s) |
|--------|--------|---------------|
| Fig 1  | scripts/figure1_pauca_tree.R        | pauca69_ML.treefile, meta.tsv, map_coffee_americas.png |
| Fig 2  | scripts/figure2_ST53_panel.R        | core78.filtered_polymorphic_sites.fasta |
| Fig 3  | scripts/figure3_introgression_tree.R | bloques32.treefile |
| Fig 4  | scripts/figure4_beast_dated.R       | apulia_2013_mcc.tree |
| Fig S1 | scripts/figureS1_midpoint.R         | pauca69_ML.treefile, meta.tsv |
| Fig S2 | scripts/figureS2_ST53_tree.R        | core78_ML.treefile |
| Fig S3 | scripts/figureS3_mag_placement.R | mag_tree.treefile |

## Data
`data/` holds the small processed inputs (masked SNP alignments, trees, metadata,
the BEAST MCC tree, the coffee-regions map). Raw sequencing reads are **not**
included; their ENA/SRA accessions are in `data/TableS1_genomes`. `pipeline/PIPELINE.md`
documents the upstream commands.

## Environment
`environment.yml` and `session_info.txt` record the exact software versions used
(R, ggtree, treeio, ape, pheatmap, patchwork, ggplotify, maps, cowplot, phangorn).

## Citation / archive
- Manuscript: [DOI pending]
- Archived release: Zenodo DOI [10.5281/zenodo.XXXXXXX]
- Repository: https://github.com/<user>/xf-pauca-balearics

## License
Code under the MIT License (see LICENSE). Data files are from the cited public
accessions and the authors' own sequencing.
