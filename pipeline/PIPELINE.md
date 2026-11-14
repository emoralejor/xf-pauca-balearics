# Upstream pipeline (commands)

Raw reads were obtained from ENA/SRA (accessions in data/TableS1_genomes). Only the
processed inputs needed to regenerate the figures are included under data/.

1. Mapping / variant calling: snippy 4.6.0 (BWA-MEM 0.7.18 + FreeBayes), defaults
   (depth >=10x, MQ >=60), reference 9a5c (AE003849.1) or De Donno (CP020870.1).
2. Core alignment: snippy-core.
3. Recombination masking: Gubbins 3.4.3.
4. ML trees: IQ-TREE 2.4.0, ModelFinder (BIC), +ASC on SNP alignments, 1000 UFBoot.
5. Introgression: fastGEAR (MATLAB Runtime R2016a); ClonalFrameML.
6. Dating: BEAST (strict clock, lognormal rate prior median 3.3e-7, constant coalescent).

See Methods for full parameters. [Pega aqui los comandos exactos ejecutados.]
