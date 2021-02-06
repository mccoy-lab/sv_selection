<!-- https://zenodo.org/record/4469976#.YBCgLXdKhs8 -->

## Code associated with "Local adaptation and archaic introgression shape global diversity at human structural variant loci"

In our paper, we used a graph-based approach to genotype long-read discovered SVs in diverse individuals from the 1000 Genomes Project. We searched for locally adaptive SVs and identified strong recent selection at the immunoglobulin locus in Southeast Asian populations. Further investigation revealed that the haplotype under positive selection was adaptively introgressed from Neanderthals.

https://www.biorxiv.org/content/10.1101/2021.01.26.428314v2

### Authors: Stephanie M. Yan, Rachel M. Sherman, Dylan J. Taylor, Divya R. Nair, Andrew N. Bortvin, Michael C. Schatz, Rajiv C. McCoy

Code for genotyping structural variants with [Paragraph](https://github.com/Illumina/paragraph) is available in `paragraph_genotyping`.

Code for metanalysis of SV genotypes (filtering based on genotyping rate and Hardy-Weinberg equilibrium, allele frequency calculation) is in `sv_postprocessing`.

Code for identifying putative adaptive SVs with [Ohana](https://github.com/jade-cheng/ohana) is in `ohana_adaptive_svs`.

Analysis and code for reproducing main and supplementary figures is available in `figs`.
