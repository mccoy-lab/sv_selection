#!/bin/bash

### Finds the genotyping callrate for each SV in the input VCF.
### Calculates Hardy-Weinberg p-values for samples in each population
### and outputs the # of populations in which an SV violates one-sided HWE.

############################################################################

### DATA

# final VCF with all 1KGP genotypes
INPUT_VCF=eicherSVs_1KGP_pgGTs_noseq.vcf.gz
# manifest file of all samples genotyped (same as manifest used to run Paragraph)
MANIFEST=1KGP_samplesCovSexPop.txt

############################################################################


### genotyping callrate
# find the genotyping callrate for each SV. Callrates for SVs on the Y chromosome are
# calculated only for male samples. Writes output to `genotyping_callrate.txt`
./find_gt_callrate.py $INPUT_VCF $MANIFEST

### one-sided Hardy-Weinberg
mkdir -p gt_counts
# create genotype counts files for each population in the `gt_counts` directory, with HWE results
./hwe_count_gts.py $INPUT_VCF $MANIFEST
echo "Finished generating genotype counts files"

# find the # of populations in which an SV violates HWE, write output to `HWE_pvals.txt`
# SVs on the Y chromosome are ignored because they are haploid
./03_HWE_one_sided.R
echo "Finished calculating HWE p-values for SVs"