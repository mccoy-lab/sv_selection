#!/bin/bash

### Replaces full SV sequences in VCF with symbolic representations (<INS> or <DEL>).
### (Symbolic representations are required by PLINK, which has a character limits for
### ref and alt alleles in VCF.)
### Replaces the `FORMAT/GT` field of Y chromosome variants with "." (no genotype
### call) for all female samples, because this is not automatically done by Paragraph.

############################################################################

# DATA

# manifest file of all samples genotyped (same as manifest used to run Paragraph)
MANIFEST=../1KGP_samplesCovSexPop.txt
# input: merged VCF
INPUT_VCF=eichlerSVs_1KGP_pgGTs.vcf.gz
# name for output VCF without SV sequences
OUTPUT_VCF=eichlerSVs_1KGP_pgGTs_noseq.vcf

############################################################################


# create a modified VCF called noseq.vcf, with <INS> and <DEL> instead of SV sequences
# and no genotype call for Y chromosome variants in female samples
./make_symbolic_alleles.R $INPUT_VCF $MANIFEST
echo "Finished removing SV sequences from VCF"

# get original vcf header
bcftools view -h $INPUT_VCF > header.txt
# add header to blank.vcf
cat header.txt noseq.vcf > $OUTPUT_VCF
bgzip $OUTPUT_VCF
tabix ${OUTPUT_VCF}.gz
echo "Finished creating new VCF"

# remove tmp files
rm noseq.vcf header.txt
echo "Finished removing tmp files"