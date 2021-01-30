#!/bin/bash

### Calculates allele frequency of SVs in each population using PLINK.

############################################################################


# VCF output from Paragraph; no full SV sequences
INPUT_VCF=eichlerSVs_1KGP_pgGTs_noseq.vcf.gz
# file specifying groups (i.e. populations) to calculate AF in
WITHIN_FILE=1KGP_within.txt

plink --vcf $INPUT_VCF \
	--within $WITHIN_FILE \
	--freq gz \
	--keep-allele-order \
	--out eichlerSVs_af_allpops_unfolded