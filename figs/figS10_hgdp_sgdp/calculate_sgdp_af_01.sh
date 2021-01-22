#!/bin/bash

### Gets subsetted VCF from SGDP for calculating allele frequencies.

# subset VCF to the region of the IGHG4 insertion
tabix -h ~/data/resources/sgdp-263-hs38DH/263.all.vcf.gz \
	chr14:105500000-106500000 \
	> sgdp.vcf

# add variant IDs to VCF
bcftools annotate -O z \
	-o sgdp_chr14_IDs.vcf.gz \
	--set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' \
	sgdp.vcf

# calculate population-specific allele frequencies with python `cyvcf2` package
./calculate_sgdp_af_02.py