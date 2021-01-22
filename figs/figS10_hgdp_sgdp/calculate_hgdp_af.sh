#!/bin/bash

### Calculates ref allele frequencies for an HGDP VCF using plink.

# subset VCF to the region of the IGHG4 insertion
tabix -h ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr14.vcf.gz \
	chr14:105422999-105822999 \
	> hgdp.vcf

# add variant IDs to VCF
bcftools annotate -O z \
	-o hgdp_chr14_IDs.vcf.gz \
	--set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' \
	hgdp.vcf

# calculate population-specific allele frequencies with plink
# `within` file is generated at the beginning of the `hgdp_sgdp_AF.R` script
plink --vcf hgdp_chr14_IDs.vcf.gz \
	--freq --keep-allele-order \
	--within hgdp_within.txt \
	--out hgdp_chr14