#!/bin/bash

### Downsample SNPs genotyped in 1000 Genomes, to create a smaller variant set
### for inferring sample admixture proportions with Ohana.
### We downsampled 100,000 variants from chr21 of the 1000 Genomes SNV/indel
### data re-called on GRCh38.

############################################################################

# DATA

# path to VCF of 1000 Genomes small variants - we used the chr21 SNP VCF
# generated for calculating LD between SNPs and SVs
SNP_VCF=chr21_2504samples_chrPrefix.vcf.gz

############################################################################


# add variant IDs to VCF
bcftools annotate -O z \
	-o 1KGP_SNV_indel_chr21_IDs.vcf.gz \
	--set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' \
	$SNP_VCF

# add in genotypes for the NA18498 sample (copying genotypes from NA18499)
# which is missing from the 1000 Genomes hg38 small variant calls
Rscript add_NA18498.R 1KGP_SNV_indel_chr21_IDs.vcf.gz
# add header to output VCF
zcat 1KGP_SNV_indel_chr21_IDs.vcf.gz | grep "##" > header.txt
cat header.txt NA18498.vcf > 1KGP_SNV_indel_chr21_IDs_NA19498.vcf
bgzip 1KGP_SNV_indel_chr21_IDs_NA19498.vcf
tabix 1KGP_SNV_indel_chr21_IDs_NA19498.vcf

# create .map and .ped files for plink
plink --vcf 1KGP_SNV_indel_chr21_IDs_NA19498.vcf.gz \
	--recode \
	--out 1KGP_SNV_indel_chr21_IDs_NA19498.vcf.gz

# prune SNPs
plink --file 1KGP_SNV_indel_chr21_IDs_NA19498.vcf.gz \
	--indep 100 10 1.5 \
	--out chr21

# downsample list of pruned variants to only ~100,000 variants total
awk 'NR % 7 == 0' chr21.prune.in > chr21.prune.in.downsampled.list
# select downsampled variants from VCF
~/code/gatk-4.1.7.0/gatk SelectVariants \
	-O chr21_pruned.vcf \
	-ids chr21.prune.in.downsampled.list \
	-V 1KGP_SNV_indel_chr21_IDs_NA19498.vcf.gz
bgzip chr21_pruned.vcf
rm chr21_pruned.vcf.idx