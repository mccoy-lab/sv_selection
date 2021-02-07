#!/bin/bash

### Runs Ohana to identify SVs under selection in 1000 Genomes.
### 1. Using downsampled SNPs, infer admixture with `k` ancestry components
### and a covariance matrix describing the relationship between the samples.
### 2. Using the inferred population structure, scan for SVs under selection.
### Ohana does this by looking for variants whose AFs are better explained
### by allowing variation in one ancestry component.

############################################################################

### DATA

# VCF with downsampled variants for admixture inference
DOWNSAMPLED_VCF=chr21_pruned.vcf.gz

# path to Ohana's `bin` folder
OHANAPATH=/progs/ohana/bin/

# VCF of SV genotypes
SV_VCF=eichlerSVs_1KGP_pgGTs_noseq.vcf.gz

############################################################################

### Infer population structure on downsampled variants

# convert downsampled VCF into .ped file with recoded genotypes for Ohana
plink --vcf $DOWNSAMPLED_VCF \
	--out $DOWNSAMPLED_VCF \
	--recode 12 \
	--geno 0.0 \
	--tab

# convert .ped file into .dgm file ("G matrix") for qpas
$OHANAPATH/convert ped2dgm chr21_pruned.vcf.gz.ped \
	chr21_v2_pruned.dgm

# determine dataset population structure on downsampled variants with qpas
# k: number of ancestry components
# mi: max # of iterations
$OHANAPATH/qpas chr21_v2_pruned.dgm \
	-k 8 \
	-qo chr21_v2_pruned_50_Q.matrix \
	-fo chr21_v2_pruned_50_F.matrix \
	-mi 50
echo "Finished running qpas on downsampled variants"

# infer covariance matrix of AFs from Q matrix
$OHANAPATH/nemeco chr21_v2_pruned.dgm \
	chr21_v2_pruned_50_F.matrix \
	-co chr21_v2_pruned_50_C.matrix \
	--max-time 1800
echo "Finished inferring covariance matrix for downsampled variants"


############################################################################

### Scan for SVs under selection in each ancestry component,
### using population structure determined with downsampled variants

# reorder samples in SV VCF to match order in downsampled VCF
bcftools query -l $DOWNSAMPLED_VCF > samples.txt
bcftools view -O z \
	-o eichlerSVs_resort.vcf.gz \
	-S samples.txt \
	$SV_VCF

# convert SV VCF into .ped file with recoded genotypes for Ohana
plink --vcf eichlerSVs_resort.vcf.gz \
	--out eichlerSVs_resort.vcf.gz \
	--recode 12 \
	--tab

# convert .ped file into .dgm file (a.k.a. G matrix) for qpas
$OHANAPATH/convert ped2dgm eichlerSVs_resort.vcf.gz.ped \
	SVs_all.dgm

# produce admixture-corrected AFs for SVs using the downsampled Q matrix
$OHANAPATH/qpas SVs_all.dgm \
	-k 8 \
	-qi chr21_v2_pruned_50_Q.matrix \
	-fo SVs_all_50_F.matrix \
	-e 0.0001 \
	-fq \
	-mi 50
echo "Finished running qpas on all variants"

# generate covariance matrices for each ancestry component
mkdir -p k8/c_matrices
Rscript make_c_matrix.R \
	chr21_v2_pruned_50_C.matrix \
	8

# test for SVs whose AFs are better explained by a selection covariance matrix
# than the genome-wide "neutral" covariance matrix
mkdir -p k8/selscan
for i in {1..8}
do
	$OHANAPATH/selscan SVs_all.dgm \
		SVs_all_50_F.matrix \
		chr21_pruned_50_C.matrix \
		-cs k8/c_matrices/chr21_pruned_50_C_p${i}.matrix \
		> k8/selscan/selscan_50_k8_p${i}.out
	echo "Finished running selscan on ancestry component "$i
done