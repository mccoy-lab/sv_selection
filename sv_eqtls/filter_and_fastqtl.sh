#!/bin/bash

### Extracts and pre-processes Geuvadis RNA-seq data from recount3.
### Performs cis-eQTL mapping for genotyped SVs using
### a modified version of fastQTL (https://github.com/hall-lab/fastqtl).

############################################################################

### DATA

# VCF of SV genotypes
SV_VCF=eichlerSVs_1KGP_pgGTs_noseq.vcf.gz

# path to bin folder of fastQTL install
FASTQTL_PATH=/fastqtl/bin/

############################################################################


### Data pre-processing

# get covariates for 1000 Genomes samples
./compute_covariates.sh

# both of these files are read in by the `process_geuvadis_data.R` script below
# names of samples in SV VCF
bcftools query -l $SV_VCF \
	> samples.txt
# names of samples included in Geuvadis - TSV with "run accession" and "sample alias"
# from https://www.ebi.ac.uk/ena/browser/view/PRJEB3366
# filereport_read_run_PRJEB3366_tsv.txt

# R script to process Geuvadis RNA-seq data and output phenotype data
Rscript process_geuvadis_data.R
# zip and index Geuvadis phenotype data
sort -k1,1 -k2,2n phenotypeData.bed > phenotypeDataSorted.bed
bgzip phenotypeDataSorted.bed
tabix -p bed phenotypeDataSorted.bed.gz

# subset SV VCF to have the same samples as Geuvadis data
bcftools view -O z \
	-o filteredGenotypeData.vcf.gz \
	-S intersectedSamples.txt \
	$SV_VCF


############################################################################

### Running fastqtl

# Running fastqtl nominal pass with covariates file: 
~/code/fastqtl/bin/fastQTL -V filteredGenotypeData.vcf.gz \
	-B phenotypeDataSorted.bed.gz \
	-O qtls1.out.newcov \
	--cov covariates2.cov \
	--chunk 1 25 \
	--exclude-sites low_freq_hwe_gtrate_sv.exc 
# Running fastqtl permutation pass with covariates file: 
~/code/fastqtl/bin/fastQTL -V filteredGenotypeData.vcf.gz \
	-B phenotypeDataSorted.bed.gz \
	--permute 1000 \
	-O qtls1.out.permcovnew \
	--cov covariates2.cov \
	--chunk 1 25 \
	--exclude-sites low_freq_hwe_gtrate_sv.exc