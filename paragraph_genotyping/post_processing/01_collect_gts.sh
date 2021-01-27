#!/bin/bash

### Merges all 1000 Genomes samples genotyped by Paragraph into one combined VCF.
### The merging takes place in multiple batches to minimize memory load.

############################################################################

### DATA

# manifest file of all samples genotyped (same as manifest used to run Paragraph)
MANIFEST=1KGP_samplesCovSexPop.txt
# path to directory with Paragraph output directories for all samples
VCFPATH=../paragraph_output/output
# name of final merged VCF
OUTPUT_VCF=eichlerSVs_1KGP_pgGTs.vcf.gz

############################################################################


# get sample IDs from manifest file
cut -f 1 $MANIFEST | \
	cut -f 1 -d "." | \
	cut -f 2 -d "/" \
	> sample_ids.tmp
# split up 2504 samples into 8 batches to merge separately
split -l 313 \
	--suffix-length=1 \
	--additional-suffix=.tmp \
	--numeric-suffixes=1 \
	sample_ids_batch

# merge sample VCFs for each of the 8 batches
for i in {1..8}; do
	while read line; do
                # add path to VCF to a master file
                echo ${VCFPATH}/${line}/genotypes.vcf.gz >> sample_vcfs_batch${i}.tmp
        done < sample_ids_batch${i}.tmp
        echo "Finished indexing sample VCFs for batch "$i

	# merge sample VCFs for that batch
	bcftools merge -O z \
		-o merge_tmp_batch${i}.vcf.gz \
		-l sample_vcfs_batch${i}.tmp \
		--force-samples \
		--threads 48
	tabix merge_tmp_batch${i}.vcf.gz
	echo "Finished merging sample VCFs for batch "$i
done

# merge sets of two batches together
for i in $(seq 1 2 8); do
	(( j=$i+1 )) # second batch VCF
	bcftools merge -O z \
		-o merge_tmp_bigbatch${i}.vcf.gz \
		--threads 48 \
		merge_tmp_batch${i}.vcf.gz merge_tmp_batch${j}.vcf.gz
	echo merge_tmp_bigbatch${i}.vcf.gz >> sample_vcfs_bigbatch.tmp
	echo "Finished merging batches "$i "and "$j
done

# merge all batches into one complete VCF
bcftools merge -O z \
	-o $OUTPUT_VCF \
	--threads 48 \
	-l sample_vcfs_bigbatch.tmp
echo "Finished merging all batches"

# remove tmp files
echo "Removing tmp files..."
rm sample_ids.tmp sample_ids_batch* sample_vcfs_batch* merge_tmp_batch* merge_tmp_bigbatch* sample_vcfs_bigbatch.tmp 

echo "Done!"