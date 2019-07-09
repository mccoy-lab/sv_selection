#!/bin/bash

# computes pairwise FST between populations
# usage: ./01-fst.sh population_1 population_2 input_vcf output_fst
# usage example: ./01-fst.sh WestEurasia Africa ../vcf/EEE_SV-Pop_1.ALL.genotypes.20181204.vcf.gz WestEurasia_Africa.fst 

target_pop=$1
background_pop=$2
input_file=$3
output_file=$4

target_indices=`cat sample_index_manualassign.txt | grep ${target_pop} | cut -f1 | tr '\n' ','`
background_indices=`cat sample_index_manualassign.txt | grep ${background_pop} | cut -f1 | tr '\n' ','`

~/work/progs/vcflib/bin/wcFst \
  --target ${target_indices} \
  --background ${background_indices} \
  --type GL \
  --file ${input_file} \
  > ${output_file}
