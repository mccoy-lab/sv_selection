#!/bin/bash
#SBATCH --job-name=mixcr
#SBATCH --time=36:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH--array=0-231%24

ml java
ml samtools
ml htslib

cd /scratch/users/rmccoy22@jhu.edu/mixcr

cram_path=`sed "${SLURM_ARRAY_TASK_ID}q;d" sample_list_paths.txt`
sample_id=`basename ${cram_path} .cram`

samtools fastq \
--reference /work-zfs/rmccoy22/resources/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
-@12 \
-1 ${sample_id}.R1.fq \
-2 ${sample_id}.R2.fq \
${cram_path}

bgzip -@ 12 ${sample_id}.R1.fq
bgzip -@ 12 ${sample_id}.R2.fq

~/code/mixcr-3.0.13/mixcr align \
-f \
-OvParameters.geneFeatureToAlign=VGeneWithP \
-s hsa \
--report ${sample_id}_analysis.report \
-t 12 \
${sample_id}.R1.fq.gz ${sample_id}.R2.fq.gz \
~/scratch/${sample_id}_alignments.vdjca

rm ${sample_id}.R1.fq.gz
rm ${sample_id}.R2.fq.gz
