#!/bin/bash
#SBATCH --job-name=mixcr
#SBATCH --time=48:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH--array=1-231%24

ml java
ml samtools
ml htslib

cd /scratch/users/rmccoy22@jhu.edu/mixcr

cram_path=`sed "${SLURM_ARRAY_TASK_ID}q;d" sample_list_paths.txt`
sample_id=`basename ${cram_path} .cram`

if [ ! -f ${sample_id}.R1.fq.gz ]; then
samtools fastq \
-@ 12 \
-c 6 \
--reference /work-zfs/rmccoy22/resources/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
-1 ${sample_id}.R1.fq.gz \
-2 ${sample_id}.R2.fq.gz \
${cram_path}
fi

~/code/mixcr-3.0.13/mixcr analyze shotgun \
-s hsa \
--starting-material dna \
--report ${sample_id}_analysis.report \
-t 12 \
${sample_id}.R1.fq.gz ${sample_id}.R2.fq.gz \
${sample_id}_mixcr_out

rm ${sample_id}.R1.fq.gz
rm ${sample_id}.R2.fq.gz
