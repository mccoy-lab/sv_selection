#!/bin/bash
#SBATCH --partition=shared,parallel

# SAMPLE_COV_LIST=$sampleList
# PARAGRAPH_INPUT=$vcf
# THREADS=$threads
# OUT=$out
# --array=$n-$m%40
# --job-name=para1KGP
# -c $threads
# --time=$minutes

module load python
module load gcc/6.4.0
module load cmake/3.12.1

REFERENCE=GRCh38_full_analysis_set_plus_decoy_hla.fa # reference genome
SCRATCH_DIR=tmp # directory to put tmp files in
PATH_PREFIX_1KGP=1KGP_crams # path to 1KGP crams
PARAGRAPH_SCRIPT=paragraph-static/bin/multigrmpy.py # path to Paragraph `multigrmpy.py` script

export TMPDIR=$SCRATCH_DIR

mkdir -p $OUT/manifests
mkdir -p $OUT/output

i=0
# "line" refers to lines in the sample list
while read -a line; do
    samplePath=$PATH_PREFIX_1KGP/${line[0]}
    cov=${line[1]%.*} # Round down (fine - doesn't need to be precise)
    sex="unknown"
    if [[ ${#line[@]} > 2 ]]; then
	sex=${line[2]}
    fi
    if [[ $i -eq $SLURM_ARRAY_TASK_ID ]]; then
        sample=${samplePath%.final.cram}
        sample=${sample##*/}
        if [[ ! -e $OUT/manifests/manifest_$sample.txt ]]; then
            # Create manifest file
            echo -e 'id\tpath\tdepth\tread length\tsex' > $OUT/manifests/manifest_$sample.txt
            echo -e $sample'\t'$samplePath'\t'$cov'\t150\t'$sex >> $OUT/manifests/manifest_$sample.txt
        fi
        
        # Make output directory and seperate scratch folder (to remove temp files even while other samples run concurrently)
        mkdir -p $OUT/output/${sample}
        mkdir -p $SCRATCH_DIR/$sample

	export TMPDIR=$SCRATCH_DIR/$sample

        echo 'Running '$sample
        ((M_setting= 20 * cov))
        # Run paragraph. -M 20*cov recommended by Sai to help speed up runs (don't try too hard in high depth regions)
        $PARAGRAPH_SCRIPT -m $OUT/manifests/manifest_$sample.txt -i $PARAGRAPH_INPUT -M $M_setting -o $OUT/output/${sample} -r $REFERENCE --threads $THREADS --scratch-dir $SCRATCH_DIR/$sample
	
	echo 'Completed run for '$sample	

        # Remove temp files so we don't fill the scratch directory
        rm -r $SCRATCH_DIR/$sample
	
	echo 'Finished removing tmp files for '$sample
    fi
    ((i=i+1))
done < $SAMPLE_COV_LIST