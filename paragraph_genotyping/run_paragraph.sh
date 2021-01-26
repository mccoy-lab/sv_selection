### Wrapper script for submitting SLURM array jobs to run Paragraph on all 1000 Genomes samples.

t=false
out="."

while getopts hs:v:p:t:b:o:n option; do
case "${option}"
in
    h) echo "Usage:
       -s list of samples and covs, tab-separated (1KGP_samplesCovSex.txt for all 1KGP)
       -v vcf, paragraph compliant
       -p threads
       -t time in minutes per single job (2 seconds per SV in vcf / num threads per Paragraph docs, but we recommend approx doubling this or more to be safe)
       -b batch size (num jobs for slurm to submit at once; recommended to do ~500 or so not all 2504 at once)
       -o output folder path. Default is '.'
       -n Optional. Include -n for a new run to use test mode and only run a single sample.
       "
       exit 0 ;;
    s) sampleList=${OPTARG};;
    v) vcf=${OPTARG};;
    p) threads=${OPTARG};;
    t) minutes=${OPTARG};;
    b) skip=${OPTARG};;
    o) out=${OPTARG};;
    n) test=true;;
esac
done

if [[ ! $sampleList || ! $vcf || ! $threads || ! $minutes || ! $skip ]]; then
    echo "Missing mandatory argument. -s, -v, -p, -t, -b are all required."
    echo "Usage:                                                    
       -s list of samples and covs, tab-separated (1KGP_samplesCovSex.txt for all 1KGP)
       -v vcf, paragraph compliant
       -p threads
       -t time in minutes per single job (2 seconds per SV in vcf / num threads per Paragraph docs, but we recommend approx doubling this or more to be safe)
       -b batch size (num jobs for slurm to submit at once; recommended to do ~500 or so not all 2504 at once)
       -o output folder path. Default is '.'
       -n Optional. Include -n for a new run to use test mode and only run a single sample.                                                      
       "
    exit 0
fi
    
echo "Sample list:" $sampleList
echo "Input vcf:" $vcf
echo "Threads:" $threads
echo "Time (minutes):" $minutes
echo "Batch Size:" $skip
echo "Output:" $out

if [[ "$test" == true ]]; then
    echo "Test mode; only processing one sample."
fi

((numLines= $(cat $sampleList | wc -l) - 1))  # max index. Num lines -1.

if [[ "$test" == true ]]; then
    ((numLines=0))
    echo "test mode"
fi

for n in `seq 0 $skip $numLines`; do
    (( m= n + skip - 1 ))
    if [[ $m -gt $numLines ]]; then
	m=$numLines
    fi
    
    if [[ $n -eq 0 ]]; then #No dependencies the first time
	sbatch --export=ALL,SAMPLE_COV_LIST=$sampleList,PARAGRAPH_INPUT=$vcf,THREADS=$threads,OUT=$out \
		--array=$n-$m%40 --job-name=para1KGP -c $threads --time=$minutes --mail-user=syan11@jhu.edu --mail-type=end \
			slurmJobParagraph.sh

    else #After that, wait for all to finish, then launch more.
	sbatch --export=ALL,SAMPLE_COV_LIST=$sampleList,PARAGRAPH_INPUT=$vcf,THREADS=$threads,OUT=$out \
		--dependency=singleton --job-name=para1KGP -c $threads --time=$minutes --array=$n-$m%40 --mail-user=syan11@jhu.edu --mail-type=end \
			slurmJobParagraph.sh
    fi
    
done
