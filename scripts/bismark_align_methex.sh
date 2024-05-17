#!/bin/bash

#SBATCH -c 3
#SBATCH --mem=12G
#SBATCH -t 1:30:00
#SBATCH -p short
#SBATCH -o methylation_extractor_job_%A.out
#SBATCH --account=chen_fec176

read1=$1
read2=$2
out_dir=$3
ref_dir=$4
bam=$5

module load gcc/9.2.0
module load bowtie2/2.5.1

/home/meb521/Bismark-0.24.2/bismark --pbat --score_min L,0,-.2 $ref_dir -1 $read1 -2 $read2 -o $out_dir

/home/meb521/Bismark-0.24.2/bismark_methylation_extractor $bam \
--comprehensive --merge_non_CpG --gzip --buffer_size 1G \
--include_overlap --paired-end --no_header --output_dir $out_dir
