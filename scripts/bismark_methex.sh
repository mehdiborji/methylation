#!/bin/bash

#SBATCH -c 3
#SBATCH --mem=1G
#SBATCH -t 0:30:00
#SBATCH -p short
#SBATCH -o methylation_extractor_job_%A.out

bam=$1
outdir=$1

/home/meb521/Bismark-0.24.2/bismark_methylation_extractor $1 \
--comprehensive --merge_non_CpG --gzip --buffer_size 1G \
--include_overlap --paired-end --no_header --output_dir $2
