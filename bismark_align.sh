#!/bin/bash

#SBATCH -c 4
#SBATCH --mem=20G
#SBATCH -t 0:50:00
#SBATCH -p short
#SBATCH -o bismark_align_job_%A.out

read1=$1
read2=$2
out_dir=$3
ref_dir=$4

module load gcc/9.2.0
module load bowtie2/2.5.1

/home/meb521/Bismark-0.24.2/bismark --pbat $4 -1 $read1 -2 $read2 -o $out_dir
