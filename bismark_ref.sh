#!/bin/bash

#SBATCH -c 20
#SBATCH --mem=20G
#SBATCH -t 1:25:00
#SBATCH -p priority
#SBATCH -o bismark_genome_preparation_job_%A.out

ref_dir=$1

echo $ref_dir

module load gcc/9.2.0
module load bowtie2/2.5.1

/home/meb521/Bismark-0.24.2/bismark_genome_preparation --verbose $1 --parallel 10
