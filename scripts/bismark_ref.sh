#!/bin/bash

#SBATCH -c 20
#SBATCH --mem=20G
#SBATCH -t 1:25:00
#SBATCH -p priority
#SBATCH -o bismark_genome_preparation_job_%A.out
#SBATCH --account=chen_fec176

ref_dir=$1

echo $ref_dir

/home/meb521/Bismark-0.24.2/bismark_genome_preparation --verbose $ref_dir --parallel 10
