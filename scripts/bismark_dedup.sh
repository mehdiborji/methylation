#!/bin/bash

#SBATCH -c 4
#SBATCH --mem=20G
#SBATCH -t 0:50:00
#SBATCH -p short
#SBATCH -o bismark_dedup_%A.out
#SBATCH --account=chen_fec176

/home/meb521/Bismark-0.24.2/deduplicate_bismark $1 --paired
