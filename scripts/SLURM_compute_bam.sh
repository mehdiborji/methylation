#!/bin/bash
#SBATCH -c 12
#SBATCH --mem=8G
#SBATCH -t 0:30:00
#SBATCH -p priority
#SBATCH -o compute_bam_job_%A.out
#SBATCH --account=chen_fec176

echo 'indir =' $1
echo 'sample =' $2
echo 'reference_index =' $3

python ~/methylation/compute_bam.py -c 20 -i $1 -s $2 -r $3
#--limit
