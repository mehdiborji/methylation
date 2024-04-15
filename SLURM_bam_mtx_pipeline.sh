#!/bin/bash
#SBATCH -c 10
#SBATCH --mem=64G
#SBATCH -t 0:40:00
#SBATCH -p short
#SBATCH -o methyl_bam_mtx_job_%A.out

echo 'indir =' $1
echo 'sample =' $2
echo 'window =' $3

python ~/methylation/methyl_bam_mtx_pipeline.py -c 10 -i $1 -s $2 -w $3
#--limit
