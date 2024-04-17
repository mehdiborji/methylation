#!/bin/bash
#SBATCH -c 20
#SBATCH --mem=20G
#SBATCH -t 0:10:00
#SBATCH -p short
#SBATCH -o methyl_bam_mtx_job_%A.out

echo 'indir =' $1
echo 'sample =' $2
echo 'window =' $3

python ~/methylation/methyl_bam_mtx_pipeline.py -c 20 -i $1 -s $2 -w $3
#--limit
