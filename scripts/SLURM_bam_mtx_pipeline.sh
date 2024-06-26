#!/bin/bash
#SBATCH -c 12
#SBATCH --mem=10G
#SBATCH -t 0:10:00
#SBATCH -p short
#SBATCH -o methyl_bam_mtx_job_%A.out
#SBATCH --account=chen_fec176

echo 'indir =' $1
echo 'sample =' $2
echo 'window =' $3
echo 'context =' $4
echo 'reference_index =' $5

python ~/methylation/methyl_bam_mtx_pipeline.py -c 14 -i $1 -s $2 -w $3 -m $4 -r $5 --limit
