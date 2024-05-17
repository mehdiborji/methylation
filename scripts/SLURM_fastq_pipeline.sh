#!/bin/bash
#SBATCH -c 12
#SBATCH --mem=6G
#SBATCH -t 0:20:00
#SBATCH -p priority
#SBATCH -o methyl_fastq_job_%A.out
#SBATCH --account=chen_fec176

echo 'indir =' $1
echo 'sample =' $2

python ~/methylation/methyl_fastq_pipeline.py -c 12 -i $1 -s $2
#--limit
