#!/bin/bash
#SBATCH -c 20
#SBATCH --mem=10G
#SBATCH -t 0:40:00
#SBATCH -p short
#SBATCH -o methyl_fastq_job_%A.out

echo 'indir =' $1
echo 'sample =' $2

python ~/methylation/methyl_fastq_pipeline.py -c 20 -i $1 -s $2
#--limit
