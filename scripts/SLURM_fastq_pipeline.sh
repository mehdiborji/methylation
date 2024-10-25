#!/bin/bash
#SBATCH -c 20
#SBATCH --mem=20G
#SBATCH -t 6:30:00
#SBATCH -p short
#SBATCH -o methyl_fastq_job_%A.out
#SBATCH --account=chen_fec176

echo 'indir =' $1
echo 'sample =' $2
echo 'threshold =' $3

python ~/methylation/methyl_fastq_pipeline.py -c 20 -i $1 -s $2 -d $3
#--limit