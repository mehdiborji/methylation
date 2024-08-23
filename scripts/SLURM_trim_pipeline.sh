#!/bin/bash
#SBATCH -c 10
#SBATCH --mem=10G
#SBATCH -t 2:20:00
#SBATCH -p priority
#SBATCH -o methyl_trim_job_%A.out
#SBATCH --account=chen_fec176

echo 'indir =' $1
echo 'sample =' $2

python ~/methylation/methyl_trim_pipeline.py -i $1 -s $2
#--limit
