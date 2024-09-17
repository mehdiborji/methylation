#!/bin/bash

#SBATCH -c 1
#SBATCH --mem=500M
#SBATCH -t 0:40:00
#SBATCH -p short
#SBATCH -o filter_non_conversion_job_%A_%a.out
#SBATCH --array=3-228
#SBATCH --account=chen_fec176

echo 'indir =' $1
echo 'sample =' $2

python ~/methylation/filter_non_conversion.py -i $1 -s $2 -p ${SLURM_ARRAY_TASK_ID}
