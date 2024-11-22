#!/bin/bash

#SBATCH -c 5
#SBATCH --mem=12G
#SBATCH -t 3:30:00
#SBATCH -p priority
#SBATCH -o align_parts_job_%A_%a.out
#SBATCH --account=chen_fec176

#wu_cjw1
#chen_fec176

echo 'indir =' $1
echo 'sample =' $2
echo 'reference =' $3

python ~/methylation/align_parts.py -i $1 -s $2 -r $3 -p ${SLURM_ARRAY_TASK_ID}
