#!/bin/bash

#SBATCH -c 20
#SBATCH --mem=16G
#SBATCH -t 2:30:00
#SBATCH -p priority
#SBATCH -o align_parts_job_%A_%a.out
#SBATCH --array=1
#SBATCH --account=chen_fec176

echo 'indir =' $1
echo 'sample =' $2
echo 'reference =' $3

python ~/methylation/align_parts_minimap.py -i $1 -s $2 -r $3 -p ${SLURM_ARRAY_TASK_ID}
