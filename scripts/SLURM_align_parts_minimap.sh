#!/bin/bash

#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 3:59:00
#SBATCH -p short
#SBATCH -o align_parts_minimap_job_%A_%a.out
#SBATCH --array=2-46
#SBATCH --account=chen_fec176

echo 'indir =' $1
echo 'sample =' $2
echo 'reference =' $3

python ~/methylation/align_parts_minimap.py -i $1 -s $2 -r $3 -p ${SLURM_ARRAY_TASK_ID}
