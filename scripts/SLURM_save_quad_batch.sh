#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=20G
#SBATCH -t 0:60:00
#SBATCH -p short
#SBATCH -o save_quad_job_%A_%a.out
#SBATCH --array=1-228
#SBATCH --account=chen_fec176

echo 'indir =' $1
echo 'sample =' $2

python ~/methylation/save_quad_batch.py -c 1 -i $1 -s $2 -p ${SLURM_ARRAY_TASK_ID}
# --limit
