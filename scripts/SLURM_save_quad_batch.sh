#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=18G
#SBATCH -t 1:00:00
#SBATCH -p short
#SBATCH -o save_quad_job_%A_%a.out
#SBATCH --account=wu_cjw1

echo 'indir =' $1
echo 'sample =' $2

python ~/methylation/save_quad_batch.py -i $1 -s $2 -p ${SLURM_ARRAY_TASK_ID}
#--limit