#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=20G
#SBATCH -t 3:00:00
#SBATCH -p short
#SBATCH -o save_quad_job_%A_%a.out
#SBATCH --array=1-13
#SBATCH --account=chen_fec176

echo 'indir =' $1
echo 'sample =' $2

python ~/methylation/save_quad_batch.py -i $1 -s $2 -p ${SLURM_ARRAY_TASK_ID}
#--limit