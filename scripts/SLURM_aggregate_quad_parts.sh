#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=40G
#SBATCH -t 1:10:00
#SBATCH -p short
#SBATCH -o aggregate_quad_parts_job_%A_%a.out
#SBATCH --array=1-20
#SBATCH --account=chen_fec176

echo 'indir =' $1
echo 'sample =' $2
echo 'context =' $3

python ~/methylation/aggregate_quad_parts.py -i $1 -s $2 -m $3 -p ${SLURM_ARRAY_TASK_ID}
