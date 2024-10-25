#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -t 0:30:00
#SBATCH -p short
#SBATCH -o aggregate_quad_parts_CpG_job_%A_%a.out
#SBATCH --array=1-3
#SBATCH --account=chen_fec176

echo 'indir =' $1
echo 'sample =' $2

python ~/methylation/aggregate_quad_parts.py -i $1 -s $2 -m CpG_context -p ${SLURM_ARRAY_TASK_ID}
