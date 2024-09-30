#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=40G
#SBATCH -t 1:50:00
#SBATCH -p short
#SBATCH -o aggregate_quad_parts_Non_CpG_job_%A_%a.out
#SBATCH --array=1-6
#SBATCH --account=chen_fec176

echo 'indir =' $1
echo 'sample =' $2

python ~/methylation/aggregate_quad_parts.py -i $1 -s $2 -m Non_CpG_context -p ${SLURM_ARRAY_TASK_ID}
