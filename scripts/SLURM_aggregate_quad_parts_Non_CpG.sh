#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=30G
#SBATCH -t 0:40:00
#SBATCH -p short
#SBATCH -o aggregate_quad_parts_Non_CpG_job_%A_%a.out
#SBATCH --account=chen_fec176

echo 'indir =' $1
echo 'sample =' $2

python ~/methylation/aggregate_quad_parts.py -i $1 -s $2 -m Non_CpG_context -p ${SLURM_ARRAY_TASK_ID}
