#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=4G
#SBATCH -t 1:40:00
#SBATCH -p short
#SBATCH -o allc_agg_sort_job_%A_%a.out
#SBATCH --account=wu_cjw1

echo 'indir =' $1
echo 'sample =' $2

python ~/methylation/allc_agg_sort_tab.py -c 4 -i $1 -s $2 -p ${SLURM_ARRAY_TASK_ID}