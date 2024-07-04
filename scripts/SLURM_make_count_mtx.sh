#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=48G
#SBATCH -t 1:30:00
#SBATCH -p short
#SBATCH -o make_count_mtx_job_%A_%a.out
#SBATCH --array=1-10
#SBATCH --account=chen_fec176

echo 'indir =' $1
echo 'sample =' $2
echo 'window =' $3
echo 'context =' $4
echo 'reference_index =' $5

python ~/methylation/make_count_mtx.py -i $1 -s $2 -w $3 -m $4 -r $5 -p ${SLURM_ARRAY_TASK_ID}
