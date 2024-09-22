#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=40G
#SBATCH -t 3:30:00
#SBATCH -p short
#SBATCH -o make_count_mtx_windows_Non_CpG_job_%A_%a.out
#SBATCH --array=1-27
#SBATCH --account=chen_fec176

echo 'indir =' $1
echo 'sample =' $2
echo 'window =' $3
echo 'reference_index =' $4

python ~/methylation/make_count_mtx_windows.py -i $1 -s $2 -w $3 -m Non_CpG_context -r $4 -p ${SLURM_ARRAY_TASK_ID}
