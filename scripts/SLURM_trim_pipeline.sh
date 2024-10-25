#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=6G
#SBATCH -t 1:30:00
#SBATCH -p short
#SBATCH -o methyl_trim_job_%A_%a.out
#SBATCH --array=1-6
#SBATCH --account=chen_fec176

echo 'indir =' $1
echo 'sample =' $2

python ~/methylation/methyl_trim_pipeline.py -i $1 -s $2 -p ${SLURM_ARRAY_TASK_ID}
