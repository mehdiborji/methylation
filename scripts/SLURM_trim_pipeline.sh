#!/bin/bash
#SBATCH -c 5
#SBATCH --mem=5G
#SBATCH -t 0:30:00
#SBATCH -p short
#SBATCH -o methyl_trim_job_%A_%a.out
#SBATCH --array=3-23
#SBATCH --account=chen_fec176

echo 'indir =' $1
echo 'sample =' $2

python ~/methylation/methyl_trim_pipeline.py -i $1 -s $2 -p ${SLURM_ARRAY_TASK_ID}
