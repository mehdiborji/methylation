#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=40G
#SBATCH -t 3:30:00
#SBATCH -p short
#SBATCH -o make_allc_Non_CpG_job_%A_%a.out
#SBATCH --array=1-27
#SBATCH --account=chen_fec176

echo 'indir =' $1
echo 'sample =' $2

python ~/methylation/make_allc.py -i $1 -s $2 -m Non_CpG_context -p ${SLURM_ARRAY_TASK_ID}
