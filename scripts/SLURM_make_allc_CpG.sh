#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -t 0:30:00
#SBATCH -p short
#SBATCH -o make_allc_CpG_job_%A_%a.out
#SBATCH --account=chen_fec176

echo 'indir =' $1
echo 'sample =' $2

python ~/methylation/make_allc.py -i $1 -s $2 -m CpG_context -p ${SLURM_ARRAY_TASK_ID}
