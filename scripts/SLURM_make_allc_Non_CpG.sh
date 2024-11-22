#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=30G
#SBATCH -t 2:30:00
#SBATCH -p short
#SBATCH -o make_allc_Non_CpG_job_%A_%a.out
#SBATCH --account=wu_cjw1

echo 'indir =' $1
echo 'sample =' $2

python ~/methylation/make_allc.py -i $1 -s $2 -m Non_CpG_context -p ${SLURM_ARRAY_TASK_ID}
