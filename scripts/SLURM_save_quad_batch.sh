#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=40G
#SBATCH -t 3:00:00
#SBATCH -p short
#SBATCH -o methyl_bam_mtx_job_%A_%a.out
#SBATCH --array=5-19
#SBATCH --account=chen_fec176

echo 'indir =' $1
echo 'sample =' $2

python ~/methylation/save_quad_batch.py -c 2 -i $1 -s $2 -p ${SLURM_ARRAY_TASK_ID}
#--limit
