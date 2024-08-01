#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 0:30:00
#SBATCH -p short
#SBATCH -o make_count_mtx_genes_job_%A_%a.out
#SBATCH --array=1-10
#SBATCH --account=chen_fec176

echo 'indir =' $1
echo 'sample =' $2
echo 'methylation_context =' $3
echo 'reference_gencode =' $4

python ~/methylation/make_count_mtx_genes.py -i $1 -s $2 -m $3 -r $4 -p ${SLURM_ARRAY_TASK_ID}
