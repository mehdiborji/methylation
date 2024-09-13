#!/bin/bash
#SBATCH -c 10
#SBATCH --mem=10G
#SBATCH -t 0:10:00
#SBATCH -p short
#SBATCH -o stack_mtx_genes_job_%A.out
#SBATCH --account=chen_fec176

echo 'indir =' $1
echo 'sample =' $2
echo 'methylation_context =' $3
echo 'reference_gencode =' $4

python ~/methylation/stack_mtx_genes.py -c 20 -i $1 -s $2 -m $3 -r $4
