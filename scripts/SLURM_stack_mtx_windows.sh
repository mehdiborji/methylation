#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=8G
#SBATCH -t 0:20:00
#SBATCH -p priority
#SBATCH -o stack_mtx_windows_job_%A.out
#SBATCH --account=chen_fec176

echo 'indir =' $1
echo 'sample =' $2
echo 'window =' $3
echo 'methylation_context =' $4
echo 'reference_index =' $5

python ~/methylation/stack_mtx_windows.py -c 4 -i $1 -s $2 -w $3 -m $4 -r $5
