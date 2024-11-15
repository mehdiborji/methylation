#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=10G
#SBATCH -t 0:30:00
#SBATCH -p priority
#SBATCH -o allc_agg_sort_job_%A_%a.out
#SBATCH --account=chen_fec176

echo 'indir =' $1
echo 'sample =' $2

python ~/methylation/allc_agg_sort_tab.py -c 20 -i $1 -s $2