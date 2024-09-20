#!/bin/bash
#SBATCH -c 20
#SBATCH --mem=8G
#SBATCH -t 1:40:00
#SBATCH -p short
#SBATCH -o tag_bam_parts_job_%A.out
#SBATCH --account=chen_fec176

echo 'indir =' $1
echo 'sample =' $2

python -u ~/methylation/tag_bam_parts.py -c 20 -i $1 -s $2 --genomic
#--limit
