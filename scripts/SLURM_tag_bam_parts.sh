#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=4G
#SBATCH -t 1:40:00
#SBATCH -p priority
#SBATCH -o tag_bam_parts_job_%A.out
#SBATCH --account=wu_cjw1

echo 'indir =' $1
echo 'sample =' $2

python -u ~/methylation/tag_bam_parts.py -c 16 -i $1 -s $2 --genomic
#--limit
#
#
