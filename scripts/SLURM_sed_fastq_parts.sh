#!/bin/bash
#SBATCH -c 20
#SBATCH --mem=1G
#SBATCH -t 0:60:00
#SBATCH -p priority
#SBATCH -o sed_fastq_job_%A.out
#SBATCH --account=chen_fec176

echo 'indir =' $1
echo 'sample =' $2

python ~/methylation/sed_fastq_parts.py -c 20 -i $1 -s $2
