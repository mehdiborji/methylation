#!/bin/bash

#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 5:25:00
#SBATCH -p short
#SBATCH -o bsbolt_index_job_%A.out
#SBATCH --account=chen_fec176

ref_fasta=$1
ref_dir=$2

echo $ref_fasta
echo $ref_dir

bsbolt Index -G $ref_fasta -DB $ref_dir


bsbolt Align -DB ../GRCm39_bolt/ -t 12 -OT 8 -F1 r1.fq -F2 r2.fq -O bolt_out

bsbolt Align -DB ../GRCm39_bolt/ -t 12 -OT 8 -F1 r2.fq -F2 r1.fq -O bolt_out