#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 5:00:00
#SBATCH -p priority
#SBATCH -o merge_bam_piped_job_%A.out
#SBATCH --account=wu_cjw1

samtools merge -@20 -u - $1/$2/split/*_tagged.bam | \
samtools fixmate -@20 -m - - | \
samtools sort -@20 -u - | \
samtools markdup -@20 --barcode-tag BC - $1/$2/$2_markdup.bam
samtools index -@20 $1/$2/$2_markdup.bam
samtools view -@20 -c $1/$2/$2_markdup.bam