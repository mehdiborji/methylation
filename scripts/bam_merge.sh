#!/bin/bash
#SBATCH -c 20
#SBATCH --mem=20G
#SBATCH -t 0:40:00
#SBATCH -p priority
#SBATCH -o merge_bam_job_%A.out
#SBATCH --account=chen_fec176


# -f : force overwrite, -n : sort by name
samtools merge -@20 -n -f -o $1/$2/$2_name_sorted.bam $1/$2/split/*_tagged.bam

# -m : add ms (mate score) tags, markdup uses select the best reads
samtools fixmate -@20 -m $1/$2/$2_name_sorted.bam $1/$2/$2_name_sorted_fixmate.bam

# sort by position
samtools sort -@20 -o $1/$2/$2_pos_sorted.bam $1/$2/$2_name_sorted_fixmate.bam
samtools index -@20 $1/$2/$2_pos_sorted.bam

# markdup and use barcode tag in the field BC
samtools markdup -@20 --barcode-tag BC $1/$2/$2_pos_sorted.bam $1/$2/$2_markdup.bam
samtools index -@20 $1/$2/$2_markdup.bam
