#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 2:40:00
#SBATCH -p priority
#SBATCH -o merge_bam_job_%A.out
#SBATCH --account=wu_cjw1


# -f : force overwrite, -n : sort by name
samtools merge -@20 -n -f -o $1/$2/$2_name_sorted.bam $1/$2/split/*_tagged.bam

# -m : add ms (mate score) tags, markdup uses select the best reads
samtools fixmate -@20 -m $1/$2/$2_name_sorted.bam $1/$2/$2_name_sorted_fixmate.bam

# sort by position
samtools sort -@20 -o $1/$2/$2_pos_sorted.bam $1/$2/$2_name_sorted_fixmate.bam
samtools index -@20 $1/$2/$2_pos_sorted.bam

# markdup and use barcode tag in the field BC
samtools markdup -@12 --barcode-tag BC $1/$2/$2_pos_sorted.bam $1/$2/$2_markdup.bam
samtools index -@20 $1/$2/$2_markdup.bam
samtools view -@20 -c $1/$2/$2_markdup.bam
