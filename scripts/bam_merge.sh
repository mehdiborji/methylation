#!/bin/bash
#SBATCH -c 20
#SBATCH --mem=10G
#SBATCH -t 0:20:00
#SBATCH -p priority
#SBATCH -o merge_bam_job_%A.out
#SBATCH --account=chen_fec176

#samtools merge -@20 -n -f -o $1/$2/$2_name_sorted.bam $1/$2/split/output_*/*_tagged.bam
#samtools fixmate -@20 -m $1/$2/$2_name_sorted.bam $1/$2/$2_name_sorted_fixmate.bam
#samtools sort -@20 -o $1/$2/$2_pos_sorted.bam $1/$2/$2_name_sorted_fixmate.bam
#samtools index -@20 $1/$2/$2_pos_sorted.bam

samtools markdup -@20 --json --barcode-tag BC $1/$2/$2_pos_sorted.bam $1/$2/$2_markdup.bam

#samtools sort -@20 -o $1/$2_pos_sorted.deduplicated.bam $1/$2_name_sorted.deduplicated.bam
#samtools index -@20 $1/$2_pos_sorted.deduplicated.bam

#samtools index -@16 $1.bam
#samtools view -@16 -f 4 -b -h $1.bam > $1_unmapped.bam
#samtools index -@16 $1_unmapped.bam
#samtools bam2fq -@16 $1_unmapped.bam > $1_unmapped.fastq
#pigz $1_unmapped.fastq

