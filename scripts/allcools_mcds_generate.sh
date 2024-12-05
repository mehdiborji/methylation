#!/bin/bash
#SBATCH -c 20
#SBATCH --mem=48G
#SBATCH -t 10:00:00
#SBATCH -p priority
#SBATCH -o allcools_mcds_generate_job_%A.out
#SBATCH --account=wu_cjw1

allcools generate-dataset \
    --allc_table $1/$2/allc_table.tsv \
    --output_path $1/$2/mcds \
    --chrom_size_path ~/methylation/data/GRCm39_v34_sizes \
    --obs_dim cell --cpu 20 --chunk_size 150 \
    --regions chrom25k 25000 \
    --regions chrom100k 100000 \
    --regions chrom250k 250000 \
    --regions genes ~/methylation/data/gencode.vM35.genes.2k.bed \
    --regions TSS_1k ~/methylation/data/gencode.vM35.TSS.1k.bed \
    --regions TSS_2k ~/methylation/data/gencode.vM35.TSS.2k.bed \
    --quantifiers chrom25k count CG,CH \
    --quantifiers chrom100k count CG,CH \
    --quantifiers chrom250k count CG,CH \
    --quantifiers genes count CG,CH \
    --quantifiers TSS_1k count CG,CH \
    --quantifiers TSS_2k count CG,CH
    
#--regions chrom5k 5000 \
#--quantifiers chrom5k count CG,CH \