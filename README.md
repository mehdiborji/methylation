


- The first step is to do some preprocessing on the input FASTQ files.
the pipeline assumes the reads are in a directory arrange in the following format

```
/input_directory/sample_name_R1_001.fastq.gz # Read1 of DNA fragment
/input_directory/sample_name_R2_001.fastq.gz # Index (barcode) read
/input_directory/sample_name_R3_001.fastq.gz # Read2 of DNA fragment

```

For example:
```
/n/scratch/users/m/meb521/methyl_seq/nextseq/xBO87_ATAC_S1_R1_001.fastq.gz
/n/scratch/users/m/meb521/methyl_seq/nextseq/xBO87_ATAC_S1_R2_001.fastq.gz
/n/scratch/users/m/meb521/methyl_seq/nextseq/xBO87_ATAC_S1_R3_001.fastq.gz
```

The following script wraps the `methyl_bam_mtx_pipeline.py` into a SLURM job with two input variables required for input:
input directory where the triplet of FASTQ files sit and second 
```
~/methylation/scripts/SLURM_fastq_pipeline.sh /n/scratch/users/m/meb521/methyl_seq/nextseq xBO87_ATAC_S1
```

The pipeline currently is harcoded with the assumption that R2 is 24nt long and has 8nt of splint adapter CAGACGCG at the beginning and reverse compliment of 10x ATAC barcodes from 9-24. It is also harcoded to clip first 15nt and last 2nt of both DNA fragment reads. These options can be modified by modifying `extract_clean_fastq` function within `methyl_utils.py` script

The pipeline does several steps including splitting, trimming, barcode transfer from index reads into cDNA reads, potentially quality filtering reads, and collecting raw barcodes for barcode matching.


- After this step we run a simple script which submits MANY jobs to the HPC for each chunk of fastq


```
python ~/methylation/scripts/methyl_alignment_pipeline.py \
    /n/scratch/users/m/meb521/methyl_seq/nextseq \
    xBO87_ATAC_S1 \
    /n/scratch/users/m/meb521/GRCh38_v44/ \
    job_submit
```

- After alignment postprocessing and count matrix generation is done with the second SLURM pipeline:

```
~/methylation/scripts/SLURM_bam_mtx_pipeline.sh /n/scratch/users/m/meb521/methyl_seq/nextseq xBO87_ATAC_S1 50000
```