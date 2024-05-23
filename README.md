


- The first step is to do some preprocessing on the input FASTQ files.
the pipeline assumes the reads are in a directory arranged in the following format:

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

The following script wraps the `methyl_fastq_pipeline.py` into a SLURM job with two input variables required for input:
`input_directory` and `sample_name`
```
~/methylation/scripts/SLURM_fastq_pipeline.sh /n/scratch/users/m/meb521/methyl_seq/nextseq xBO87_ATAC_S1
```
Another example submitted with sbatch
```
sbatch ~/methylation/scripts/SLURM_fastq_pipeline.sh /n/scratch/users/m/meb521/xBO140/fastqs xBO140a_S1
```

The pipeline currently is harcoded with the assumption that R2 is 24nt long and has 8nt of splint adapter CAGACGCG at the beginning and reverse compliment of 10x ATAC barcodes from 9-24. It is also harcoded to clip first 11nt and last 2nt of R1 first 2nt and last 2nt of R3 reads. These options can be modified by modifying `extract_clean_fastq` function within `methyl_utils.py` script

The pipeline does several steps including splitting, trimming, barcode transfer from index reads into cDNA reads, potentially quality filtering reads, and collecting raw barcodes for barcode matching.


- After this step we run a simple script which submits MANY jobs to the HPC for each chunk of fastq


```
python ~/methylation/methyl_alignment_pipeline.py \
        -r /n/scratch/users/m/meb521/GRCm39_full \
        -i /n/scratch/users/m/meb521/xBO140/fastqs \
        -s xBO140a_S1 -b -j
```

To monitor the state of each aligment job we can look at last line of the log which contains total number of reads processed so far
```
find . -type f -name 'methylation_extractor_job_*' -exec tail -n 1 {} \;
```

- After alignment postprocessing and count matrix generation is done with the second SLURM pipeline:
Required arguments are window_size for binning, context which can be two values `Non_CpG_context` and `CpG_context` and fasta index of the reference used in alignment two such indices are available in data folder human `GRCh38_v44_chrs.fasta` and mouse `GRCm39_v34_allcontigs.fasta.fai`

```
~/methylation/scripts/SLURM_bam_mtx_pipeline.sh /n/scratch/users/m/meb521/methyl_seq/nextseq xBO87_ATAC_S1 50000 

~/methylation/scripts/SLURM_bam_mtx_pipeline.sh /n/scratch/users/m/meb521/xBO140/fastqs xBO140a_S1 100000 Non_CpG_context ~/methylation/data/GRCm39_v34_allcontigs.fasta.fai

sbatch ~/methylation/scripts/SLURM_bam_mtx_pipeline.sh \
        /n/scratch/users/m/meb521/xBO140/fastqs \
        xBO140a_S1 \
        200000 \
        CpG_context \
        ~/methylation/data/GRCm39_v34_allcontigs.fasta.fai
```