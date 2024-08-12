# methylranger

This is a computational workflow for bioinformatics analysis of single-cell DNA methylation sequencing data, based on modification of 10x Genomics Multiome Kit.


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
~/methylation/scripts/SLURM_fastq_pipeline.sh \
        /n/scratch/users/m/meb521/methyl_seq/nextseq \
        xBO87_ATAC_S1
```
Other examples submitted with sbatch
```
sbatch ~/methylation/scripts/SLURM_fastq_pipeline.sh \
        /n/scratch/users/m/meb521/xBO140/fastqs \
        xBO140a_S1
```

```
sbatch ~/methylation/scripts/SLURM_fastq_pipeline.sh \
        /n/scratch/users/m/meb521/xBO140_nova \
        xBO140_novaseq
```

```
sbatch ~/methylation/scripts/SLURM_fastq_pipeline.sh \
        /n/scratch/users/m/meb521/xBO153 \
        xBO153_ATAC_240606_S1
```

```
sbatch ~/methylation/scripts/SLURM_fastq_pipeline.sh /n/scratch/users/m/meb521/A22KHFFLT3_out xBO173
```

The pipeline currently is harcoded with the assumption that R2 is 24nt long and has 8nt of splint adapter CAGACGCG at the beginning and reverse compliment of 10x ATAC barcodes from 9-24. It is also harcoded to clip first 11nt and last 2nt of R1 first 2nt and last 2nt of R3 reads. These options can be modified by modifying `extract_clean_fastq` function within `methyl_utils.py` script

The pipeline does several steps including splitting, trimming, barcode transfer from index reads into cDNA reads, potentially quality filtering reads, and collecting raw barcodes for barcode matching.


- After this step we run a simple script which submits MANY jobs to the HPC for each chunk of fastq


```
python ~/methylation/methyl_alignment_pipeline.py \
        -r /n/scratch/users/m/meb521/GRCm39_full \
        -i /n/scratch/users/m/meb521/xBO140/fastqs \
        -s xBO140a_S1 -b -j
        

python ~/methylation/methyl_alignment_pipeline.py \
        -r /n/scratch/users/m/meb521/GRCm39_full \
        -i /n/scratch/users/m/meb521/xBO140_nova \
        -s xBO140_novaseq -b -j
```

To monitor the state of each aligment job we can look at last line of the log which contains total number of reads processed so far
```
find . -type f -name 'methylation_extractor_job_*' -exec tail -n 1 {} \;
```

- After alignment postprocessing and count matrix generation is done with the second SLURM pipeline:
Required arguments are window_size for binning, context which can be two values `Non_CpG_context` and `CpG_context` and fasta index of the reference used in alignment two such indices are available in data folder human `GRCh38_v44_chrs.fasta` and mouse `GRCm39_v34_allcontigs.fasta.fai`

```
~/methylation/scripts/SLURM_bam_mtx_pipeline.sh \
        /n/scratch/users/m/meb521/xBO140/fastqs \
        xBO140a_S1 \
        100000 \
        Non_CpG_context \
        ~/methylation/data/GRCm39_v34_allcontigs.fasta.fai

sbatch ~/methylation/scripts/SLURM_bam_mtx_pipeline.sh \
        /n/scratch/users/m/meb521/xBO140/fastqs \
        xBO140a_S1 \
        200000 \
        CpG_context \
        ~/methylation/data/GRCm39_v34_allcontigs.fasta.fai
```

```
sbatch ~/methylation/scripts/SLURM_bam_mtx_pipeline.sh \
        /n/scratch/users/m/meb521/xBO140_nova \
        xBO140_novaseq \
        100000 \
        CpG_context \
        ~/methylation/data/GRCm39_v34_allcontigs.fasta.fai
```




- For very large datasets we split the pipeline into pieces and submit them as array jobs:


This is an array job `SLURM_ARRAY_TASK_ID` is embedded as an input argument to `save_quad_batch.py`. Each  which subsequently splits the parts into batches of 12 and processess them in pools of two using 2 cores. `TASK_ID` will determine which 12-part batch of parts the task will process. `TASK_ID` 1 will process parts `000` to `011`, `TASK_ID` 4 will process parts `036` to `047` and so on.

```
sbatch ~/methylation/scripts/SLURM_save_quad_batch.sh \
        /n/scratch/users/m/meb521/xBO140_nova \
        xBO140_novaseq
```

```
sbatch ~/methylation/scripts/SLURM_aggregate_quad_parts.sh \
        /n/scratch/users/m/meb521/xBO140_nova \
        xBO140_novaseq \
        CpG_context
```



- To build count matrices from batches, we first make matrix from each batch and then stack them

First we make count matrix from methylation calls for barcodes in each batch. For each context and window size three different matices are built:
z-scored methylation levels
counts of methylated bases
counts of nonmethylated bases
```
sbatch ~/methylation/scripts/SLURM_make_count_mtx_windows.sh \
        /n/scratch/users/m/meb521/xBO140_nova \
        xBO140_novaseq \
        100000 \
        CpG_context \
        ~/methylation/data/GRCm39_v34_allcontigs.fasta.fai
```

Then we stack all count matricies to make one set of final matrices for each of three in above, we also make coverage matrix, these final matrices are stored in AnnData format.
```
sbatch ~/methylation/scripts/SLURM_stack_mtx_windows.sh \
        /n/scratch/users/m/meb521/xBO140_nova \
        xBO140_novaseq \
        100000 \
        CpG_context \
        ~/methylation/data/GRCm39_v34_allcontigs.fasta.fai
```

- Another way to build count matrices is using gene intervals instead of windows, for this we have a preprocessed version of gencode gtf which is essenetially a bed file with intervals and ensembl gene id form `https://www.gencodegenes.org/mouse/release_M35.html`

```
sbatch ~/methylation/scripts/SLURM_make_count_mtx_genes.sh \
        /n/scratch/users/m/meb521/xBO140_nova \
        xBO140_novaseq \
        CpG_context \
        ~/methylation/data/gencode.vM35.csv.gz
```
Then stacking method also needs some modifications
```
sbatch ~/methylation/scripts/SLURM_stack_mtx_genes.sh \
        /n/scratch/users/m/meb521/xBO140_nova \
        xBO140_novaseq \
        CpG_context \
        ~/methylation/data/gencode.vM35.csv.gz
```


- To build final bam with duplications and barcodes marked
We first add barcode tag into each part and filter out all reads and did not match to whitelist:
```
sbatch ~/methylation/scripts/SLURM_tag_bam_parts.sh \
        /n/scratch/users/m/meb521/xBO140_nova \
        xBO140_novaseq
```

After tagging all bam parts we can aggregate and mark duplicates using the following:
```
sbatch ~/methylation/scripts/bam_merge.sh \
        /n/scratch/users/m/meb521/xBO140_nova \
        xBO140_novaseq
```
Finally we can compute statistics from the entire bam. We do this in a per chromosome level and in parallel:
```
sbatch ~/methylation/scripts/SLURM_compute_bam.sh \
        /n/scratch/users/m/meb521/xBO140_nova \
        xBO140_novaseq
        ~/methylation/data/GRCm39_v34_allcontigs.fasta.fai
        
```