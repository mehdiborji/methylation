# methylranger

This is a computational workflow for bioinformatics analysis of single-cell DNA and DNA methylation sequencing data, 
based on a custom molecular biology modification of 10x Genomics Chromium Single Cell Multiome ATAC + Gene Expression assay.

- The first step is to do some preprocessing on the input FASTQ files.
the pipeline assumes the reads are in a directory arranged in the following format:

```
/input_directory/sample_R1_001.fastq.gz # Read1 of DNA fragment
/input_directory/sample_R2_001.fastq.gz # Index (barcode) read
/input_directory/sample_R3_001.fastq.gz # Read2 of DNA fragment
```

The following script wraps the `methyl_fastq_pipeline.py` into a SLURM job with three input variables required for input
`input_directory` and `sample` and a third variable which is twice the estimate of expected cells and is used for calling:

```
sbatch ~/methylation/scripts/SLURM_fastq_pipeline.sh input_directory sample 30000
```
or
```
cd input_directory
sbatch ~/methylation/scripts/SLURM_fastq_pipeline.sh . sample 30000
```

The pipeline currently is harcoded with the assumption that R2 is 24nt long and has 8nt of splint adapter CAGACGCG at the beginning and reverse compliment of 10x ATAC barcodes from 9-24. These options can be modified by modifying `extract_clean_fastq` function within `methyl_utils.py` script

The pipeline does several steps including splitting, trimming, barcode transfer from index reads into cDNA reads, potentially quality filtering reads, and collecting raw barcodes for barcode matching.

- Adapater trimming is done using fastp

```
sbatch ~/methylation/scripts/SLURM_trim_pipeline.sh . sample
```

- Alignments are done with Bismark and array jobs, one task for each part
```
sbatch ~/methylation/scripts/SLURM_align_parts.sh . sample ../GRCm39_full
sbatch ~/methylation/scripts/SLURM_align_parts.sh . sample ../GRCh38_v44/
```

- Alignments for genomic data is done with minimap2, multiple parts per task:
```
sbatch ~/methylation/scripts/SLURM_align_parts_minimap.sh . sample  ../GRCh38_v44/GRCh38_v44_chrs.mmi
```

- There are partially converted reads which need to be removed
```
sbatch ~/methylation/scripts/SLURM_filter_non_conv.sh . sample
```

- For very large datasets we split the pipeline into pieces and submit them as array jobs:

This is an array job `SLURM_ARRAY_TASK_ID` is embedded as an input argument to `save_quad_batch.py`. 

Each  which subsequently splits the parts into batches of 12 and processess them in pools of two using 2 cores. `TASK_ID` will determine which 12-part batch of parts the task will process. `TASK_ID` 1 will process parts `000` to `011`, `TASK_ID` 4 will process parts `036` to `047` and so on.
Each bam of roughly 8m paired-end reads with fragment average of 100nt will need 18GB RAM to produce jsons

```
sbatch ~/methylation/scripts/SLURM_save_quad_batch.sh . sample
```

Aggregation of batches into single files is done separately for CpG an Non_CpG contexts:

```
sbatch ~/methylation/scripts/SLURM_aggregate_quad_parts_CpG.sh . sample
sbatch ~/methylation/scripts/SLURM_aggregate_quad_parts_Non_CpG.sh . sample
```

- To build count matrices from batches, we first make matrix from each batch and then stack them

First we make count matrix from methylation calls for barcodes in each batch. For each context and window size three different matices are built:
z-scored methylation levels
counts of methylated bases
counts of nonmethylated bases

A mouse example with `GRCm39` reference
```
sbatch ~/methylation/scripts/SLURM_make_count_mtx_windows_CpG.sh . sample 100000 ~/methylation/data/GRCm39_v34_allcontigs.fasta.fai
sbatch ~/methylation/scripts/SLURM_make_count_mtx_windows_Non_CpG.sh . sample 100000 ~/methylation/data/GRCm39_v34_allcontigs.fasta.fai
```

A human example with `GRCh38` reference
```
sbatch ~/methylation/scripts/SLURM_make_count_mtx_windows.sh . sample 100000 CpG_context ~/methylation/data/GRCh38_v44_chrs.fasta.fai 
```

Then we stack all count matricies to make one set of final matrices for each of three in above, we also make coverage matrix, these final matrices are stored in AnnData format.
```
sbatch ~/methylation/scripts/SLURM_stack_mtx_windows.sh . sample 100000 CpG_context ~/methylation/data/GRCm39_v34_allcontigs.fasta.fai
```

- Another way to build count matrices is using gene intervals instead of windows, for this we have a preprocessed version of gencode gtf which is essenetially a bed file with intervals and ensembl gene id. 

For mouse we can obtain `Comprehensive gene annotation` from: `https://www.gencodegenes.org/mouse/release_M35.html`
For human we can obtain one from: `https://www.gencodegenes.org/human/release_46.html`

```
sbatch ~/methylation/scripts/SLURM_make_count_mtx_genes.sh . sample CpG_context ~/methylation/data/gencode.vM35.csv.gz
sbatch ~/methylation/scripts/SLURM_make_count_mtx_genes.sh . sample CpG_context ~/methylation/data/gencode.v46.csv.gz 
```

Then stacking method also needs some modifications
```
sbatch ~/methylation/scripts/SLURM_stack_mtx_genes.sh . sample CpG_context ~/methylation/data/gencode.vM35.csv.gz
```

- To build final bam with duplications marked and barcodes add as tag `BC`
We first add barcode tag into each part and filter out all reads and did not match to whitelist:

```
sbatch ~/methylation/scripts/SLURM_tag_bam_parts.sh . sample
```

After tagging all bam parts we can aggregate and mark duplicates using the following:
```
sbatch ~/methylation/scripts/bam_merge.sh . sample
```
Finally we can compute statistics from the entire bam. We do this in a per chromosome level and in parallel:
```
sbatch ~/methylation/scripts/SLURM_compute_bam.sh . sample ~/methylation/data/GRCm39_v34_allcontigs.fasta.fai
sbatch ~/methylation/scripts/SLURM_compute_bam.sh . sample ~/methylation/data/GRCh38_v44_chrs.fasta.fai
        
```
